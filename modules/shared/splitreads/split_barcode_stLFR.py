#!/usr/bin/env python3
"""Parallel stLFR barcode splitter.

Python replacement for the legacy stLFR/cLFR split Perl scripts.  The CLI,
output file names, and summary logs are kept compatible so existing Snakemake
rules and downstream parsers continue to work.
"""

import argparse
import gzip
import multiprocessing as mp
import os
import sys
from collections import defaultdict


BARCODE_HASH = {}
BARCODE_HASH_STRING = {}
OPTIONS = None


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--barcode")
    parser.add_argument("--r1", required=True)
    parser.add_argument("--r2", required=True)
    parser.add_argument("--read_len", type=int, required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--bc_start", type=int, required=True)
    parser.add_argument("--gdna_start", type=int, required=True)
    parser.add_argument("--additional_bc_start", type=int, required=True)
    parser.add_argument("--additional_bc_len", type=int, required=True)
    parser.add_argument("--gdna_start_r1", type=int, required=True)
    parser.add_argument("--read_len_r1", type=int, required=True)
    parser.add_argument("--adapter_len", type=int, default=0)
    parser.add_argument("--cbc_len", type=int, default=0)
    parser.add_argument("--bc_len_redundant", type=int, default=0)
    parser.add_argument("--additional_bc_len_r1", type=int, default=0)
    parser.add_argument("--reverse_complement", action="store_true")
    parser.add_argument("--swap", choices=("none", "ac", "gc"), default="none")
    parser.add_argument(
        "--output_mode",
        choices=("single", "separate", "stratified", "cwgs"),
        default="single",
    )
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument(
        "--chunk_reads",
        type=int,
        default=10000,
        help="Read pairs per worker chunk [default: 10000].",
    )
    args = parser.parse_args()
    if not args.barcode and not args.cbc_len:
        parser.error("Either --barcode (stLFR) or --cbc_len (cLFR) is required.")
    return args


def reverse_complement(seq):
    return seq.translate(str.maketrans("ATGCatgc", "TACGtacg"))[::-1]


def mutate_one_mismatch(seq):
    bases = ("A", "G", "C", "T")
    chars = list(seq)
    for idx in range(10):
        original = chars[idx]
        for base in bases:
            chars[idx] = base
            yield "".join(chars)
        chars[idx] = original


def load_barcodes(path, need_string_hash):
    barcode_hash = {}
    barcode_hash_string = {}
    barcode_count = 0

    with open(path) as handle:
        for line in handle:
            fields = line.split()
            if len(fields) < 2:
                continue
            barcode_count += 1
            barcode, barcode_id = fields[0], fields[1]
            for mutated in mutate_one_mismatch(barcode):
                barcode_hash[mutated] = barcode_id
                if need_string_hash:
                    barcode_hash_string[mutated] = barcode

    return barcode_hash, barcode_hash_string, barcode_count


def init_worker(barcode_hash, barcode_hash_string, options):
    global BARCODE_HASH, BARCODE_HASH_STRING, OPTIONS
    BARCODE_HASH = barcode_hash
    BARCODE_HASH_STRING = barcode_hash_string
    OPTIONS = options


def open_text_gzip(path, mode):
    return gzip.open(path, mode, compresslevel=6)


def read_fastq_record(handle):
    header = handle.readline()
    if not header:
        return None
    seq = handle.readline()
    plus = handle.readline()
    qual = handle.readline()
    if not qual:
        raise ValueError("Truncated FASTQ record")
    return (
        header.rstrip("\n"),
        seq.rstrip("\n"),
        plus.rstrip("\n"),
        qual.rstrip("\n"),
    )


def fastq_chunks(r1_path, r2_path, read_len_r1, chunk_reads):
    with open_text_gzip(r2_path, "rt") as r2_handle:
        r1_handle = open_text_gzip(r1_path, "rt") if read_len_r1 > 0 else None
        try:
            chunk = []
            chunk_index = 0
            while True:
                r2 = read_fastq_record(r2_handle)
                if r2 is None:
                    break
                r1 = read_fastq_record(r1_handle) if r1_handle is not None else None
                chunk.append((r1, r2))
                if len(chunk) >= chunk_reads:
                    yield chunk_index, chunk
                    chunk_index += 1
                    chunk = []
            if chunk:
                yield chunk_index, chunk
        finally:
            if r1_handle is not None:
                r1_handle.close()


def first_id(header):
    return header.split()[0].split("/")[0]


def translate(seq, swap):
    if swap == "ac":
        return seq.translate(str.maketrans("ACac", "CAca"))
    if swap == "gc":
        return seq.translate(str.maketrans("GCgc", "CGcg"))
    return seq


def translate_r2_for_barcode(seq, swap):
    if swap in ("ac", "gc"):
        return seq.translate(str.maketrans("ACac", "CAca"))
    return seq


def substr(seq, start, length):
    return seq[start:start + length]


def make_record(header, seq, plus, qual):
    return "%s\n%s\n%s\n%s\n" % (header, seq, plus, qual)


def empty_outputs(output_mode):
    outputs = {
        "main1": [],
        "main2": [],
        "nobc1": [],
        "nobc2": [],
        "nobc2bc3_1": [],
        "nobc2bc3_2": [],
        "noall3_1": [],
        "noall3_2": [],
        "nobc3_1": [],
        "nobc3_2": [],
        "nobc1_1": [],
        "nobc1_2": [],
    }
    if output_mode in ("single", "cwgs"):
        return {key: outputs[key] for key in ("main1", "main2")}
    if output_mode == "separate":
        return {key: outputs[key] for key in ("main1", "main2", "nobc1", "nobc2")}
    return outputs


def process_chunk(item):
    chunk_index, records = item
    opts = OPTIONS
    outputs = empty_outputs(opts["output_mode"])
    stats = {
        "reads_num": 0,
        "split_reads_num": 0,
        "bc1_cnt": 0,
        "bc2_cnt": 0,
        "bc3_cnt": 0,
        "barcode_reads": defaultdict(int),
        "barcode_order": [],
    }

    bc_start = opts["bc_start"] - 1
    gdna_start = opts["gdna_start"] - 1
    additional_bc_start = opts["additional_bc_start"] - 1
    gdna_start_r1 = opts["gdna_start_r1"] - 1
    read_len = opts["read_len"]
    read_len_r1 = opts["read_len_r1"]
    additional_bc_len = opts["additional_bc_len"]
    output_mode = opts["output_mode"]
    swap = opts["swap"]
    adapter_len = opts["adapter_len"]

    n1, n2, n3, n4, n5 = 10, 6, 10, 6, 10

    for r1, r2 in records:
        stats["reads_num"] += 1
        r2_header, r2_seq_raw, r2_plus, r2_qual = r2
        id_str = first_id(r2_header)
        r2_seq = translate_r2_for_barcode(r2_seq_raw, swap)

        read = substr(r2_seq, gdna_start, read_len)
        b1 = substr(r2_seq, bc_start, n1)
        b2 = substr(r2_seq, bc_start + n1 + n2, n3)
        b3 = substr(r2_seq, bc_start + n1 + n2 + n3 + n4, n5)
        bc_str = substr(r2_seq, bc_start, n1 + n2 + n3 + n4 + n5)
        additional_bc = (
            substr(r2_seq, additional_bc_start, additional_bc_len)
            if additional_bc_len > 0 else ""
        )

        bc1 = BARCODE_HASH.get(b1)
        bc2 = BARCODE_HASH.get(b2)
        bc3 = BARCODE_HASH.get(b3)
        bc1_ok = bc1 is not None
        bc2_ok = bc2 is not None
        bc3_ok = bc3 is not None
        stats["bc1_cnt"] += 1 if bc1_ok else 0
        stats["bc2_cnt"] += 1 if bc2_ok else 0
        stats["bc3_cnt"] += 1 if bc3_ok else 0

        if bc1_ok and bc2_ok and bc3_ok:
            barcode = "%s_%s_%s" % (bc1, bc2, bc3)
            stats["split_reads_num"] += 1
            if barcode not in stats["barcode_reads"]:
                stats["barcode_order"].append(barcode)
            stats["barcode_reads"][barcode] += 1
            hash_string = ""
            if output_mode not in ("stratified", "cwgs"):
                hash_string = (
                    BARCODE_HASH_STRING[b1]
                    + BARCODE_HASH_STRING[b2]
                    + BARCODE_HASH_STRING[b3]
                )

            if read_len_r1 > 0 and r1 is not None:
                r1_header, r1_seq, r1_plus, r1_qual = r1
                r1_seq = translate(r1_seq, swap)
                out_header = "%s#%s#/1\tBX:Z:%s" % (id_str, barcode, barcode)
                if output_mode == "cwgs":
                    out_header = "%s#%s/1" % (id_str, barcode)
                elif output_mode != "stratified":
                    out_header = "%s#%s/1\tBX:Z:%s\tBC:Z:%s" % (
                        id_str, barcode, barcode, hash_string
                    )
                    if additional_bc_len > 0:
                        out_header += "\tMI:Z:%s" % additional_bc
                outputs["main1"].append(make_record(
                    out_header,
                    substr(r1_seq, gdna_start_r1, read_len_r1),
                    r1_plus,
                    substr(r1_qual, gdna_start_r1, read_len_r1),
                ))

            out_header = "%s#%s#/2\tBX:Z:%s" % (id_str, barcode, barcode)
            if output_mode == "cwgs":
                out_header = "%s#%s/2" % (id_str, barcode)
            elif output_mode == "stratified" and additional_bc_len > 0:
                out_header = "%s#%s_%s#/2\tBX:Z:%s" % (
                    id_str, barcode, additional_bc, barcode
                )
            elif output_mode != "stratified":
                out_header = "%s#%s/2\tBX:Z:%s\tBC:Z:%s" % (
                    id_str, barcode, barcode, hash_string
                )
                if additional_bc_len > 0:
                    out_header += "\tMI:Z:%s" % additional_bc
            outputs["main2"].append(make_record(
                out_header,
                read,
                r2_plus,
                substr(r2_qual, gdna_start, read_len),
            ))
            continue

        stats["barcode_reads"]["0_0_0"] += 1
        hash_string = b1 + b2 + b3

        if output_mode in ("single", "separate", "cwgs"):
            out1_key, out2_key = ("main1", "main2")
            if output_mode == "separate":
                out1_key, out2_key = ("nobc1", "nobc2")

            if read_len_r1 > 0 and r1 is not None:
                r1_header, r1_seq, r1_plus, r1_qual = r1
                r1_seq = translate(r1_seq, swap) if output_mode in ("single", "cwgs") else r1_seq
                if output_mode == "cwgs":
                    out_header = "%s#0_0_0/1" % id_str
                else:
                    out_header = "%s#0_0_0/1\tBX:Z:0_0_0\tBC:Z:%s" % (
                        id_str, hash_string
                    )
                if additional_bc_len > 0 and output_mode != "cwgs":
                    out_header += "\tMI:Z:%s" % additional_bc
                outputs[out1_key].append(make_record(
                    out_header,
                    substr(r1_seq, gdna_start_r1, read_len_r1),
                    r1_plus,
                    substr(r1_qual, gdna_start_r1, read_len_r1),
                ))

            if output_mode == "cwgs":
                out_header = "%s#0_0_0/2" % id_str
            else:
                out_header = "%s#0_0_0/2\tBX:Z:0_0_0\tBC:Z:%s" % (
                    id_str, hash_string
                )
            if additional_bc_len > 0 and output_mode != "cwgs":
                out_header += "\tMI:Z:%s" % additional_bc
            outputs[out2_key].append(make_record(
                out_header,
                read,
                r2_plus,
                substr(r2_qual, gdna_start, read_len),
            ))
            continue

        full_len = n1 + n2 + n3 + n4 + n5 + adapter_len + read_len
        full_read = substr(r2_seq, bc_start, full_len)
        partial_barcode = "%s_%s_%s" % (bc1 or "", bc2 or "", bc3 or "")

        if bc1_ok and not bc2_ok and not bc3_ok:
            out1_key, out2_key = "nobc2bc3_1", "nobc2bc3_2"
        elif (not bc1_ok) and (not bc2_ok) and (not bc3_ok):
            out1_key, out2_key = "noall3_1", "noall3_2"
        elif (not bc1_ok) and bc2_ok and bc3_ok:
            out1_key, out2_key = "nobc1_1", "nobc1_2"
        elif bc1_ok and bc2_ok and not bc3_ok:
            out1_key, out2_key = "nobc3_1", "nobc3_2"
        else:
            continue

        if read_len_r1 > 0 and r1 is not None:
            r1_header, r1_seq, r1_plus, r1_qual = r1
            outputs[out1_key].append(make_record(
                "%s#%s#/1\tBX:Z:%s" % (id_str, bc_str, partial_barcode),
                substr(r1_seq, gdna_start_r1, read_len_r1),
                r1_plus,
                substr(r1_qual, gdna_start_r1, read_len_r1),
            ))
        out_header = "%s#%s#/2\tBX:Z:%s" % (id_str, bc_str, partial_barcode)
        if additional_bc_len > 0:
            out_header = "%s#%s_%s#/2\tBX:Z:%s" % (
                id_str, bc_str, additional_bc, partial_barcode
            )
        outputs[out2_key].append(make_record(
            out_header,
            full_read,
            r2_plus,
            substr(r2_qual, bc_start, full_len),
        ))

    return chunk_index, {key: "".join(value) for key, value in outputs.items()}, stats


def process_clfr_chunk(item):
    chunk_index, records = item
    opts = OPTIONS
    outputs = {"main1": [], "main2": []}
    stats = {
        "reads_num": 0,
        "split_reads_num": 0,
        "barcode_reads": defaultdict(int),
        "barcode_order": [],
    }

    bc_start = opts["bc_start"] - 1
    gdna_start = opts["gdna_start"] - 1
    additional_bc_start = opts["additional_bc_start"] - 1
    gdna_start_r1 = opts["gdna_start_r1"] - 1
    read_len = opts["read_len"]
    read_len_r1 = opts["read_len_r1"]
    cbc_len = opts["cbc_len"]
    bc_len_redundant = opts["bc_len_redundant"]
    additional_bc_len = opts["additional_bc_len"]
    additional_bc_len_r1 = opts["additional_bc_len_r1"]
    use_reverse_complement = opts["reverse_complement"]

    for r1, r2 in records:
        stats["reads_num"] += 1
        r2_header, r2_seq, r2_plus, r2_qual = r2
        read_id = first_id(r2_header)

        if use_reverse_complement:
            bc = reverse_complement(substr(r2_seq, bc_start, bc_len_redundant))
            bx = reverse_complement(substr(r2_seq, bc_start, cbc_len))
        else:
            bc = substr(r2_seq, bc_start, bc_len_redundant)
            bx = substr(r2_seq, bc_start, cbc_len)

        if bx not in stats["barcode_reads"]:
            stats["barcode_order"].append(bx)
        stats["barcode_reads"][bx] += 1
        stats["split_reads_num"] += 1

        r2_read = substr(r2_seq, gdna_start, read_len)
        r2_qual_out = substr(r2_qual, gdna_start, read_len)

        if read_len_r1 > 0 and r1 is not None:
            r1_header, r1_seq, r1_plus, r1_qual = r1
            if additional_bc_len > 0:
                additional_bc = substr(r2_seq, additional_bc_start, additional_bc_len)
                header1 = "%s#%s#%s/1\tBX:Z:%s" % (read_id, bc, additional_bc, bx)
                header2 = "%s#%s#%s/2\tBX:Z:%s" % (read_id, bc, additional_bc, bx)
            elif additional_bc_len_r1 > 0:
                additional_bc = substr(r1_seq, additional_bc_start, additional_bc_len_r1)
                header1 = "%s#%s#%s/1\tBX:Z:%s" % (read_id, bc, additional_bc, bx)
                header2 = "%s#%s#%s/2\tBX:Z:%s" % (read_id, bc, additional_bc, bx)
            else:
                header1 = "%s#%s/1\tBX:Z:%s" % (read_id, bc, bx)
                header2 = "%s#%s/2\tBX:Z:%s" % (read_id, bc, bx)

            outputs["main1"].append(make_record(
                header1,
                substr(r1_seq, gdna_start_r1, read_len_r1),
                r1_plus,
                substr(r1_qual, gdna_start_r1, read_len_r1),
            ))
            outputs["main2"].append(make_record(header2, r2_read, r2_plus, r2_qual_out))
        else:
            if additional_bc_len > 0:
                additional_bc = substr(r2_seq, additional_bc_start, additional_bc_len)
                header2 = "%s#%s#%s/2\tBX:Z:%s" % (read_id, bc, additional_bc, bx)
            else:
                header2 = "%s#%s/2\tBX:Z:%s" % (read_id, bc, bx)
            outputs["main2"].append(make_record(header2, r2_read, r2_plus, r2_qual_out))

    return chunk_index, {key: "".join(value) for key, value in outputs.items()}, stats


def output_paths(output_prefix, output_mode):
    paths = {
        "main1": output_prefix + ".1.fq.gz",
        "main2": output_prefix + ".2.fq.gz",
    }
    if output_mode == "separate":
        paths.update({
            "nobc1": "data/noBC.1.fq.gz",
            "nobc2": "data/noBC.2.fq.gz",
        })
    elif output_mode == "stratified":
        paths.update({
            "nobc2bc3_1": output_prefix + "_nobc2bc3.1.fq.gz",
            "nobc2bc3_2": output_prefix + "_nobc2bc3.2.fq.gz",
            "noall3_1": output_prefix + "_noall3.1.fq.gz",
            "noall3_2": output_prefix + "_noall3.2.fq.gz",
            "nobc3_1": output_prefix + "_nobc3.1.fq.gz",
            "nobc3_2": output_prefix + "_nobc3.2.fq.gz",
            "nobc1_1": output_prefix + "_nobc1.1.fq.gz",
            "nobc1_2": output_prefix + "_nobc1.2.fq.gz",
        })
    return paths


def write_logs(stats, barcode_reads, barcode_order, barcode_count):
    barcode_types = barcode_count * barcode_count * barcode_count
    split_barcode_num = len([bc for bc in barcode_reads if bc != "0_0_0"])
    reads_num = stats["reads_num"]
    split_reads_num = stats["split_reads_num"]

    with open("split_stat_read1.log", "w") as out:
        out.write("Barcode_types = %s * %s * %s = %s\n" % (
            barcode_count, barcode_count, barcode_count, barcode_types
        ))
        barcode_rate = 0 if barcode_types == 0 else 100 * split_barcode_num / barcode_types
        reads_rate = 0 if reads_num == 0 else 100 * split_reads_num / reads_num
        out.write("Real_Barcode_types = %s (%s %%)\n" % (split_barcode_num, barcode_rate))
        out.write("Reads_pair_num  = %s \n" % reads_num)
        out.write("Reads_pair_num(after split) = %s (%s %%)\n" % (
            split_reads_num, reads_rate
        ))
        for index, barcode in enumerate(barcode_order, 1):
            out.write("%s\t%s\t%s\n" % (index, barcode_reads[barcode], barcode))

    with open("bc.individual.count.log", "w") as out:
        total_reads = reads_num if stats["read_len_r1"] == 0 else reads_num * 2
        out.write("%s\n" % total_reads)
        out.write("bc1_cnt=%s \n" % stats["bc1_cnt"])
        out.write("bc2_cnt=%s \n" % stats["bc2_cnt"])
        out.write("bc3_cnt=%s \n" % stats["bc3_cnt"])
        out.write("%s\t%s\t%s\n" % (stats["bc1_cnt"], stats["bc2_cnt"], stats["bc3_cnt"]))
        if reads_num == 0:
            out.write("0\t0\t0\n")
        else:
            out.write("%s\t%s\t%s\n" % (
                stats["bc1_cnt"] / reads_num,
                stats["bc2_cnt"] / reads_num,
                stats["bc3_cnt"] / reads_num,
            ))


def write_clfr_logs(stats, barcode_reads, barcode_order, cbc_len):
    reads_num = stats["reads_num"]
    split_reads_num = stats["split_reads_num"]
    split_barcode_num = len(barcode_order)
    barcode_types = 4 ** cbc_len

    with open("split_stat_read1.log", "w") as out:
        out.write("Barcode_types = 4** %s = %s\n" % (cbc_len, barcode_types))
        barcode_rate = 0 if barcode_types == 0 else 100 * split_barcode_num / barcode_types
        reads_rate = 0 if reads_num == 0 else 100 * split_barcode_num / reads_num
        out.write("Real_Barcode_types = %s (%s %%)\n" % (split_barcode_num, barcode_rate))
        out.write("Reads_pair_num  = %s \n" % reads_num)
        out.write("Reads_pair_num(after split) = %s (%s %%)\n" % (
            split_reads_num, reads_rate
        ))
        for index, barcode in enumerate(barcode_order, 1):
            out.write("%s\t%s\t%s\n" % (index, barcode_reads[barcode], barcode))

    with open("bc.diversity.count.log", "w") as out:
        total_reads = reads_num if stats["read_len_r1"] == 0 else reads_num * 2
        diversity = 0 if reads_num == 0 else split_barcode_num / reads_num
        out.write("%s\n" % total_reads)
        out.write("%s\n" % diversity)


def main():
    args = parse_args()
    threads = max(1, args.threads)
    is_clfr = bool(args.cbc_len)
    barcode_count = 0
    if is_clfr:
        barcode_hash, barcode_hash_string = {}, {}
    else:
        need_string_hash = args.output_mode not in ("stratified", "cwgs")
        barcode_hash, barcode_hash_string, barcode_count = load_barcodes(
            args.barcode, need_string_hash
        )
    options = vars(args).copy()
    output_handles = {}
    paths = output_paths(args.output, args.output_mode)
    for key, path in paths.items():
        directory = os.path.dirname(path)
        if directory:
            os.makedirs(directory, exist_ok=True)
        output_handles[key] = open_text_gzip(path, "wt")

    aggregate = {
        "reads_num": 0,
        "split_reads_num": 0,
        "bc1_cnt": 0,
        "bc2_cnt": 0,
        "bc3_cnt": 0,
        "read_len_r1": args.read_len_r1,
    }
    barcode_reads = defaultdict(int)
    barcode_order = []
    barcode_seen = set()
    pool = None

    try:
        chunks = fastq_chunks(args.r1, args.r2, args.read_len_r1, args.chunk_reads)
        worker = process_clfr_chunk if is_clfr else process_chunk
        if threads == 1:
            init_worker(barcode_hash, barcode_hash_string, options)
            iterator = map(worker, chunks)
        else:
            pool = mp.Pool(
                processes=threads,
                initializer=init_worker,
                initargs=(barcode_hash, barcode_hash_string, options),
            )
            iterator = pool.imap(worker, chunks, chunksize=1)

        for chunk_index, outputs, stats in iterator:
            if chunk_index % 100 == 0:
                print("chunks processed %s ..." % chunk_index, file=sys.stderr)
            for key, payload in outputs.items():
                output_handles[key].write(payload)
            for key in ("reads_num", "split_reads_num", "bc1_cnt", "bc2_cnt", "bc3_cnt"):
                if key in stats:
                    aggregate[key] += stats[key]
            for barcode, count in stats["barcode_reads"].items():
                barcode_reads[barcode] += count
            for barcode in stats["barcode_order"]:
                if barcode not in barcode_seen:
                    barcode_seen.add(barcode)
                    barcode_order.append(barcode)

        if threads > 1:
            pool.close()
            pool.join()
            pool = None
    finally:
        if pool is not None:
            pool.terminate()
        for handle in output_handles.values():
            handle.close()

    if is_clfr:
        write_clfr_logs(aggregate, barcode_reads, barcode_order, args.cbc_len)
    else:
        write_logs(aggregate, barcode_reads, barcode_order, barcode_count)
    print("all done!")


if __name__ == "__main__":
    main()
