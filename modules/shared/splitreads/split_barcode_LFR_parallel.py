#!/usr/bin/env python3
"""Run split_barcode_LFR.py (cLFR mode) on several slices of the input in parallel.

split_barcode_LFR.py is bounded by one main process (gzip read + fastq parse + merge of per-barcode
counts), and that process slows down as its barcode tables grow, so extra --threads barely help.
This wrapper cuts R2 (and R1 when --read_len_r1 > 0) into --parts read-aligned slices, runs one
split_barcode_LFR.py per slice at the same time, then rebuilds exactly what a single run produces:
  <output>.1.fq.gz / <output>.2.fq.gz  (slice outputs concatenated in input order)
  split_stat_read1.log, bc.diversity.count.log  (per-barcode counts summed across slices)

Slice i is started as soon as its input is written, so the first slices are already running while
the rest of the input is still being cut.  Unknown options are passed through to split_barcode_LFR.py.

usage: split_barcode_LFR_parallel.py --parts 4 [--total_reads N] --threads 40 \
           <the usual split_barcode_LFR.py options: --r1 --r2 --output --cbc_len ...>
"""

import argparse
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
SPLIT_SCRIPT = os.path.join(HERE, "split_barcode_LFR.py")
BLOCK = 1 << 26
BARCODE_ORDER_STRIDE = 10 ** 12  # part * stride + per-part index = global first-seen order


def parse_args():
    # no abbreviations: "--read_len" must pass through to the child, not match "--read_len_r1"
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0], allow_abbrev=False)
    parser.add_argument("--parts", type=int, default=4)
    parser.add_argument("--total_reads", type=int, default=0,
                        help="Read pairs in the input; 0 = count them (one extra pass over R2).")
    parser.add_argument("--r1", required=True)
    parser.add_argument("--r2", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--read_len_r1", type=int, required=True)
    parser.add_argument("--cbc_len", type=int, required=True)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--tmpdir")
    parser.add_argument("--pigz", default="pigz")
    parser.add_argument("--sort_mem", default="4G")
    args, passthrough = parser.parse_known_args()
    if args.parts < 1:
        parser.error("--parts must be >= 1")
    return args, passthrough


def tool(name, fallback):
    return name if shutil.which(name) else fallback


def decompress_cmd(args, path):
    return [tool(args.pigz, "gzip"), "-dc", path]


def compress_cmd(args, threads):
    if shutil.which(args.pigz):
        return [args.pigz, "-p", str(max(1, threads)), "-c"]
    return ["gzip", "-c"]


def count_reads(args, path):
    lines = 0
    proc = subprocess.Popen(decompress_cmd(args, path), stdout=subprocess.PIPE, bufsize=BLOCK)
    while True:
        block = proc.stdout.read(BLOCK)
        if not block:
            break
        lines += block.count(b"\n")
    if proc.wait() != 0:
        sys.exit("ERROR: failed to read %s" % path)
    if lines % 4:
        sys.exit("ERROR: %s has %d lines, not a multiple of 4" % (path, lines))
    return lines // 4


def cut_fastq(args, path, lines_per_part, part_path, on_part_done, pigz_threads):
    """Write consecutive slices of `path` (lines_per_part lines each) to part_path(i)."""
    reader = subprocess.Popen(decompress_cmd(args, path), stdout=subprocess.PIPE, bufsize=BLOCK)
    writer = out_handle = None
    part = -1
    part_lines = 0

    def open_next():
        nonlocal writer, out_handle, part, part_lines
        part += 1
        part_lines = 0
        out_handle = open(part_path(part), "wb")
        writer = subprocess.Popen(compress_cmd(args, pigz_threads), stdin=subprocess.PIPE,
                                  stdout=out_handle)

    def close_current():
        nonlocal writer, out_handle
        writer.stdin.close()
        if writer.wait() != 0:
            sys.exit("ERROR: compressing slice %d of %s failed" % (part, path))
        out_handle.close()
        writer = out_handle = None
        on_part_done(part)

    while True:
        block = reader.stdout.read(BLOCK)
        if not block:
            break
        pos = 0
        while pos < len(block):
            if writer is None:
                open_next()
            need = lines_per_part - part_lines
            available = block.count(b"\n", pos)
            if available < need:
                writer.stdin.write(block[pos:])
                part_lines += available
                pos = len(block)
            else:
                end = pos - 1
                for _ in range(need):
                    end = block.find(b"\n", end + 1)
                writer.stdin.write(block[pos:end + 1])
                pos = end + 1
                close_current()
    if writer is not None:
        close_current()
    if reader.wait() != 0:
        sys.exit("ERROR: failed to read %s" % path)
    return part + 1


def merge_gz(inputs, output):
    """Concatenated gzip members are a valid gzip file."""
    with open(output, "wb") as out:
        for path in inputs:
            with open(path, "rb") as handle:
                shutil.copyfileobj(handle, out, length=1 << 22)
            os.unlink(path)


def read_log_header(path):
    reads = split = None
    with open(path) as handle:
        for line in handle:
            match = re.match(r"Reads_pair_num\s+=\s+(\d+)", line)
            if match:
                reads = int(match.group(1))
            match = re.match(r"Reads_pair_num\(after split\)\s+=\s+(\d+)", line)
            if match:
                split = int(match.group(1))
            if reads is not None and split is not None:
                break
    if reads is None or split is None:
        sys.exit("ERROR: cannot parse %s" % path)
    return reads, split


def sort_supports_parallel():
    return subprocess.call(["sort", "--parallel=2"], stdin=subprocess.DEVNULL,
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL) == 0


def merge_barcode_counts(part_logs, workdir, sort_mem):
    """Sum per-barcode read counts over slices with on-disk sorts (no in-memory table).

    Slice log lines (after the 4 header lines) are "<index>\t<reads>\t<barcode>", barcodes in
    first-seen order.  The merged file keeps global first-seen order and renumbers from 1.
    Returns (distinct_barcodes, path of merged body).
    """
    parallel = " --parallel=8" if sort_supports_parallel() else ""
    sort = "LC_ALL=C sort -S %s -T %s%s" % (sort_mem, workdir, parallel)
    body = os.path.join(workdir, "merged_barcodes.tsv")
    feed = "; ".join(
        "awk -F'\\t' -v p=%d 'NR > 4 {printf \"%%s\\t%%d\\t%%.0f\\n\", $3, $2, p * %d + $1}' %s"
        % (i, BARCODE_ORDER_STRIDE, path)
        for i, path in enumerate(part_logs)
    )
    # barcode, reads, order -> one row per barcode: order, reads, barcode -> sorted by order
    group = (
        "awk -F'\\t' 'BEGIN {OFS=\"\\t\"} "
        "$1 != cur {if (NR > 1) printf \"%.0f\\t%d\\t%s\\n\", first, total, cur; cur=$1; total=0; first=$3} "
        "{total += $2; if ($3 < first) first = $3} "
        "END {if (NR > 0) printf \"%.0f\\t%d\\t%s\\n\", first, total, cur}'"
    )
    renumber = "awk -F'\\t' 'BEGIN {OFS=\"\\t\"} {print NR, $2, $3}'"
    pipeline = "set -o pipefail; ({feed}) | {sort} -t $'\\t' -k1,1 | {group} | {sort} -t $'\\t' -k1,1n | {renumber} > {body}".format(
        feed=feed, sort=sort, group=group, renumber=renumber, body=body)
    subprocess.check_call(["bash", "-c", pipeline])
    with open(body, "rb") as handle:
        distinct = sum(1 for _ in handle)
    return distinct, body


def write_merged_logs(args, part_logs, workdir):
    reads_num = split_reads_num = 0
    for path in part_logs:
        reads, split = read_log_header(path)
        reads_num += reads
        split_reads_num += split
    distinct, body = merge_barcode_counts(part_logs, workdir, args.sort_mem)

    barcode_types = 4 ** args.cbc_len
    barcode_rate = 0 if barcode_types == 0 else 100 * distinct / barcode_types
    reads_rate = 0 if reads_num == 0 else 100 * distinct / reads_num
    with open("split_stat_read1.log", "w") as out:
        out.write("Barcode_types = 4** %s = %s\n" % (args.cbc_len, barcode_types))
        out.write("Real_Barcode_types = %s (%s %%)\n" % (distinct, barcode_rate))
        out.write("Reads_pair_num  = %s \n" % reads_num)
        out.write("Reads_pair_num(after split) = %s (%s %%)\n" % (split_reads_num, reads_rate))
        out.flush()
        with open(body, "rb") as handle:
            shutil.copyfileobj(handle, out.buffer if hasattr(out, "buffer") else out)
    with open("bc.diversity.count.log", "w") as out:
        total_reads = reads_num if args.read_len_r1 == 0 else reads_num * 2
        diversity = 0 if reads_num == 0 else distinct / reads_num
        out.write("%s\n" % total_reads)
        out.write("%s\n" % diversity)


def main():
    args, passthrough = parse_args()
    output = os.path.abspath(args.output)
    r1 = os.path.abspath(args.r1)
    r2 = os.path.abspath(args.r2)
    use_r1 = args.read_len_r1 > 0
    os.makedirs(os.path.dirname(output), exist_ok=True)
    workdir = tempfile.mkdtemp(prefix=".split_parallel.",
                               dir=args.tmpdir or os.path.dirname(output))
    final_cwd = os.getcwd()

    total_reads = args.total_reads or count_reads(args, r2)
    if total_reads == 0:
        sys.exit("ERROR: %s is empty" % r2)
    lines_per_part = math.ceil(total_reads / args.parts) * 4
    n_parts = math.ceil(total_reads * 4 / lines_per_part)
    threads_per_part = max(1, args.threads // n_parts)
    print("split_parallel: %d read pairs -> %d slices of %d reads, %d threads each"
          % (total_reads, n_parts, lines_per_part // 4, threads_per_part), file=sys.stderr)

    def part_dir(i):
        return os.path.join(workdir, "p%d" % i)

    def part_in(i, mate):
        return os.path.join(part_dir(i), "in.%d.fq.gz" % mate)

    for i in range(n_parts):
        os.makedirs(part_dir(i))

    procs = {}

    def start_part(i):
        cmd = [sys.executable, SPLIT_SCRIPT] + passthrough + [
            "--r1", part_in(i, 1) if use_r1 else r1,
            "--r2", part_in(i, 2),
            "--output", os.path.join(part_dir(i), "split_read"),
            "--read_len_r1", str(args.read_len_r1),
            "--cbc_len", str(args.cbc_len),
            "--threads", str(threads_per_part),
            "--pigz", args.pigz,
        ]
        err = open(os.path.join(part_dir(i), "split.err"), "w")
        procs[i] = subprocess.Popen(cmd, cwd=part_dir(i), stderr=err)
        print("split_parallel: slice %d started" % i, file=sys.stderr)

    def fail(message):
        for proc in procs.values():
            if proc.poll() is None:
                proc.terminate()
        sys.exit("ERROR: %s (kept %s for debugging)" % (message, workdir))

    def check_children():
        for i, proc in procs.items():
            if proc.poll() not in (None, 0):
                tail = open(os.path.join(part_dir(i), "split.err")).read()[-2000:]
                fail("slice %d failed:\n%s" % (i, tail))

    if use_r1:
        cut_fastq(args, r1, lines_per_part, lambda i: part_in(i, 1), lambda i: None, threads_per_part)

    def r2_path(i):
        if i >= n_parts:
            fail("R2 has more reads than --total_reads=%d" % total_reads)
        return part_in(i, 2)

    def r2_done(i):
        check_children()
        start_part(i)

    cut = cut_fastq(args, r2, lines_per_part, r2_path, r2_done, threads_per_part)
    if cut != n_parts:
        fail("expected %d slices, wrote %d" % (n_parts, cut))

    for i in range(n_parts):
        if procs[i].wait() != 0:
            check_children()
            fail("slice %d failed" % i)
        for mate in (1, 2):
            path = part_in(i, mate)
            if os.path.exists(path):
                os.unlink(path)
        print("split_parallel: slice %d done" % i, file=sys.stderr)
        with open(os.path.join(part_dir(i), "split.err")) as handle:
            sys.stderr.write(handle.read())

    for mate in (1, 2):
        merge_gz([os.path.join(part_dir(i), "split_read.%d.fq.gz" % mate) for i in range(n_parts)],
                 "%s.%d.fq.gz" % (output, mate))
    os.chdir(final_cwd)
    write_merged_logs(args, [os.path.join(part_dir(i), "split_stat_read1.log") for i in range(n_parts)],
                      workdir)
    shutil.rmtree(workdir, ignore_errors=True)
    print("all done!")


if __name__ == "__main__":
    main()
