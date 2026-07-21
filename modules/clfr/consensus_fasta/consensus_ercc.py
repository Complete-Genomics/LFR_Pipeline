#!/usr/bin/env python3
"""Summarize ERCC consensus FASTA length recovery and variants.

The script compares cLFR consensus FASTA records against ERCC truth sequences.
It reports how much of the original ERCC length each consensus recovers and
whether the consensus isoform carries SNPs/indels relative to the ERCC reference.
"""

from __future__ import print_function

import argparse
import csv
import os
import re
import sys


ERCC_RE = re.compile(r"(ERCC-\d+)")


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--fasta", required=True, help="Consensus FASTA file.")
    parser.add_argument(
        "--ercc_ref",
        default=None,
        help="ERCC truth table with columns 'ERCC ID' and 'Sequence'.",
    )
    parser.add_argument(
        "--summary",
        default="consensus/consensus_ercc_summary.tsv",
        help="Output per-consensus summary TSV.",
    )
    parser.add_argument(
        "--variants",
        default="consensus/consensus_ercc_variants.tsv",
        help="Output per-variant TSV.",
    )
    parser.add_argument("--match", type=int, default=2)
    parser.add_argument("--mismatch", type=int, default=-1)
    parser.add_argument("--gap", type=int, default=-2)
    return parser.parse_args()


def default_ercc_ref():
    here = os.path.dirname(os.path.abspath(__file__))
    return os.path.abspath(
        os.path.join(here, "..", "calc_frag_len", "ercc_truth.txt")
    )


def read_fasta(path):
    name = None
    seq_chunks = []
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq_chunks).upper()
                name = line[1:].split()[0]
                seq_chunks = []
            else:
                seq_chunks.append(line)
    if name is not None:
        yield name, "".join(seq_chunks).upper()


def read_ercc_truth(path):
    refs = {}
    with open(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if "ERCC ID" not in reader.fieldnames or "Sequence" not in reader.fieldnames:
            raise ValueError("ERCC reference must contain 'ERCC ID' and 'Sequence'")
        for row in reader:
            ercc_id = row["ERCC ID"].strip()
            seq = row["Sequence"].strip().upper()
            if ercc_id and seq:
                refs[ercc_id] = seq
    return refs


def ercc_from_header(header):
    match = ERCC_RE.search(header)
    return match.group(1) if match else None


def kmer_score(query, ref, k=15):
    if len(query) < k or len(ref) < k:
        return 0
    kmers = set(query[i:i + k] for i in range(0, len(query) - k + 1, k))
    return sum(1 for i in range(0, len(ref) - k + 1, k) if ref[i:i + k] in kmers)


def choose_reference(header, query, refs):
    ercc_id = ercc_from_header(header)
    if ercc_id in refs:
        return ercc_id
    return max(refs, key=lambda ref_id: kmer_score(query, refs[ref_id]))


def smith_waterman(ref, query, match_score, mismatch_score, gap_score):
    """Local alignment of query consensus to reference.

    Trace codes:
      0 stop, 1 diag, 2 up/ref gap in query, 3 left/query insertion.
    """
    n_ref = len(ref)
    n_query = len(query)
    prev = [0] * (n_query + 1)
    trace = [bytearray(n_query + 1) for _ in range(n_ref + 1)]
    max_score = 0
    max_i = 0
    max_j = 0

    for i in range(1, n_ref + 1):
        curr = [0] * (n_query + 1)
        ref_base = ref[i - 1]
        for j in range(1, n_query + 1):
            query_base = query[j - 1]
            diag = prev[j - 1] + (match_score if ref_base == query_base else mismatch_score)
            up = prev[j] + gap_score
            left = curr[j - 1] + gap_score
            best = max(0, diag, up, left)
            curr[j] = best
            if best == 0:
                trace[i][j] = 0
            elif best == diag:
                trace[i][j] = 1
            elif best == up:
                trace[i][j] = 2
            else:
                trace[i][j] = 3
            if best > max_score:
                max_score = best
                max_i = i
                max_j = j
        prev = curr

    ops = []
    i = max_i
    j = max_j
    while i > 0 and j > 0 and trace[i][j] != 0:
        code = trace[i][j]
        if code == 1:
            ops.append(("M", i - 1, j - 1))
            i -= 1
            j -= 1
        elif code == 2:
            ops.append(("D", i - 1, None))
            i -= 1
        elif code == 3:
            ops.append(("I", None, j - 1))
            j -= 1
    ops.reverse()
    return {
        "score": max_score,
        "ops": ops,
        "ref_start": i + 1,
        "ref_end": max_i,
        "query_start": j + 1,
        "query_end": max_j,
    }


def summarize_alignment(header, ercc_id, ref, query, aln):
    variants = []
    snp_count = 0
    ins_count = 0
    del_count = 0
    n_count = query.count("N")
    aligned_ref_bases = 0
    aligned_query_bases = 0

    for op, ref_idx, query_idx in aln["ops"]:
        if op == "M":
            aligned_ref_bases += 1
            aligned_query_bases += 1
            ref_base = ref[ref_idx]
            query_base = query[query_idx]
            if ref_base != query_base:
                if query_base == "N":
                    variant_type = "ambiguous_N"
                else:
                    variant_type = "snp"
                    snp_count += 1
                variants.append({
                    "consensus_id": header,
                    "ercc_id": ercc_id,
                    "variant_type": variant_type,
                    "ref_pos_1based": ref_idx + 1,
                    "ref_base": ref_base,
                    "consensus_base": query_base,
                })
        elif op == "D":
            aligned_ref_bases += 1
            del_count += 1
            variants.append({
                "consensus_id": header,
                "ercc_id": ercc_id,
                "variant_type": "deletion",
                "ref_pos_1based": ref_idx + 1,
                "ref_base": ref[ref_idx],
                "consensus_base": "-",
            })
        elif op == "I":
            aligned_query_bases += 1
            ins_count += 1
            variants.append({
                "consensus_id": header,
                "ercc_id": ercc_id,
                "variant_type": "insertion",
                "ref_pos_1based": "",
                "ref_base": "-",
                "consensus_base": query[query_idx],
            })

    ref_len = len(ref)
    query_len = len(query)
    summary = {
        "consensus_id": header,
        "ercc_id": ercc_id,
        "ref_len": ref_len,
        "consensus_len": query_len,
        "consensus_len_pct_ref": 100.0 * query_len / ref_len if ref_len else 0,
        "aligned_ref_bases": aligned_ref_bases,
        "aligned_ref_pct": 100.0 * aligned_ref_bases / ref_len if ref_len else 0,
        "aligned_query_bases": aligned_query_bases,
        "alignment_score": aln["score"],
        "ref_start_1based": aln["ref_start"],
        "ref_end_1based": aln["ref_end"],
        "query_start_1based": aln["query_start"],
        "query_end_1based": aln["query_end"],
        "snp_count": snp_count,
        "insertion_count": ins_count,
        "deletion_count": del_count,
        "ambiguous_N_count": n_count,
        "has_snp": "yes" if snp_count > 0 else "no",
        "has_indel": "yes" if (ins_count + del_count) > 0 else "no",
    }
    return summary, variants


def write_table(path, rows, fieldnames):
    directory = os.path.dirname(path)
    if directory:
        os.makedirs(directory, exist_ok=True)
    with open(path, "w") as out:
        writer = csv.DictWriter(out, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def main():
    args = parse_args()
    ercc_ref = args.ercc_ref or default_ercc_ref()
    refs = read_ercc_truth(ercc_ref)
    summaries = []
    variants = []

    for header, query in read_fasta(args.fasta):
        ercc_id = choose_reference(header, query, refs)
        aln = smith_waterman(
            refs[ercc_id],
            query,
            match_score=args.match,
            mismatch_score=args.mismatch,
            gap_score=args.gap,
        )
        summary, record_variants = summarize_alignment(header, ercc_id, refs[ercc_id], query, aln)
        summaries.append(summary)
        variants.extend(record_variants)

    summary_fields = [
        "consensus_id", "ercc_id", "ref_len", "consensus_len", "consensus_len_pct_ref",
        "aligned_ref_bases", "aligned_ref_pct", "aligned_query_bases", "alignment_score",
        "ref_start_1based", "ref_end_1based", "query_start_1based", "query_end_1based",
        "snp_count", "insertion_count", "deletion_count", "ambiguous_N_count",
        "has_snp", "has_indel",
    ]
    variant_fields = [
        "consensus_id", "ercc_id", "variant_type", "ref_pos_1based",
        "ref_base", "consensus_base",
    ]
    write_table(args.summary, summaries, summary_fields)
    write_table(args.variants, variants, variant_fields)
    sys.stderr.write("Wrote %s consensus summaries and %s variants.\n" % (
        len(summaries), len(variants)
    ))


if __name__ == "__main__":
    main()
