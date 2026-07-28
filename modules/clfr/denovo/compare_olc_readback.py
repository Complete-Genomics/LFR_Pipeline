#!/usr/bin/env python3
"""Read-back support for primary OLC contigs lengthened by a new run.

This deliberately does not treat MEGAHIT as truth.  A read supports a contig
only after a unique-orientation placement with >=3 exact 17-mer anchors
spanning >=40 bp; support metrics are then computed directly from the input
reads.
"""

import argparse
import csv
from collections import defaultdict

from analyze_olc_shortfall import read_first_barcode_groups, kmer_index, place_either


def load_fasta(path):
    records = defaultdict(list)
    barcode = None
    chunks = []
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if line.startswith(">"): 
                if barcode is not None:
                    records[barcode].append("".join(chunks))
                barcode = line[1:].split()[0].split(">", 1)[0]
                chunks = []
            elif line:
                chunks.append(line)
    if barcode is not None:
        records[barcode].append("".join(chunks))
    return records


def support(contig, reads):
    index = kmer_index(contig, 17)
    diff = [0] * (len(contig) + 1)
    placed = left_end = right_end = 0
    for read in set(reads):
        hit = place_either(read, contig, index)
        if hit is None:
            continue
        start, end, _, _ = hit
        lo, hi = max(0, start), min(len(contig), end)
        if lo >= hi:
            continue
        placed += 1
        diff[lo] += 1
        diff[hi] -= 1
        if start <= 0 and end >= 20:
            left_end += 1
        if start <= len(contig) - 20 and end >= len(contig):
            right_end += 1
    coverage = []
    value = 0
    for delta in diff[:-1]:
        value += delta
        coverage.append(value)
    return {
        "placed_reads": placed,
        "breadth_1x": sum(x >= 1 for x in coverage) / len(contig),
        "breadth_2x": sum(x >= 2 for x in coverage) / len(contig),
        "left_end_reads": left_end,
        "right_end_reads": right_end,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--r2", required=True)
    parser.add_argument("--baseline", required=True)
    parser.add_argument("--candidate", required=True)
    parser.add_argument("--n", type=int, default=20000)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    reads = read_first_barcode_groups(args.r2, args.n)
    baseline = load_fasta(args.baseline)
    candidate = load_fasta(args.candidate)
    rows = []
    for barcode, umi_reads in reads.items():
        old = max(baseline.get(barcode, []), key=len, default="")
        new = max(candidate.get(barcode, []), key=len, default="")
        if len(new) <= len(old):
            continue
        row = {
            "barcode": barcode,
            "old_primary_len": len(old),
            "new_primary_len": len(new),
            "gain_bp": len(new) - len(old),
            "old_component_count": len(baseline.get(barcode, [])),
            "new_component_count": len(candidate.get(barcode, [])),
            "umi_reads_raw": len(umi_reads),
            "umi_reads_unique": len(set(umi_reads)),
        }
        for prefix, contig in (("old", old), ("new", new)):
            for key, value in support(contig, umi_reads).items():
                row["{}_{}".format(prefix, key)] = value
        row["new_two_end_supported"] = int(
            row["new_left_end_reads"] >= 2 and row["new_right_end_reads"] >= 2)
        rows.append(row)

    with open(args.out, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]) if rows else [], delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print("lengthened_primary_contigs={}".format(len(rows)))
    print("new_two_end_supported={}".format(sum(r["new_two_end_supported"] for r in rows)))


if __name__ == "__main__":
    main()
