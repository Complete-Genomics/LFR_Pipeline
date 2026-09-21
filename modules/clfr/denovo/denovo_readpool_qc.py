#!/usr/bin/env python3
"""Pre-assembly read-pool consistency QC (reference-free).

For each barcode, BEFORE assembly, check whether its raw reads collapse into
one connected group under a STRICT overlap requirement, rather than the
assembler's own (looser) merge threshold.

Why stricter than the assembler: denovo.md sec 28 found that the cross-
species chimeras in the ZymoBIOMICS control formed via short, highly
conserved 16S anchor regions -- exactly the kind of overlap the assembler's
own (~20bp minimum) overlap criterion will happily accept. A read-pool QC
check reusing that same threshold would have the identical blind spot: two
reads from different source molecules that merely share a conserved primer-
adjacent patch would still look "connected". Requiring a much longer
matching span (default 150bp) before treating two reads as same-molecule
makes coincidental conserved-region bridging between different species much
less likely to pass by chance, while genuine same-molecule reads (true
overlapping fragments of one physical source DNA) should support long,
high-identity matches, not just a short conserved island.

This is a *pre*-assembly signal: it flags barcodes whose read pool itself
looks mixed, independent of what the assembler eventually does with those
reads -- unlike post-assembly read-back QC (readback_qc_single.py), which
cannot see this problem because a chimera built from a genuinely mixed pool
is, by construction, well-supported by that same pool (denovo.md sec 30).
"""
import argparse
import os
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + "/test")
from analyze_olc_shortfall import read_first_barcode_groups, kmer_index, rc  # noqa: E402


def strict_overlap(read_a, index_b, index_b_rc, k, min_anchors, min_span):
    """True if read_a shares a long, high-confidence overlap with the read
    that index_b/index_b_rc were built from, in either orientation."""
    for index in (index_b, index_b_rc):
        offsets = defaultdict(list)
        for qpos in range(len(read_a) - k + 1):
            hits = index.get(read_a[qpos:qpos + k], ())
            if len(hits) > 12:
                continue
            for tpos in hits:
                offsets[tpos - qpos].append(qpos)
        for positions in offsets.values():
            unique = sorted(set(positions))
            if len(unique) >= min_anchors and unique[-1] - unique[0] >= min_span:
                return True
    return False


def connected_components(reads, k, min_anchors, min_span, max_reads):
    reads = reads[:max_reads]
    n = len(reads)
    indexes = [kmer_index(r, k) for r in reads]
    indexes_rc = [kmer_index(rc(r), k) for r in reads]

    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(x, y):
        rx, ry = find(x), find(y)
        if rx != ry:
            parent[rx] = ry

    for i in range(n):
        for j in range(i + 1, n):
            if find(i) == find(j):
                continue
            if strict_overlap(reads[i], indexes[j], indexes_rc[j], k, min_anchors, min_span):
                union(i, j)

    groups = defaultdict(int)
    for i in range(n):
        groups[find(i)] += 1
    return sorted(groups.values(), reverse=True), n


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--r2", required=True)
    ap.add_argument("--n", type=int, default=100000000, help="barcodes to scan (uncapped by default)")
    ap.add_argument("--k", type=int, default=17)
    ap.add_argument("--min-anchors", type=int, default=5)
    ap.add_argument("--min-span", type=int, default=150,
                     help="much longer than assembly's own overlap minimum -- "
                          "see module docstring for why")
    ap.add_argument("--max-reads-per-barcode", type=int, default=40,
                     help="cap for O(reads^2) cost; a subsample is enough to "
                          "detect substantial mixing")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    groups = read_first_barcode_groups(args.r2, args.n)

    import csv
    n_suspect = 0
    n_total = 0
    with open(args.out, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["barcode", "n_reads_checked", "n_components",
                          "largest_component_frac", "suspect_mixed_source"])
        for barcode, reads in groups.items():
            reads = list(dict.fromkeys(reads))
            if len(reads) < 2:
                continue
            n_total += 1
            components, n_checked = connected_components(
                reads, args.k, args.min_anchors, args.min_span, args.max_reads_per_barcode)
            largest_frac = components[0] / n_checked if n_checked else 1.0
            suspect = int(len(components) >= 2 and largest_frac < 0.9)
            n_suspect += suspect
            writer.writerow([barcode, n_checked, len(components),
                              f"{largest_frac:.3f}", suspect])

    print(f"barcodes_checked={n_total}")
    if n_total:
        print(f"suspect_mixed_source={n_suspect} ({100*n_suspect/n_total:.2f}%)")


if __name__ == "__main__":
    main()