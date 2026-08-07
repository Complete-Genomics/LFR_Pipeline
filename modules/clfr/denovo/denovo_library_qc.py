#!/usr/bin/env python3
"""Pre-assembly library-level mixing QC (reference-free, with built-in control).

Question this answers: before spending hours assembling, does this library
have a read-pool mixing problem at all (barcodes carrying DNA from more than
one source molecule)? It reports one library-level rate, not a per-barcode
verdict -- denovo.md sec 30 showed per-barcode calls from this kind of data
are dominated by depth/coverage noise, but the aggregate rate is stable.

The trick is that the negative control is inside the data. Take pairs of
reads that overlap each other:
  - WITHIN one barcode: a clean library gives pairs off the same source
    molecule, so overlap identity should sit near the sequencing-error floor
    (~99%).
  - ACROSS different barcodes: these are, by construction, different source
    molecules. Any overlap between them is coincidental similarity (16S
    conserved regions), so their identity distribution is exactly what
    "different molecules that happen to overlap" looks like in THIS library.

If the library is clean, the within-barcode distribution is almost entirely
in the high-identity mode. The more the within-barcode distribution borrows
mass from the cross-barcode mode, the more mixing there is. No reference
database is involved, and both distributions come from the same run, so
chemistry/error-rate/species-composition effects cancel out.
"""
import argparse
import os
import random
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + "/test")
from analyze_olc_shortfall import read_first_barcode_groups, kmer_index, place, rc  # noqa: E402


def overlap_identity(read_a, read_b, k=17, min_overlap=100):
    """Identity over the overlapping stretch of two reads, or None if they do
    not overlap unambiguously. Ungapped, matching the placement model used
    elsewhere in this module."""
    index = kmer_index(read_b, k)
    best = None
    for oriented in (read_a, rc(read_a)):
        hit = place(oriented, read_b, index)
        if hit is None:
            continue
        start = hit[0]
        lo = max(0, start)
        hi = min(len(read_b), start + len(oriented))
        if hi - lo < min_overlap:
            continue
        matches = sum(1 for pos in range(lo, hi)
                      if oriented[pos - start] == read_b[pos])
        ident = matches / (hi - lo)
        if best is None or ident > best:
            best = ident
    return best


def sample_pairs(groups, n_pairs, within, rng, max_reads, min_overlap):
    """Collect overlap identities for within-barcode or cross-barcode pairs."""
    barcodes = list(groups.keys())
    idents = []
    attempts = 0
    max_attempts = n_pairs * 40
    while len(idents) < n_pairs and attempts < max_attempts:
        attempts += 1
        if within:
            bc = rng.choice(barcodes)
            reads = list(dict.fromkeys(groups[bc]))[:max_reads]
            if len(reads) < 2:
                continue
            a, b = rng.sample(reads, 2)
        else:
            bc1, bc2 = rng.sample(barcodes, 2)
            r1 = list(dict.fromkeys(groups[bc1]))[:max_reads]
            r2 = list(dict.fromkeys(groups[bc2]))[:max_reads]
            if not r1 or not r2:
                continue
            a, b = rng.choice(r1), rng.choice(r2)
        if len(a) < min_overlap or len(b) < min_overlap:
            continue
        ident = overlap_identity(a, b, min_overlap=min_overlap)
        if ident is not None:
            idents.append(ident)
    return idents


def summarise(idents, label, cut):
    n = len(idents)
    if n == 0:
        print(f"{label}: no overlapping pairs found")
        return None
    hi = sum(1 for x in idents if x >= cut)
    print(f"{label}: n={n}  median={sorted(idents)[n//2]:.4f}  "
          f"frac>={cut:.2f}: {hi/n:.4f}")
    return hi / n


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--r2", required=True)
    ap.add_argument("--n-barcodes", type=int, default=3000,
                     help="barcodes to sample the library from")
    ap.add_argument("--pairs", type=int, default=3000,
                     help="overlapping read pairs to collect per group")
    ap.add_argument("--max-reads", type=int, default=40)
    ap.add_argument("--min-overlap", type=int, default=100)
    ap.add_argument("--same-molecule-cut", type=float, default=0.97,
                     help="identity above which a pair is called same-molecule")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    rng = random.Random(args.seed)
    groups = read_first_barcode_groups(args.r2, args.n_barcodes)

    within = sample_pairs(groups, args.pairs, True, rng, args.max_reads, args.min_overlap)
    cross = sample_pairs(groups, args.pairs, False, rng, args.max_reads, args.min_overlap)

    w = summarise(within, "within-barcode ", args.same_molecule_cut)
    c = summarise(cross, "cross-barcode  ", args.same_molecule_cut)

    with open(args.out, "w") as fh:
        fh.write("metric\tvalue\n")
        fh.write(f"within_barcode_pairs\t{len(within)}\n")
        fh.write(f"cross_barcode_pairs\t{len(cross)}\n")
        if w is not None and c is not None:
            fh.write(f"within_frac_same_molecule\t{w:.4f}\n")
            fh.write(f"cross_frac_same_molecule\t{c:.4f}\n")
            # Cross-barcode pairs are different source molecules by
            # construction, so their same-molecule-looking rate is this test's
            # false-positive floor. How far the within-barcode rate rises above
            # that floor is the barcoding signal: a clean library concentrates
            # same-molecule pairs inside barcodes, a mixed one does not.
            pseudo = 0.5 / max(1, len(cross))
            enrichment = w / max(c, pseudo)
            fh.write(f"same_molecule_enrichment\t{enrichment:.2f}\n")
            print(f"\nsame-molecule enrichment (within vs cross): {enrichment:.1f}x")
            print("  Higher = cleaner barcoding. This is a RELATIVE indicator for")
            print("  comparing libraries/runs, not an absolute mixing percentage:")
            print("  the cross-barcode floor rises in low-diversity samples (two")
            print("  barcodes often hold the same species), which deflates the")
            print("  ratio for mock communities relative to diverse samples.")
            print("  Measured (denovo.md sec 31): ZymoBIOMICS mock 7.5x (known")
            print("  heavily mixed) vs soil 236x (known far cleaner).")


if __name__ == "__main__":
    main()