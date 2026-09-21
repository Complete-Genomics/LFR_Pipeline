#!/usr/bin/env python3
"""Cheap pre-assembly probe: decide whether the read filter is worth running,
and pick its conflict threshold from this sample's own data.

Two decisions, both made from a small barcode subset so the probe costs
minutes rather than the hour+ the full read filter costs at 1.5-3M UMI:

1. SKIP OR RUN. If almost no reads get dropped on the probe subset, the
   library has no meaningful within-barcode mixing and the full filter is
   wasted compute. Skipping is safe in the sense that it only ever means
   "keep all reads" -- it cannot lose data.

2. DIAGNOSTICS, not auto-tuning. The probe reports this sample's
   cross-barcode identity distribution (reads from different barcodes are
   different source molecules by construction, so this is a direct
   measurement of what "different molecule" looks like here). It does NOT
   set the conflict threshold from it.

   Deriving the threshold automatically was implemented and then rejected on
   evidence (denovo.md sec 33): placing the cutoff at the cross-barcode 99th
   percentile picked 0.95 for the ZymoBIOMICS control, and running the full
   filter at 0.95 measurably degraded the assembly versus the hand-validated
   0.90 (mean identity 94.42 vs 95.04). On the one sample where ground truth
   exists, the "principled" derivation chose the worse value, so the
   threshold stays a validated constant and the distribution is reported for
   a human to look at instead.

Barcode sampling deliberately skips low-complexity barcodes (AAAA...-style
artifacts) and does not read from the head of the sorted TSV, because sorted
order concentrates those artifacts at the front and biased the chimera rate
badly in an earlier analysis (denovo.md sec 29).
"""
import argparse
import itertools
import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + "/test")
from analyze_olc_shortfall import kmer_index  # noqa: E402
from denovo_read_filter import (iter_barcode_groups, overlap_identity,  # noqa: E402
                                 find_contaminants)


def is_low_complexity(barcode, max_frac=0.8):
    if not barcode:
        return True
    return max(barcode.count(c) for c in "ACGT") / len(barcode) >= max_frac


def collect_subset(path, n_barcodes, min_reads, rng, reservoir_factor=20):
    """Reservoir-sample barcodes across the WHOLE file, skipping low-complexity
    ones, so the probe is not dominated by the artifact barcodes that sort to
    the front."""
    keep = []
    seen = 0
    target = n_barcodes * reservoir_factor
    for barcode, _lines, seqs in iter_barcode_groups(path):
        if is_low_complexity(barcode) or len(seqs) < min_reads:
            continue
        seen += 1
        if len(keep) < target:
            keep.append((barcode, seqs))
        else:
            j = rng.randrange(seen)
            if j < target:
                keep[j] = (barcode, seqs)
    rng.shuffle(keep)
    return keep[:n_barcodes]


def cross_barcode_identities(subset, n_pairs, min_overlap, rng):
    """Identity distribution for reads from DIFFERENT barcodes == known
    different source molecules == this sample's conflict reference."""
    idents = []
    attempts = 0
    while len(idents) < n_pairs and attempts < n_pairs * 40:
        attempts += 1
        (b1, r1), (b2, r2) = rng.sample(subset, 2)
        a, b = rng.choice(r1), rng.choice(r2)
        if len(a) < min_overlap or len(b) < min_overlap:
            continue
        v = overlap_identity(a, kmer_index(b, 17), b, min_overlap)
        if v is not None:
            idents.append(v)
    return idents


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--r2", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--n-barcodes", type=int, default=400)
    ap.add_argument("--min-reads", type=int, default=8)
    ap.add_argument("--pairs", type=int, default=3000)
    ap.add_argument("--min-overlap", type=int, default=150)
    ap.add_argument("--max-reads", type=int, default=25)
    ap.add_argument("--skip-below-drop-rate", type=float, default=0.02,
                     help="if the probe drops less than this fraction of reads, "
                          "recommend skipping the full read filter")
    ap.add_argument("--same-molecule-id", type=float, default=0.90,
                     help="conflict threshold used for the probe and reported "
                          "for the full run. Validated constant, not derived "
                          "from the data -- see module docstring")
    ap.add_argument("--cross-percentile", type=float, default=0.99,
                     help="cross-barcode percentile reported as a diagnostic")
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    rng = random.Random(args.seed)
    subset = collect_subset(args.r2, args.n_barcodes, args.min_reads, rng)
    if len(subset) < 20:
        with open(args.out, "w") as fh:
            fh.write("metric\tvalue\n")
            fh.write("status\tinsufficient_data\n")
            fh.write("run_read_filter\tTrue\n")
            fh.write(f"same_molecule_id\t{args.same_molecule_id}\n")
        print("probe: too few usable barcodes; defaulting to run filter at 0.90")
        return

    cross = cross_barcode_identities(subset, args.pairs, args.min_overlap, rng)
    if cross:
        cross_sorted = sorted(cross)
        idx = min(len(cross_sorted) - 1,
                   int(args.cross_percentile * len(cross_sorted)))
        derived = cross_sorted[idx]
        cross_median = cross_sorted[len(cross_sorted) // 2]
    else:
        derived, cross_median = float("nan"), float("nan")
    threshold = args.same_molecule_id

    dropped = checked = 0
    for _bc, seqs in subset:
        drop = find_contaminants(seqs, threshold, args.min_overlap, args.max_reads)
        dropped += len(drop)
        # denominator is ALL reads in the barcode, not just the examined cap,
        # so this rate is directly comparable to denovo_read_filter's own
        # reads_dropped percentage (reads past --max-reads are always kept)
        checked += len(seqs)
    drop_rate = dropped / checked if checked else 0.0
    run_filter = drop_rate >= args.skip_below_drop_rate

    with open(args.out, "w") as fh:
        fh.write("metric\tvalue\n")
        fh.write("status\tok\n")
        fh.write(f"probe_barcodes\t{len(subset)}\n")
        fh.write(f"cross_barcode_pairs\t{len(cross)}\n")
        fh.write(f"cross_median_identity\t{cross_median:.4f}\n")
        fh.write(f"cross_p{int(args.cross_percentile*100)}_identity\t{derived:.4f}\n")
        fh.write(f"same_molecule_id\t{threshold:.4f}\n")
        fh.write(f"probe_drop_rate\t{drop_rate:.4f}\n")
        fh.write(f"run_read_filter\t{run_filter}\n")

    print(f"probe barcodes={len(subset)}  cross-barcode median={cross_median:.3f} "
          f"(diagnostic only, p{int(args.cross_percentile*100)}={derived:.3f})")
    print(f"conflict threshold={threshold:.3f} (fixed, validated -- not derived)")
    # Cross-barcode median identity separates sample regimes, and the regimes
    # want opposite settings (denovo.md sec 35). Reads from different barcodes
    # are different organisms by construction, so a high median means the
    # community is dominated by closely related organisms (mock, gut), where
    # co-barcoded reads CAN merge into chimeras and filtering pays; a low
    # median means a diverse community (soil), where they cannot merge and
    # filtering mostly costs contig length.
    if cross_median == cross_median:  # not NaN
        if cross_median >= 0.70:
            suggestion = "balanced"
            why = "low-diversity community (mock/gut-like); chimeras do form here"
        elif cross_median <= 0.60:
            suggestion = "high_div"
            why = "high-diversity community (soil-like); filtering costs more than it fixes"
        else:
            suggestion = "balanced"
            why = "intermediate diversity; balanced is the safer default"
        print(f"suggested qc_preset={suggestion}  ({why})")
        with open(args.out, "a") as fh:
            fh.write(f"suggested_qc_preset\t{suggestion}\n")

    print(f"probe drop rate={100*drop_rate:.2f}%  -> run_read_filter={run_filter}")


if __name__ == "__main__":
    main()