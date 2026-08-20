#!/usr/bin/env python3
"""Paired identity comparison against a "local ground truth" reference.

Implements the validation method used throughout denovo.md (sec 59/62/67) to
compare read-filter arms without a known-truth reference: assemble a NEUTRAL
arm (no read filtering at all) and use vsearch to assign each barcode its own
best-hit reference from a broad 16S database. That per-barcode reference is
then held FIXED while comparing other arms' contigs against it, so the
comparison is paired (same barcode, same target, only the query differs) and
does not depend on the reference being a close match in absolute terms --
only the paired difference matters. Limitation carried over from that same
history: this compresses true effect size 3-5x and cannot rank two
close-performing arms (denovo.md sec 55/59), so treat p-values as directional
evidence, not a precise estimate.

Inputs are vsearch --userout files (query, target, id columns):
  --reference-hits   nofilter arm vs the big reference db, TOP HIT ONLY
                      (assigns each barcode's fixed reference)
  --arm-a-hits        arm A (e.g. plain) vs the small assigned-refs subset,
                      ALL accepted hits (may include targets other than the
                      one assigned -- this script filters to the assigned one)
  --arm-b-hits        arm B (e.g. ml), same shape as arm A
"""
import argparse
import csv
import sys
from collections import defaultdict


def load_top_hit(path):
    """query -> target, from a top-hit-only vsearch userout file."""
    out = {}
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 2:
                continue
            out[f[0]] = f[1]
    return out


def load_identity_to_target(path):
    """query -> {target: identity}, from a (possibly multi-hit) userout file."""
    out = defaultdict(dict)
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 3:
                continue
            query, target, ident = f[0], f[1], float(f[2])
            out[query][target] = ident
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--reference-hits", required=True)
    ap.add_argument("--arm-a-hits", required=True)
    ap.add_argument("--arm-b-hits", required=True)
    ap.add_argument("--arm-a-name", default="arm_a")
    ap.add_argument("--arm-b-name", default="arm_b")
    ap.add_argument("--out", required=True, help="per-barcode paired TSV")
    args = ap.parse_args()

    assigned_ref = load_top_hit(args.reference_hits)
    a_hits = load_identity_to_target(args.arm_a_hits)
    b_hits = load_identity_to_target(args.arm_b_hits)

    rows = []
    for query, ref in assigned_ref.items():
        a_id = a_hits.get(query, {}).get(ref)
        b_id = b_hits.get(query, {}).get(ref)
        if a_id is None or b_id is None:
            continue
        rows.append((query, ref, a_id, b_id))

    with open(args.out, "w", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(["contig_id", "assigned_ref", args.arm_a_name, args.arm_b_name,
                    "diff_b_minus_a"])
        for query, ref, a_id, b_id in rows:
            w.writerow([query, ref, f"{a_id:.2f}", f"{b_id:.2f}", f"{b_id - a_id:.2f}"])

    n = len(rows)
    print(f"barcodes with a reference:        {len(assigned_ref)}", file=sys.stderr)
    print(f"barcodes with BOTH arms scored:    {n}", file=sys.stderr)
    if n == 0:
        print("no paired rows -- nothing to test", file=sys.stderr)
        return

    a_vals = [r[2] for r in rows]
    b_vals = [r[3] for r in rows]
    diffs = [b - a for a, b in zip(a_vals, b_vals)]
    wins = sum(1 for d in diffs if d > 0)
    losses = sum(1 for d in diffs if d < 0)
    ties = n - wins - losses
    print(f"{args.arm_a_name} mean identity:   {sum(a_vals)/n:.4f}", file=sys.stderr)
    print(f"{args.arm_b_name} mean identity:   {sum(b_vals)/n:.4f}", file=sys.stderr)
    print(f"mean diff ({args.arm_b_name} - {args.arm_a_name}): {sum(diffs)/n:.4f}",
          file=sys.stderr)
    print(f"win/loss/tie ({args.arm_b_name} vs {args.arm_a_name}): "
          f"{wins}/{losses}/{ties}", file=sys.stderr)

    try:
        from scipy.stats import wilcoxon
        nonzero = [d for d in diffs if d != 0]
        if len(nonzero) >= 10:
            stat, p = wilcoxon(a_vals, b_vals)
            print(f"paired Wilcoxon signed-rank: statistic={stat:.1f} p={p:.4g}",
                  file=sys.stderr)
        else:
            print("too few non-tied pairs for Wilcoxon", file=sys.stderr)
    except ImportError:
        print("scipy not available -- skipping significance test "
              "(win/loss/tie and mean diff above still stand)", file=sys.stderr)

    print(f"\nwritten: {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
