#!/usr/bin/env python3
"""Pick which assembled candidate contig (k41_0, k41_1, ...) a UMI delivers.

Two modes:

  longest (default, current production behaviour, unchanged): always deliver
    k41_0 -- denovo_seed_olc.py writes each barcode's candidates
    longest-first, so this is "keep only the longest per UMI", same as the
    awk filter this rule used to be.

  gated_switch (opt-in, denovo.md sec 109-112): deliver the first
    k41_rank-ordered candidate with span_cov_ratio >= --max-span-ratio AND
    placed_reads >= --min-placed-reads, falling back to k41_0 if no
    candidate qualifies. Rule-only, no model -- validated by a one-time
    tail_raw held-out eval (sec 112): decidable-basis chimera rate
    9.64% -> 5.52%.

    Caveat worth reading before turning this on (sec 116/117/120): measured
    WITHOUT any post-selection polish, this rule's severe-loss (>=5pt
    identity drop vs k41_0) was 3.93% on that same tail_raw set -- right at
    the project's 3% acceptance line, and above it on that particular
    sample (merged-pool measured 1.92%, comfortably under; the two samples
    disagree). A racon polish step after selection (sec 100/103) is known to
    help but is NOT wired into this pipeline yet. `denovo/
    candidate_select_report.tsv`'s `switched` column is exactly the
    diagnostic to watch if this is enabled on new data: an unexpectedly high
    switch rate on a real sample is the earliest signal something differs
    from Zymo.

    A model-scored variant (GBDT+gate, argmin P(chimera) among gate
    qualifiers) was also evaluated and its production readiness was
    REVOKED (sec 116: severe-loss 15.94%, 4x this rule's own) -- do not add
    a model to the selection decision itself. The model's only validated
    role right now is shadow scoring (denovo_shadow_score.py, sec 118/120),
    which never touches this choice.

Input contract: --contigs is denovo_seed_olc.py's raw final_contigs_N.fa
(ALL candidates per barcode), not the already-reduced denovo.longest.fasta.
--candidate-qc is denovo_junction_qc.py --all-candidates output (required
for gated_switch, ignored for longest).
"""
import argparse
import csv
import sys
from collections import defaultdict


def load_candidates_fasta(path):
    """barcode -> [(k41_rank, header, seq), ...] in file order.

    k41_rank is assigned by position, not by parsing the header -- this is
    the same order-as-rank assumption denovo_junction_qc.py --all-candidates
    and the historical filterOLC_longest awk filter both already rely on
    (denovo_seed_olc.py guarantees longest-first, i.e. position 0 == k41_0).
    """
    grouped = defaultdict(list)
    header = None
    chunks = []

    def flush():
        if header is None:
            return
        barcode = header[1:].split(">", 1)[0]
        grouped[barcode].append((header, "".join(chunks)))

    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                flush()
                header = line
                chunks = []
            else:
                chunks.append(line)
        flush()

    out = {}
    n_order_mismatch = 0
    for barcode, entries in grouped.items():
        ranked = [(rank, h, s) for rank, (h, s) in enumerate(entries)]
        if ranked and not ranked[0][1].endswith(">k41_0"):
            n_order_mismatch += 1
        out[barcode] = ranked
    if n_order_mismatch:
        print(f"WARNING: {n_order_mismatch} barcode(s) whose first fasta "
              f"entry is not >k41_0 -- the longest-first order this script "
              f"relies on may be violated", file=sys.stderr)
    return out


def load_candidate_qc(path):
    """(barcode, k41_rank) -> QC row dict, from denovo_junction_qc.py
    --all-candidates output."""
    qc = {}
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            qc[(row["barcode"], int(row["k41_rank"]))] = row
    return qc


def select_barcode(barcode, entries, qc, mode, max_span_ratio, min_placed_reads):
    """entries: [(k41_rank, header, seq), ...], k41_rank-ordered (0 = k41_0).
    Returns (chosen_rank, chosen_header, chosen_seq, switched, reason)."""
    primary_rank, primary_header, primary_seq = entries[0]
    if mode == "longest":
        return primary_rank, primary_header, primary_seq, False, "mode=longest"

    for rank, header, seq in entries:
        row = qc.get((barcode, rank))
        if row is None:
            continue  # no QC row (e.g. zero placed reads) -- never a gate pass
        try:
            span_cov_ratio = float(row["span_cov_ratio"])
            placed_reads = int(row["placed_reads"])
        except (TypeError, ValueError):
            continue
        if span_cov_ratio >= max_span_ratio and placed_reads >= min_placed_reads:
            switched = rank != primary_rank
            reason = "gate_pass_switched" if switched else "gate_pass_primary"
            return rank, header, seq, switched, reason

    return primary_rank, primary_header, primary_seq, False, "no_candidate_passed_gate"


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--contigs", required=True,
                     help="final_contigs_0.fa (ALL candidates per barcode, "
                          "k41_0 first)")
    ap.add_argument("--candidate-qc",
                     help="denovo_junction_qc.py --all-candidates output; "
                          "required when --mode gated_switch")
    ap.add_argument("--mode", choices=["longest", "gated_switch"], default="longest")
    ap.add_argument("--max-span-ratio", type=float, default=0.25)
    ap.add_argument("--min-placed-reads", type=int, default=2)
    ap.add_argument("--out-fasta", required=True)
    ap.add_argument("--out-report", required=True,
                     help="per-barcode outcome: which candidate was "
                          "delivered and whether it was a switch")
    ap.add_argument("--out-decision", required=True,
                     help="single-row settings record, mirrors qc_decision.tsv")
    args = ap.parse_args()

    if args.mode == "gated_switch" and not args.candidate_qc:
        ap.error("--mode gated_switch requires --candidate-qc")

    candidates = load_candidates_fasta(args.contigs)
    qc = load_candidate_qc(args.candidate_qc) if args.candidate_qc else {}

    fields = ["barcode", "chosen_rank", "n_candidates", "switched", "reason"]
    n_switched = 0
    n_total = 0
    with open(args.out_fasta, "w") as out_fa, \
         open(args.out_report, "w", newline="") as out_rep:
        writer = csv.DictWriter(out_rep, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for barcode, entries in candidates.items():
            n_total += 1
            rank, header, seq, switched, reason = select_barcode(
                barcode, entries, qc, args.mode,
                args.max_span_ratio, args.min_placed_reads)
            n_switched += int(switched)
            out_fa.write(f"{header}\n{seq}\n")
            writer.writerow({"barcode": barcode, "chosen_rank": rank,
                              "n_candidates": len(entries),
                              "switched": int(switched), "reason": reason})

    with open(args.out_decision, "w") as fh:
        fh.write("setting\tvalue\n")
        fh.write(f"mode\t{args.mode}\n")
        fh.write(f"max_span_ratio\t{args.max_span_ratio}\n")
        fh.write(f"min_placed_reads\t{args.min_placed_reads}\n")
        fh.write(f"barcodes_total\t{n_total}\n")
        fh.write(f"barcodes_switched\t{n_switched}\n")

    print(f"barcodes_total={n_total}")
    if n_total:
        print(f"switched={n_switched} ({100 * n_switched / n_total:.2f}%)")


if __name__ == "__main__":
    main()
