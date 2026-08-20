#!/usr/bin/env python3
"""Shadow-deployment scorer for the candidate-selection chimera model.

Scores every assembled candidate and compares the rule-delivered candidate
with the model's lowest-risk candidate. Nothing this script computes is
allowed to change what gets delivered -- it is deliberately a monitoring-only
dead end, not a step in the selection path. No downstream rule reads its
output.

Why this exists (sec 118): the model's chimera-ranking genuinely
generalizes (sec 112, tail_raw within-UMI accuracy 0.85 vs 0.76 for the
single best rule feature), but every active-selection shape tried produced
a severe-loss regression too large to ship (sec 116/117/119), and
training/eval so far is Zymo-only, so there is no way yet to check this
model against ground truth on a real sample. Shadow scoring is what closes
that gap: run it against real batches, then look at whether the score
distribution and the disagreement rate with candidate_select's actual
choice resemble what Zymo predicted, before ever reconsidering an active
(veto or otherwise) role for this model.

Consumes:
  --candidate-qc  denovo_junction_qc.py --all-candidates output (every
                   candidate's QC -- needed to derive len_ratio/n_candidates
                   for the delivered candidate specifically, since those are
                   relative to that barcode's whole candidate set)
  --selection     denovo_candidate_select.py --out-report (which k41_rank
                   got delivered per barcode)
  --model         a model produced by denovo_train_candidate_model.py

Writes a candidate-level table and one per-UMI summary table.
"""
import argparse
import csv
import os
from collections import defaultdict

FEATURES = ["span_cov_ratio", "min_local_span_ratio", "placed_reads",
            "contig_len", "len_ratio", "k41_rank", "n_candidates"]


def load_candidate_qc(path):
    by_barcode = defaultdict(dict)
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            by_barcode[row["barcode"]][int(row["k41_rank"])] = row
    return by_barcode


def load_selection(path):
    chosen = {}
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            chosen[row["barcode"]] = (int(row["chosen_rank"]), int(row["switched"]))
    return chosen


def build_features(barcode, rank, candidates):
    group = candidates.get(barcode)
    if not group or rank not in group:
        return None
    row = group[rank]
    max_len = max(int(r["contig_len"]) for r in group.values())
    return {
        "span_cov_ratio": float(row["span_cov_ratio"]),
        "min_local_span_ratio": float(row["min_local_span_ratio"]),
        "placed_reads": int(row["placed_reads"]),
        "contig_len": int(row["contig_len"]),
        "len_ratio": (int(row["contig_len"]) / max_len) if max_len else 0.0,
        "k41_rank": rank,
        "n_candidates": len(group),
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--candidate-qc", required=True)
    ap.add_argument("--selection", required=True)
    ap.add_argument("--model", required=True)
    ap.add_argument("--out-candidates", required=True)
    ap.add_argument("--out-summary", required=True)
    args = ap.parse_args()

    import lightgbm as lgb
    import numpy as np

    candidates = load_candidate_qc(args.candidate_qc)
    selection = load_selection(args.selection)
    booster = lgb.Booster(model_file=args.model)

    rows = []
    feats = []
    for barcode, group in candidates.items():
        chosen_rank = selection.get(barcode, (None, 0))[0]
        for rank in sorted(group):
            f = build_features(barcode, rank, candidates)
            if f is None:
                continue
            rows.append({"barcode": barcode, "k41_rank": rank,
                         "selected_by_rule": int(rank == chosen_rank), **f})
            feats.append([f[c] for c in FEATURES])

    scores = booster.predict(np.array(feats)) if feats else []
    for row, score in zip(rows, scores):
        row["p_chimera"] = float(score)

    os.makedirs(os.path.dirname(args.out_candidates) or ".", exist_ok=True)
    os.makedirs(os.path.dirname(args.out_summary) or ".", exist_ok=True)
    candidate_fields = ["barcode", "k41_rank", "selected_by_rule",
                        "p_chimera"] + [
                            feature for feature in FEATURES if feature != "k41_rank"
                        ]
    with open(args.out_candidates, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=candidate_fields, delimiter="\t")
        writer.writeheader()
        for row in rows:
            formatted = dict(row)
            formatted["p_chimera"] = f"{row['p_chimera']:.6f}"
            writer.writerow(formatted)

    by_barcode = defaultdict(list)
    for row in rows:
        by_barcode[row["barcode"]].append(row)

    summary_fields = ["barcode", "delivered_rank", "switched",
                      "delivered_p_chimera", "lowest_p_chimera_rank",
                      "lowest_p_chimera", "gbdt_disagrees_with_rule",
                      "n_candidates", "missing_delivered_qc"]
    n_missing = 0
    with open(args.out_summary, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=summary_fields, delimiter="\t")
        writer.writeheader()
        for barcode, (rank, switched) in sorted(selection.items()):
            group = by_barcode.get(barcode, [])
            delivered = next((row for row in group if row["k41_rank"] == rank), None)
            if delivered is None:
                n_missing += 1
                writer.writerow({"barcode": barcode, "delivered_rank": rank,
                                 "switched": switched, "n_candidates": len(group),
                                 "missing_delivered_qc": 1})
                continue
            lowest = min(group, key=lambda row: (row["p_chimera"], row["k41_rank"]))
            writer.writerow({
                "barcode": barcode,
                "delivered_rank": rank,
                "switched": switched,
                "delivered_p_chimera": f"{delivered['p_chimera']:.6f}",
                "lowest_p_chimera_rank": lowest["k41_rank"],
                "lowest_p_chimera": f"{lowest['p_chimera']:.6f}",
                "gbdt_disagrees_with_rule": int(lowest["k41_rank"] != rank),
                "n_candidates": len(group),
                "missing_delivered_qc": 0,
            })

    print(f"scored_candidates={len(rows)} scored_barcodes={len(selection) - n_missing} "
          f"missing_delivered_qc={n_missing}")


if __name__ == "__main__":
    main()
