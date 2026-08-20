#!/usr/bin/env python3
"""Train the candidate-selection chimera classifier used SHADOW-ONLY (see
denovo.md sec 118/120) by denovo_shadow_score.py.

Predicts P(confident_chimera) for one assembled candidate contig (one of a
UMI's k41_0..k41_N OLC candidates) from 7 QC/geometry features computed by
denovo_junction_qc.py --all-candidates. This is a DIFFERENT model from
models/model_identity.lgb (per-read identity, drives the read-filter
conflict-graph tie-break) -- this one scores whole assembled candidate
contigs, and its production role is presently shadow-only (score, never
decide): denovo.md sec 106-120 found this model's chimera-vs-clean ranking
genuinely generalizes (tail_raw within-UMI pairwise accuracy 0.8507 vs
0.7584 for span_cov_ratio alone), but every active-selection deployment
shape tried (argmin-in-gate, length/support guardrails, veto) either
reproduced or failed to fix a severe-loss regression too large to ship (sec
116/117/119), and training/eval was Zymo-only throughout, so there is
currently no way to verify it on a real sample (sec 118). Do not wire this
model's score into a selection decision without re-reading sec 118-120 first.

Labels come from the project's quarter-split ground truth: a candidate
contig is confident_clean/confident_chimera/undecidable depending on how its
four quarters recruit hits against a known reference. This only exists for
the ZymoBIOMICS control -- there is no way to label a field sample this way,
which is exactly the domain-shift gap shadow deployment exists to probe.

Feature contract (7 columns, order matters for FEATURES below but the model
reads them by name via pandas, unlike model_identity.lgb's positional
contract):
    span_cov_ratio, min_local_span_ratio, placed_reads, contig_len,
    len_ratio, k41_rank, n_candidates
len_ratio and n_candidates are derived per barcode by this script; the rest
come straight from --candidate-qc.

GroupKFold(5, by barcode) CV is reported for reference only, matching
denovo_train_model.py's convention -- CV metrics are not a shipping
criterion by themselves (sec 50's identity-model counterexample applies
here too). The real evidence is the tail_raw one-time final eval already on
record (sec 112) for this exact recipe (n_estimators=300, max_depth=5,
learning_rate=0.05, num_leaves=31, min_child_samples=30) -- retraining with
a different recipe invalidates that evidence until re-verified.
"""
import argparse
import csv
import sys
from collections import defaultdict

FEATURES = ["span_cov_ratio", "min_local_span_ratio", "placed_reads",
            "contig_len", "len_ratio", "k41_rank", "n_candidates"]


def load_candidate_qc(path):
    rows = []
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            rows.append(row)
    return rows


def load_labels(path):
    """contig header -> label, from a quarter-split ground-truth TSV
    (columns: contig, label)."""
    labels = {}
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            labels[row["contig"]] = row["label"]
    return labels


def add_derived_features(rows):
    """Adds len_ratio (contig_len / that barcode's longest candidate) and
    n_candidates (candidate count for that barcode), in place."""
    by_barcode = defaultdict(list)
    for row in rows:
        by_barcode[row["barcode"]].append(row)
    for group in by_barcode.values():
        max_len = max(int(r["contig_len"]) for r in group)
        n = len(group)
        for r in group:
            r["len_ratio"] = (int(r["contig_len"]) / max_len) if max_len else 0.0
            r["n_candidates"] = n
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--candidate-qc", required=True, action="append",
                     help="denovo_junction_qc.py --all-candidates output "
                          "(barcode/header/k41_rank/contig_len/placed_reads/"
                          "span_cov_ratio/min_local_span_ratio); repeatable "
                          "to pool multiple runs")
    ap.add_argument("--labels", required=True, action="append",
                     help="quarter-split ground truth TSV (columns: contig, "
                          "label); repeatable, matched 1:1 with "
                          "--candidate-qc in the order given")
    ap.add_argument("--out-model", required=True)
    ap.add_argument("--out-merged", help="features+labels TSV, for inspection")
    ap.add_argument("--out-metrics", help="write CV metrics here as TSV")
    ap.add_argument("--min-labeled-candidates", type=int, default=1000,
                     help="refuse to train on fewer decidable "
                          "(confident_clean/confident_chimera) rows than this")
    ap.add_argument("--skip-cv", action="store_true")
    args = ap.parse_args()

    if len(args.candidate_qc) != len(args.labels):
        ap.error("--candidate-qc and --labels must be given the same number "
                  "of times (one label file per candidate-qc file)")

    import numpy as np
    import pandas as pd
    import lightgbm as lgb

    all_rows = []
    for qc_path, label_path in zip(args.candidate_qc, args.labels):
        rows = add_derived_features(load_candidate_qc(qc_path))
        labels = load_labels(label_path)
        for row in rows:
            row["label"] = labels.get(row["header"], "undecidable")
        all_rows.extend(rows)
        n_labeled = sum(1 for r in rows if r["header"] in labels)
        print(f"{qc_path}: {len(rows)} candidates, {n_labeled} labeled",
              file=sys.stderr)

    df = pd.DataFrame(all_rows)
    # A handful of candidates carry an empty QC value (contig too degenerate
    # for denovo_junction_qc.py's analyze() to return a result, e.g. zero
    # placed reads) -- drop rather than crash on the cast, same convention
    # as denovo_train_model.py's identity-label dropna.
    qc_cols = ["span_cov_ratio", "min_local_span_ratio", "placed_reads", "contig_len"]
    for col in qc_cols:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    before = len(df)
    df = df.dropna(subset=qc_cols).reset_index(drop=True)
    if len(df) < before:
        print(f"dropped {before - len(df)} candidates with missing/invalid "
              f"QC values", file=sys.stderr)

    df["len_ratio"] = df["len_ratio"].astype(float)
    for col in ["placed_reads", "k41_rank", "n_candidates"]:
        df[col] = df[col].astype(int)

    decidable = df[df["label"].isin(["confident_clean", "confident_chimera"])].copy()
    decidable["y"] = (decidable["label"] == "confident_chimera").astype(int)
    print(f"decidable candidates: {len(decidable)} "
          f"({int(decidable['y'].sum())} chimera, "
          f"{int((1 - decidable['y']).sum())} clean)", file=sys.stderr)
    if len(decidable) < args.min_labeled_candidates:
        sys.exit(f"only {len(decidable)} decidable candidates "
                 f"(< --min-labeled-candidates {args.min_labeled_candidates}); "
                 f"refusing to train")

    if args.out_merged:
        df.to_csv(args.out_merged, sep="\t", index=False)

    X, y, groups = decidable[FEATURES], decidable["y"], decidable["barcode"]

    cv_acc = float("nan")
    if not args.skip_cv:
        from sklearn.model_selection import GroupKFold
        accs = []
        for tr, te in GroupKFold(n_splits=5).split(X, y, groups):
            m = lgb.LGBMClassifier(n_estimators=300, max_depth=5,
                                    learning_rate=0.05, num_leaves=31,
                                    min_child_samples=30, verbose=-1)
            m.fit(X.iloc[tr], y.iloc[tr])
            accs.append(m.score(X.iloc[te], y.iloc[te]))
        cv_acc = float(np.mean(accs))
        print(f"5-fold GroupKFold accuracy: {cv_acc:.4f}", file=sys.stderr)

    model = lgb.LGBMClassifier(n_estimators=300, max_depth=5, learning_rate=0.05,
                                num_leaves=31, min_child_samples=30, verbose=-1)
    model.fit(X, y)
    model.booster_.save_model(args.out_model)
    print(f"\nmodel written: {args.out_model}", file=sys.stderr)

    imp = sorted(zip(FEATURES, model.feature_importances_), key=lambda t: -t[1])
    print("feature importance:", file=sys.stderr)
    for name, val in imp:
        print(f"  {name:<24} {val}", file=sys.stderr)

    if args.out_metrics:
        with open(args.out_metrics, "w") as fh:
            fh.write("metric\tvalue\n")
            fh.write(f"decidable_candidates\t{len(decidable)}\n")
            fh.write(f"chimera_frac\t{decidable['y'].mean():.4f}\n")
            fh.write(f"cv_accuracy\t{cv_acc:.4f}\n")
            for name, val in imp:
                fh.write(f"importance_{name}\t{val}\n")
        print(f"metrics written: {args.out_metrics}", file=sys.stderr)


if __name__ == "__main__":
    main()
