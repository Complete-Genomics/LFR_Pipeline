#!/usr/bin/env python3
"""Feature-distribution drift between a model's training data and new data.

This is the TRIGGER half of auto_retrain.smk, and it is deliberately not the
gate. Drift says the input distribution moved; it does not say a retrained
model would be better -- a model retrained on drifted data can easily be
worse (small control, overfitting, or the "drift" was just noise). Deciding
whether to actually deploy a candidate is denovo_promotion_gate.py's job, and
that one needs real identity measurements. What drift buys is only cost:
skipping a retrain+A/B cycle that would very likely end in "keep incumbent"
anyway. If it were removed entirely, the promotion gate would still keep the
pipeline safe -- it would just burn compute more often.

Metric is PSI (population stability index), the standard drift measure for
tabular features:

    PSI = sum over bins of (new_frac - base_frac) * ln(new_frac / base_frac)

with the conventional reading: < 0.10 no meaningful shift, 0.10-0.25
moderate, > 0.25 significant. PSI rather than a KS test on purpose -- at
these row counts (the shipped baseline is 221k reads) a KS test rejects on
differences far too small to matter, which is the usual large-n trap: it
answers "is there ANY difference" when the question is "is the difference big
enough to act on". PSI measures effect size against fixed thresholds instead.

Three verdicts, because "the distribution moved" is not one situation:

    no_drift      fewer than --min-flagged features moved materially.
    drift         features moved, and NOT in the direction of worse data --
                  a protocol/chemistry change. This is what retraining is for.
    degradation   features moved and the quality-directional ones moved the
                  WRONG way (see QUALITY_DIRECTION). Retraining is refused:
                  the model would learn to treat the decay as normal, hiding
                  a wet-lab problem inside the product. denovo.md sec 61 has
                  the real instance -- one sample in an otherwise fine batch
                  at 52.5% 16S content.

Note what this cannot do: one control cannot separate a bad RUN from a
permanently changed process, since both look identical in a single snapshot.
Degradation is transient and drift is the new normal, so telling them apart
needs several consecutive controls -- history this tool does not keep. The
directional check above is the single-sample approximation: it is what keeps
an unattended pipeline from retraining on decay, without claiming to know
whether the decay is permanent.

Two modes:
    --build-baseline   features TSV -> baseline profile (ships beside the model)
    (default)          baseline profile + new features TSV -> PSI report
"""
import argparse
import json
import math
import sys

# Must match denovo_read_features.py FIELDS (minus read_id/barcode) and
# denovo_train_model.py FEATURES -- the model's positional feature contract.
FEATURES = ["length", "qual_mean", "qual_min", "qual_head", "qual_tail",
            "qual_trend", "hp_frac", "hp_max_run", "mini_gap_mean",
            "mini_gap_max", "mini_gap_var", "pool_size",
            "pool_kmer_popular_frac"]

# Which direction of change means "the data got worse", for the features
# where that question has an answer. PSI itself is symmetric and cannot tell
# a quality collapse from a protocol change -- both are just "the
# distribution moved" -- so direction has to be tracked separately.
#
# This is what separates degradation from drift on a single control, and it
# matters because the two want opposite responses: a lab whose libraries are
# decaying should be told to fix the wet lab, NOT have the model quietly
# retrained to accept the decay as normal (that buries the problem in the
# product and the customer never finds out). A genuine protocol change --
# different read length, different chemistry -- moves features without
# degrading them, and that is the case retraining is for.
#
# -1 means lower is worse, +1 means higher is worse. Features absent here
# (pool_size, mini_gap_*) shift for reasons that are not quality-directional.
QUALITY_DIRECTION = {
    "qual_mean": -1, "qual_min": -1, "qual_head": -1, "qual_tail": -1,
    "length": -1,
    "hp_frac": +1, "hp_max_run": +1,
}

N_BINS = 10
# Guards ln(0) when a bin the baseline covered is empty in the new data (or
# vice versa). Small enough not to distort real bins, large enough to keep a
# single empty bin from dominating the sum.
EPS = 1e-6


def read_feature_columns(path):
    """path -> {feature: [values]}. Streams; only the model's columns are kept."""
    cols = {f: [] for f in FEATURES}
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        try:
            idx = {f: header.index(f) for f in FEATURES}
        except ValueError as exc:
            sys.exit(f"{path}: missing expected feature column ({exc})")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < len(header):
                continue
            for name, i in idx.items():
                try:
                    cols[name].append(float(f[i]))
                except ValueError:
                    pass
    return cols


def quantile_edges(values, n_bins=N_BINS):
    """Bin edges at equal-count quantiles of the baseline.

    Equal-count (not equal-width) so every baseline bin carries real mass;
    equal-width bins on skewed features (pool_size, mini_gap_max) would leave
    most bins near-empty and make PSI hostage to the eps floor. Duplicate
    edges are collapsed, so a feature with few distinct values simply gets
    fewer bins rather than degenerate ones.
    """
    s = sorted(values)
    if not s:
        return []
    edges = [s[int(round(q * (len(s) - 1)))] for q in
             (i / n_bins for i in range(1, n_bins))]
    return sorted(set(edges))


def bin_fractions(values, edges):
    """Fraction of `values` falling in each of len(edges)+1 bins."""
    counts = [0] * (len(edges) + 1)
    for v in values:
        lo, hi = 0, len(edges)
        while lo < hi:  # bisect_right
            mid = (lo + hi) // 2
            if v < edges[mid]:
                hi = mid
            else:
                lo = mid + 1
        counts[lo] += 1
    total = len(values) or 1
    return [c / total for c in counts]


def psi(base_fracs, new_fracs):
    total = 0.0
    for b, n in zip(base_fracs, new_fracs):
        b = max(b, EPS)
        n = max(n, EPS)
        total += (n - b) * math.log(n / b)
    return total


def build_baseline(features_path, out_path):
    cols = read_feature_columns(features_path)
    n_rows = len(cols[FEATURES[0]])
    profile = {"n_rows": n_rows, "n_bins": N_BINS, "features": {}}
    for name in FEATURES:
        values = cols[name]
        edges = quantile_edges(values)
        profile["features"][name] = {
            "edges": edges,
            "fracs": bin_fractions(values, edges),
            "mean": sum(values) / len(values) if values else 0.0,
        }
    with open(out_path, "w") as fh:
        json.dump(profile, fh, indent=1)
    print(f"baseline written: {out_path} ({n_rows} rows, {len(FEATURES)} features)",
          file=sys.stderr)


def relative_worsening(name, base_mean, new_mean):
    """Fractional move in the "worse" direction, or None if not applicable.

    Normalised by the baseline mean so features on different scales (Phred
    ~34, read length ~450, hp_frac ~0.17) share one threshold.
    """
    direction = QUALITY_DIRECTION.get(name)
    if direction is None or not base_mean:
        return None
    return direction * (new_mean - base_mean) / abs(base_mean)


def compare(baseline_path, features_path, out_path, threshold, min_flagged,
            degradation_tolerance=0.05):
    with open(baseline_path) as fh:
        profile = json.load(fh)
    cols = read_feature_columns(features_path)
    n_rows = len(cols[FEATURES[0]])

    rows = []
    flagged = 0
    degraded = []
    for name in FEATURES:
        base = profile["features"].get(name)
        if base is None:
            sys.exit(f"baseline {baseline_path} has no feature {name!r}")
        values = cols[name]
        new_fracs = bin_fractions(values, base["edges"])
        score = psi(base["fracs"], new_fracs)
        new_mean = sum(values) / len(values) if values else 0.0
        is_flagged = score > threshold
        flagged += is_flagged
        worse = relative_worsening(name, base["mean"], new_mean)
        # Only a feature that BOTH moved materially and moved the wrong way
        # counts as degradation; a big PSI on a neutral feature (pool_size,
        # mini_gap_*) is exactly the protocol-change case retraining is for.
        if is_flagged and worse is not None and worse > degradation_tolerance:
            degraded.append(name)
        rows.append((name, score, base["mean"], new_mean, is_flagged, worse))

    rows.sort(key=lambda r: -r[1])
    shifted = flagged >= min_flagged

    if not shifted:
        verdict = "no_drift"
    elif degraded:
        # Distribution moved, but in the direction of worse data. Retraining
        # here would teach the model that the decay is normal; the actionable
        # finding is about the wet lab instead.
        verdict = "degradation"
    else:
        verdict = "drift"

    with open(out_path, "w") as fh:
        fh.write("feature\tpsi\tbaseline_mean\tnew_mean\tflagged\trel_worsening\n")
        for name, score, bmean, nmean, is_flagged, worse in rows:
            w = "NA" if worse is None else f"{worse:+.4f}"
            fh.write(f"{name}\t{score:.4f}\t{bmean:.4f}\t{nmean:.4f}\t"
                     f"{is_flagged}\t{w}\n")
        fh.write(f"#baseline_rows\t{profile['n_rows']}\n")
        fh.write(f"#new_rows\t{n_rows}\n")
        fh.write(f"#psi_threshold\t{threshold}\n")
        fh.write(f"#features_flagged\t{flagged}\n")
        fh.write(f"#min_flagged_for_drift\t{min_flagged}\n")
        fh.write(f"#degradation_tolerance\t{degradation_tolerance}\n")
        fh.write(f"#degraded_features\t{','.join(degraded) if degraded else 'none'}\n")
        fh.write(f"#verdict\t{verdict}\n")

    for name, score, bmean, nmean, is_flagged, worse in rows:
        mark = "  <-- flagged" if is_flagged else ""
        if name in degraded:
            mark = "  <-- WORSE"
        w = "" if worse is None else f"  worse={worse:+.1%}"
        print(f"  {name:<24} PSI={score:7.4f}  base_mean={bmean:10.3f}  "
              f"new_mean={nmean:10.3f}{w}{mark}", file=sys.stderr)
    print(f"\n{flagged}/{len(FEATURES)} features above PSI {threshold}", file=sys.stderr)
    if degraded:
        print(f"degraded (moved the wrong way): {', '.join(degraded)}", file=sys.stderr)
    print(f"verdict: {verdict.upper()}", file=sys.stderr)
    if verdict == "degradation":
        print("  -> NOT a retrain trigger. The control looks worse than the "
              "model's training data, not merely different; retraining would "
              "normalise the decay. Investigate the library prep/run.",
              file=sys.stderr)
    print(f"report written: {out_path}", file=sys.stderr)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--features", required=True,
                    help="denovo_read_features.py output TSV")
    ap.add_argument("--build-baseline", metavar="OUT_JSON",
                    help="build a baseline profile from --features instead of "
                         "comparing against one")
    ap.add_argument("--baseline", help="baseline profile JSON (compare mode)")
    ap.add_argument("--out", help="PSI report TSV (compare mode)")
    ap.add_argument("--psi-threshold", type=float, default=0.25,
                    help="per-feature PSI above which a feature counts as "
                         "shifted; 0.25 is the conventional 'significant' line "
                         "[0.25]")
    ap.add_argument("--min-flagged", type=int, default=2,
                    help="how many shifted features constitute drift. >1 by "
                         "default so a single noisy feature does not trigger a "
                         "retrain cycle on its own [2]")
    ap.add_argument("--degradation-tolerance", type=float, default=0.05,
                    help="a shifted quality feature this much worse than "
                         "baseline (as a fraction of the baseline mean) makes "
                         "the verdict `degradation` rather than `drift` [0.05]")
    args = ap.parse_args()

    if args.build_baseline:
        build_baseline(args.features, args.build_baseline)
        return
    if not args.baseline or not args.out:
        ap.error("compare mode needs --baseline and --out")
    compare(args.baseline, args.features, args.out,
            args.psi_threshold, args.min_flagged, args.degradation_tolerance)


if __name__ == "__main__":
    main()
