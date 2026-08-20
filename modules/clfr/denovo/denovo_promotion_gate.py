#!/usr/bin/env python3
"""Decide whether a retrained candidate model replaces the incumbent.

This is the safety property that makes unattended retraining acceptable.
denovo.md sec 50 is the standing counterexample: a model can predict per-read
identity well (MAE 2.10 vs an 8.48 baseline) and still make a WORSE filter,
because what matters is which read the conflict graph drops, not how
accurately each read is scored. So training metrics never decide promotion --
the two models are used to actually filter and assemble the control, and the
resulting contigs are scored against the control's KNOWN reference.

Because a mock community control has real ground truth, this is a true A/B,
not the greengenes "local ground truth" proxy random_inspection.smk falls
back to for field samples (that proxy compresses effect size 3-5x and cannot
rank two close arms -- denovo.md sec 55/59; it would be the wrong instrument
for a promote/reject decision).

Promotion requires the candidate to win by --min-improvement, not merely to
tie or edge ahead: an unattended swap should need evidence, and the incumbent
is the known quantity. It also rejects an excessive severe-loss tail and,
when primary-contig/readback artifacts are provided, any loss of >=1 kb yield
or raw-read-supported bases. Everything else keeps the incumbent and says why.
"""
import argparse
import csv
import sys
from collections import defaultdict


def load_identity(path):
    """query -> {target: identity} from a vsearch --userout (query+target+id)."""
    out = defaultdict(dict)
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 3:
                continue
            out[f[0]][f[1]] = float(f[2])
    return out


def best_identity(hits):
    """query -> best identity over all targets, for a known-truth reference.

    Unlike the field-sample path, the reference here is the control's own
    small mock-community database, so the query's best hit in it IS the
    organism it came from -- no per-barcode reference assignment needed.
    """
    return {q: max(t.values()) for q, t in hits.items() if t}


def barcode_from_contig_id(contig_id):
    """Return the barcode portion of the established ``barcode>k41_0`` id."""
    return contig_id.split(">", 1)[0]


def primary_lengths(path):
    """Return longest contig length per barcode from a FASTA file."""
    lengths = {}
    name = None
    length = 0
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    barcode = barcode_from_contig_id(name)
                    lengths[barcode] = max(lengths.get(barcode, 0), length)
                name = line[1:].split()[0]
                length = 0
            else:
                length += len(line)
    if name is not None:
        barcode = barcode_from_contig_id(name)
        lengths[barcode] = max(lengths.get(barcode, 0), length)
    return lengths


def readback_summary(path):
    """Summarize raw-read support without using reference identity.

    ``supported_bp`` is the breadth-weighted length of all primary contigs.
    It rewards a candidate for extending a contig only when raw reads cover
    that extension, and is intentionally computed from the same unfiltered
    read pool for both A/B arms.
    """
    rows = 0
    two_end = 0
    supported_bp = 0.0
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            rows += 1
            two_end += int(float(row["two_end_supported"]))
            supported_bp += float(row["contig_len"]) * float(row["breadth_1x"])
    return {"contigs_checked": rows, "two_end_supported": two_end,
            "supported_bp": supported_bp}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--incumbent-hits", required=True,
                    help="vsearch userout: incumbent-model arm vs control reference")
    ap.add_argument("--candidate-hits", required=True,
                    help="vsearch userout: candidate-model arm vs control reference")
    ap.add_argument("--incumbent-model", required=True, help="path, recorded in the decision")
    ap.add_argument("--candidate-model", required=True, help="path, recorded in the decision")
    ap.add_argument("--out-decision", required=True)
    ap.add_argument("--out-promoted-model", required=True,
                    help="the model the pipeline should use going forward; this "
                         "script copies whichever side won here")
    ap.add_argument("--min-improvement", type=float, default=0.05,
                    help="mean identity points the candidate must gain to be "
                         "promoted [0.05]")
    ap.add_argument("--min-paired", type=int, default=200,
                    help="minimum barcodes scored by BOTH arms; below this the "
                         "comparison is too small to act on unattended [200]")
    ap.add_argument("--max-alpha", type=float, default=0.05,
                    help="required Wilcoxon p-value; scipy absence refuses promotion [0.05]")
    ap.add_argument("--severe-loss-points", type=float, default=5.0,
                    help="identity-point drop counted as a severe regression [5.0]")
    ap.add_argument("--max-severe-loss-rate", type=float, default=0.01,
                    help="maximum paired severe-regression rate [0.01]")
    ap.add_argument("--incumbent-contigs",
                    help="primary-contig FASTA for the incumbent arm")
    ap.add_argument("--candidate-contigs",
                    help="primary-contig FASTA for the candidate arm")
    ap.add_argument("--incumbent-readback",
                    help="raw-read readback TSV for the incumbent arm")
    ap.add_argument("--candidate-readback",
                    help="raw-read readback TSV for the candidate arm")
    ap.add_argument("--min-primary-len", type=int, default=1000,
                    help="minimum primary-contig length used for yield [1000]")
    ap.add_argument("--min-yield-ratio", type=float, default=1.0,
                    help="candidate/incumbent minimum >=min-primary-len yield [1.0]")
    ap.add_argument("--min-readback-bp-ratio", type=float, default=1.0,
                    help="candidate/incumbent minimum breadth-weighted readback bp [1.0]")
    args = ap.parse_args()

    support_args = [args.incumbent_contigs, args.candidate_contigs,
                    args.incumbent_readback, args.candidate_readback]
    if any(support_args) and not all(support_args):
        ap.error("contig-yield/readback gate requires all four --*-contigs and "
                 "--*-readback inputs")

    inc = best_identity(load_identity(args.incumbent_hits))
    cand = best_identity(load_identity(args.candidate_hits))

    shared = sorted(set(inc) & set(cand))
    n = len(shared)
    lines = []

    def record(key, value):
        lines.append(f"{key}\t{value}")

    record("incumbent_model", args.incumbent_model)
    record("candidate_model", args.candidate_model)
    record("incumbent_contigs_scored", len(inc))
    record("candidate_contigs_scored", len(cand))
    record("paired_contigs", n)

    if n < args.min_paired:
        record("decision", "keep_incumbent")
        record("reason", f"only {n} paired contigs (< min_paired {args.min_paired})")
        promote = False
        p_value = float("nan")
        delta = float("nan")
    else:
        a = [inc[q] for q in shared]
        b = [cand[q] for q in shared]
        mean_a, mean_b = sum(a) / n, sum(b) / n
        delta = mean_b - mean_a
        diffs = [y - x for x, y in zip(a, b)]
        wins = sum(1 for d in diffs if d > 0)
        losses = sum(1 for d in diffs if d < 0)
        severe_losses = sum(1 for d in diffs if d <= -args.severe_loss_points)
        severe_loss_rate = severe_losses / n
        record("incumbent_mean_identity", f"{mean_a:.4f}")
        record("candidate_mean_identity", f"{mean_b:.4f}")
        record("mean_improvement", f"{delta:.4f}")
        record("min_improvement_required", args.min_improvement)
        record("candidate_wins", wins)
        record("candidate_losses", losses)
        record("ties", n - wins - losses)
        record("severe_loss_points", args.severe_loss_points)
        record("severe_losses", severe_losses)
        record("severe_loss_rate", f"{severe_loss_rate:.6f}")
        record("max_severe_loss_rate", args.max_severe_loss_rate)

        p_value = float("nan")
        try:
            from scipy.stats import wilcoxon
            if len([d for d in diffs if d != 0]) >= 10:
                _stat, p_value = wilcoxon(a, b)
        except ImportError:
            pass
        record("wilcoxon_p", "NA" if p_value != p_value else f"{p_value:.4g}")

        big_enough = delta >= args.min_improvement
        significant = p_value == p_value and p_value <= args.max_alpha
        tail_safe = severe_loss_rate <= args.max_severe_loss_rate
        yield_safe = True
        readback_safe = True
        if all(support_args):
            incumbent_lengths = primary_lengths(args.incumbent_contigs)
            candidate_lengths = primary_lengths(args.candidate_contigs)
            incumbent_yield = sum(v >= args.min_primary_len
                                  for v in incumbent_lengths.values())
            candidate_yield = sum(v >= args.min_primary_len
                                  for v in candidate_lengths.values())
            incumbent_readback = readback_summary(args.incumbent_readback)
            candidate_readback = readback_summary(args.candidate_readback)
            yield_ratio = (candidate_yield / incumbent_yield
                           if incumbent_yield else float("inf"))
            readback_ratio = (candidate_readback["supported_bp"] /
                              incumbent_readback["supported_bp"]
                              if incumbent_readback["supported_bp"] else float("inf"))
            yield_safe = yield_ratio >= args.min_yield_ratio
            readback_safe = readback_ratio >= args.min_readback_bp_ratio
            record("min_primary_len", args.min_primary_len)
            record("incumbent_primary_yield", incumbent_yield)
            record("candidate_primary_yield", candidate_yield)
            record("primary_yield_ratio", f"{yield_ratio:.6f}")
            record("min_yield_ratio", args.min_yield_ratio)
            record("incumbent_readback_contigs", incumbent_readback["contigs_checked"])
            record("candidate_readback_contigs", candidate_readback["contigs_checked"])
            record("incumbent_two_end_supported", incumbent_readback["two_end_supported"])
            record("candidate_two_end_supported", candidate_readback["two_end_supported"])
            record("incumbent_readback_supported_bp",
                   f"{incumbent_readback['supported_bp']:.2f}")
            record("candidate_readback_supported_bp",
                   f"{candidate_readback['supported_bp']:.2f}")
            record("readback_supported_bp_ratio", f"{readback_ratio:.6f}")
            record("min_readback_bp_ratio", args.min_readback_bp_ratio)

        promote = big_enough and significant and tail_safe and yield_safe and readback_safe

        if promote:
            record("decision", "promote_candidate")
            record("reason", f"candidate gained {delta:.4f} identity points "
                             f"(>= {args.min_improvement})")
        elif not big_enough:
            record("decision", "keep_incumbent")
            record("reason", f"candidate gained only {delta:.4f} identity points "
                             f"(< {args.min_improvement})")
        else:
            record("decision", "keep_incumbent")
            if not significant:
                record("reason", f"improvement {delta:.4f} not significant "
                                 f"(p={p_value:.4g}, require <= {args.max_alpha})")
            elif not tail_safe:
                record("reason", f"severe-loss rate {severe_loss_rate:.4%} exceeds "
                                 f"{args.max_severe_loss_rate:.4%}")
            elif not yield_safe:
                record("reason", f"primary-contig yield ratio {yield_ratio:.4f} is below "
                                 f"{args.min_yield_ratio:.4f}")
            else:
                record("reason", f"readback-supported-bp ratio {readback_ratio:.4f} is below "
                                 f"{args.min_readback_bp_ratio:.4f}")

    chosen = args.candidate_model if promote else args.incumbent_model
    record("model_in_use", chosen)

    with open(args.out_decision, "w") as fh:
        fh.write("key\tvalue\n")
        fh.write("\n".join(lines) + "\n")

    import shutil
    shutil.copyfile(chosen, args.out_promoted_model)

    for line in lines:
        print(line, file=sys.stderr)
    print(f"\ndecision written: {args.out_decision}", file=sys.stderr)
    print(f"model in use:     {args.out_promoted_model} (copied from {chosen})",
          file=sys.stderr)


if __name__ == "__main__":
    main()
