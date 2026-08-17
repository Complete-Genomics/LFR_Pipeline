#!/usr/bin/env python3
"""Detect Zymo-relative read-pool drift and optionally salvage SE reads.

The input must be grouped by barcode, as denovo/data_R2_sorted.tsv is after
denovo_preprocess.smk.  Two shipped Zymo distributions define the reference:
reads per UMI and read length.  Automatic salvage is deliberately narrow: both
distributions must show large, adverse drift, and the projected retained pool
must still be deep enough.  The salvage operation keeps reads >=300 bp and at
most 300 reads per UMI by default.

Input formats:
  tsv    BX:Z:<barcode> @<read-id>\t<sequence>
  fastq  four-line FASTQ whose header contains BX:Z:<barcode> (the local
         benchmark files also support the legacy #<barcode>/2 spelling)
"""

import argparse
import gzip
import hashlib
import math
import os
import shutil
from collections import defaultdict


DEPTH_BIN_UPPERS = (24, 44, 74, 99, 149, 199, 299, 499, math.inf)
LENGTH_BIN_UPPERS = (99, 199, 299, 399, 499, 599, 699, math.inf)

# PSI >=0.25 conventionally denotes a large population shift.  A second,
# directional effect-size guard prevents a merely deeper (and healthier)
# library from being called degraded.
PSI_LARGE_DRIFT = 0.25
TOP1PCT_ADVERSE_DELTA = 0.10
SHORT_READ_ADVERSE_DELTA = 0.10
SAFE_POST_MEDIAN = 45
SAFE_POST_BELOW25_FRACTION = 0.10

DEPTH_BASELINE_NAME = "zymo_reads_per_umi_distribution.tsv"
LENGTH_BASELINE_NAME = "zymo_read_length_distribution.tsv"


def open_text(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def barcode_from_header(header):
    for field in header.rstrip("\n").split():
        if field.startswith("BX:Z:"):
            barcode = field[5:]
            if barcode:
                return barcode
    if "#" in header:
        barcode = header.rsplit("#", 1)[1].split("/", 1)[0].split()[0]
        if barcode:
            return barcode
    raise ValueError("FASTQ header has no BX:Z: tag or #<barcode>/ suffix: {}"
                     .format(header.rstrip("\n")[:120]))


def iter_fastq(path):
    with open_text(path) as fh:
        record = 0
        while True:
            header = fh.readline()
            if not header:
                return
            sequence = fh.readline()
            plus = fh.readline()
            quality = fh.readline()
            record += 1
            if not sequence or not plus or not quality:
                raise ValueError("truncated FASTQ record {} in {}".format(record, path))
            if not header.startswith("@") or not plus.startswith("+"):
                raise ValueError("malformed FASTQ record {} in {}".format(record, path))
            sequence = sequence.rstrip("\r\n")
            quality = quality.rstrip("\r\n")
            if len(sequence) != len(quality):
                raise ValueError("sequence/quality length mismatch at record {} in {}"
                                 .format(record, path))
            yield barcode_from_header(header), sequence


def iter_tsv(path):
    with open_text(path) as fh:
        for line_no, line in enumerate(fh, 1):
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) < 2:
                raise ValueError("expected at least two TSV columns at {}:{}"
                                 .format(path, line_no))
            tag = fields[0].split(None, 1)[0]
            if not tag.startswith("BX:Z:") or len(tag) == 5:
                raise ValueError("missing BX:Z:<barcode> at {}:{}"
                                 .format(path, line_no))
            yield tag[5:], fields[1]


def percentile_from_hist(histogram, fraction):
    n = sum(histogram.values())
    if not n:
        return 0
    rank = max(1, int(math.ceil(fraction * n)))
    seen = 0
    for value in sorted(histogram):
        seen += histogram[value]
        if seen >= rank:
            return value
    return max(histogram)


def _update_top_two(top, depth, short):
    top.append((depth, short))
    top.sort(reverse=True)
    del top[2:]


def _binned_counts(histogram, uppers):
    counts = [0] * len(uppers)
    for value, count in histogram.items():
        for i, upper in enumerate(uppers):
            if value <= upper:
                counts[i] += count
                break
    return counts


def population_stability_index(observed_counts, expected_counts):
    """Return PSI with a small pseudocount to keep empty bins finite."""
    if len(observed_counts) != len(expected_counts):
        raise ValueError("observed and baseline distributions have different bins")
    observed_total = sum(observed_counts)
    expected_total = sum(expected_counts)
    if not observed_total or not expected_total:
        raise ValueError("observed distribution is empty")
    pseudo = 0.5
    observed_denominator = observed_total + pseudo * len(observed_counts)
    expected_denominator = expected_total + pseudo * len(expected_counts)
    psi = 0.0
    for observed_count, expected_count in zip(observed_counts, expected_counts):
        observed = (observed_count + pseudo) / observed_denominator
        expected = (expected_count + pseudo) / expected_denominator
        psi += (observed - expected) * math.log(observed / expected)
    return psi


def _collect(reads, min_read_length, max_reads_per_umi):
    depth_hist = defaultdict(int)
    length_hist = defaultdict(int)
    post_depth_hist = defaultdict(int)
    by_depth = defaultdict(lambda: [0, 0, 0])  # UMIs, reads, short reads
    top_two = []
    total_reads = 0
    short_reads = 0
    n_umi = 0
    new_below25 = 0
    current = None
    depth = 0
    short = 0

    def finish_group(group_depth, group_short):
        nonlocal n_umi, new_below25
        if group_depth == 0:
            return
        long_depth = group_depth - group_short
        post_depth = min(long_depth, max_reads_per_umi)
        n_umi += 1
        depth_hist[group_depth] += 1
        post_depth_hist[post_depth] += 1
        bucket = by_depth[group_depth]
        bucket[0] += 1
        bucket[1] += group_depth
        bucket[2] += group_short
        _update_top_two(top_two, group_depth, group_short)
        if group_depth >= 25 and post_depth < 25:
            new_below25 += 1

    for barcode, sequence in reads:
        if current is None:
            current = barcode
        elif barcode != current:
            finish_group(depth, short)
            current = barcode
            depth = 0
            short = 0
        sequence_length = len(sequence)
        depth += 1
        total_reads += 1
        length_hist[sequence_length] += 1
        if sequence_length < min_read_length:
            short += 1
            short_reads += 1
    finish_group(depth, short)

    if total_reads == 0:
        raise ValueError("input contains no reads")

    top_n = max(1, int(math.ceil(n_umi * 0.01)))
    remaining = top_n
    top1pct_reads = 0.0
    top1pct_short = 0.0
    for group_depth in sorted(by_depth, reverse=True):
        groups, reads_at_depth, short_at_depth = by_depth[group_depth]
        take = min(remaining, groups)
        share = take / groups
        top1pct_reads += reads_at_depth * share
        top1pct_short += short_at_depth * share
        remaining -= take
        if remaining == 0:
            break

    post_reads = sum(value * count for value, count in post_depth_hist.items())
    post_below25 = sum(count for value, count in post_depth_hist.items()
                       if value < 25)
    post_below45 = sum(count for value, count in post_depth_hist.items()
                       if value < 45)
    top1_reads = top_two[0][0] if top_two else 0
    top2_reads = sum(item[0] for item in top_two)

    metrics = {
        "total_reads": total_reads,
        "n_umi": n_umi,
        "read_length_median": percentile_from_hist(length_hist, 0.50),
        "read_length_p10": percentile_from_hist(length_hist, 0.10),
        "depth_median": percentile_from_hist(depth_hist, 0.50),
        "depth_p90": percentile_from_hist(depth_hist, 0.90),
        "depth_p99": percentile_from_hist(depth_hist, 0.99),
        "depth_max": max(depth_hist),
        "top1_read_fraction": top1_reads / total_reads,
        "top2_read_fraction": top2_reads / total_reads,
        "top1pct_read_fraction": top1pct_reads / total_reads,
        "short_read_fraction": short_reads / total_reads,
        "top1pct_short_fraction": (top1pct_short / top1pct_reads
                                     if top1pct_reads else 0.0),
        "projected_postfilter_reads": post_reads,
        "projected_read_retention": post_reads / total_reads,
        "projected_depth_median": percentile_from_hist(post_depth_hist, 0.50),
        "projected_umi_below25": post_below25,
        "projected_umi_below25_fraction": post_below25 / n_umi,
        "projected_new_umi_below25": new_below25,
        "projected_new_umi_below25_fraction": new_below25 / n_umi,
        "projected_umi_below45": post_below45,
        "projected_umi_below45_fraction": post_below45 / n_umi,
    }
    return metrics, dict(depth_hist), dict(length_hist)


def summarize(reads, min_read_length=300, max_reads_per_umi=300,
              baseline=None):
    """Summarize reads and, when supplied, compare against a baseline."""
    metrics, depth_hist, length_hist = _collect(
        reads, min_read_length, max_reads_per_umi)
    if baseline is None:
        return metrics

    depth_ref = baseline["depth"]
    length_ref = baseline["length"]
    depth_psi = population_stability_index(
        _binned_counts(depth_hist, depth_ref["uppers"]), depth_ref["counts"])
    length_psi = population_stability_index(
        _binned_counts(length_hist, length_ref["uppers"]),
        length_ref["counts"])
    top_delta = (metrics["top1pct_read_fraction"] -
                 depth_ref["top1pct_read_fraction"])
    short_delta = (metrics["short_read_fraction"] -
                   length_ref["short_read_fraction"])
    depth_drift = depth_psi >= PSI_LARGE_DRIFT
    length_drift = length_psi >= PSI_LARGE_DRIFT
    depth_adverse = depth_drift and top_delta >= TOP1PCT_ADVERSE_DELTA
    length_adverse = length_drift and short_delta >= SHORT_READ_ADVERSE_DELTA
    postfilter_safe = (
        metrics["projected_depth_median"] >= SAFE_POST_MEDIAN and
        metrics["projected_umi_below25_fraction"] <=
        SAFE_POST_BELOW25_FRACTION
    )
    if depth_adverse and length_adverse and postfilter_safe:
        verdict = "salvage_candidate"
    elif not depth_adverse and not length_adverse:
        verdict = "pass_through"
    else:
        verdict = "drift_no_salvage"

    metrics.update({
        "depth_psi": depth_psi,
        "read_length_psi": length_psi,
        "baseline_top1pct_read_fraction": depth_ref["top1pct_read_fraction"],
        "top1pct_read_fraction_delta": top_delta,
        "baseline_short_read_fraction": length_ref["short_read_fraction"],
        "short_read_fraction_delta": short_delta,
        "depth_distribution_drift": depth_drift,
        "read_length_distribution_drift": length_drift,
        "depth_adverse_drift": depth_adverse,
        "read_length_adverse_drift": length_adverse,
        "projected_postfilter_safe": postfilter_safe,
        "candidate_verdict": verdict,
    })
    return metrics


def _render_upper(upper):
    return "inf" if math.isinf(upper) else str(int(upper))


def write_distribution(path, name, histogram, uppers, metrics, source_path):
    counts = _binned_counts(histogram, uppers)
    total = sum(counts)
    digest = hashlib.sha256()
    with open(source_path, "rb") as source:
        for chunk in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(chunk)
    with open(path, "w") as out:
        out.write("# model={}\n".format(name))
        out.write("# source={}\n".format(os.path.basename(source_path)))
        out.write("# source_sha256={}\n".format(digest.hexdigest()))
        out.write("# n={}\n".format(total))
        out.write("# top1pct_read_fraction={:.12f}\n".format(
            metrics["top1pct_read_fraction"]))
        out.write("# short_read_fraction={:.12f}\n".format(
            metrics["short_read_fraction"]))
        out.write("bin_lower\tbin_upper\tcount\tfraction\n")
        lower = 0
        for upper, count in zip(uppers, counts):
            out.write("{}\t{}\t{}\t{:.12f}\n".format(
                lower, _render_upper(upper), count, count / total))
            lower = int(upper) + 1 if not math.isinf(upper) else lower


def build_baseline(reads, output_dir, input_path, min_read_length,
                   max_reads_per_umi):
    # Baseline construction is a one-time operation and intentionally accepts
    # the original unsorted FASTQ.  Runtime input remains barcode-grouped so
    # sample QC can stream in O(number of depth/length values) memory.
    barcode_depth = defaultdict(int)
    length_hist = defaultdict(int)
    total_reads = 0
    short_reads = 0
    for barcode, sequence in reads:
        barcode_depth[barcode] += 1
        length_hist[len(sequence)] += 1
        total_reads += 1
        if len(sequence) < min_read_length:
            short_reads += 1
    if not total_reads:
        raise ValueError("input contains no reads")
    depth_hist = defaultdict(int)
    for depth in barcode_depth.values():
        depth_hist[depth] += 1
    top_n = max(1, int(math.ceil(len(barcode_depth) * 0.01)))
    top1pct_reads = sum(sorted(barcode_depth.values(), reverse=True)[:top_n])
    metrics = {
        "total_reads": total_reads,
        "n_umi": len(barcode_depth),
        "top1pct_read_fraction": top1pct_reads / total_reads,
        "short_read_fraction": short_reads / total_reads,
    }
    os.makedirs(output_dir, exist_ok=True)
    write_distribution(os.path.join(output_dir, DEPTH_BASELINE_NAME),
                       "zymo_reads_per_umi", depth_hist, DEPTH_BIN_UPPERS,
                       metrics, input_path)
    write_distribution(os.path.join(output_dir, LENGTH_BASELINE_NAME),
                       "zymo_read_length", length_hist, LENGTH_BIN_UPPERS,
                       metrics, input_path)
    return metrics


def load_distribution(path):
    metadata = {}
    uppers = []
    counts = []
    fractions = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("# "):
                key, value = line[2:].rstrip("\n").split("=", 1)
                metadata[key] = value
                continue
            if line.startswith("bin_lower"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 4:
                raise ValueError("malformed baseline row in {}".format(path))
            uppers.append(math.inf if fields[1] == "inf" else int(fields[1]))
            counts.append(int(fields[2]))
            fractions.append(float(fields[3]))
    if not uppers or abs(sum(fractions) - 1.0) > 1e-6:
        raise ValueError("invalid distribution in {}".format(path))
    return {
        "uppers": tuple(uppers),
        "counts": counts,
        "fractions": fractions,
        "top1pct_read_fraction": float(metadata["top1pct_read_fraction"]),
        "short_read_fraction": float(metadata["short_read_fraction"]),
        "source_sha256": metadata.get("source_sha256", ""),
    }


def load_baseline(depth_path, length_path):
    return {
        "depth": load_distribution(depth_path),
        "length": load_distribution(length_path),
    }


def _copy_passthrough(source, destination):
    temp = destination + ".tmp"
    try:
        os.link(source, temp)
    except OSError:
        shutil.copyfile(source, temp)
    os.replace(temp, destination)


def filter_tsv(source, destination, min_read_length, max_reads_per_umi):
    """Apply deterministic length/cap salvage to a barcode-grouped TSV."""
    temp = destination + ".tmp"
    current = None
    kept_in_group = 0
    kept = 0
    with open_text(source) as inp, open(temp, "w") as out:
        for line_no, line in enumerate(inp, 1):
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) < 2:
                raise ValueError("expected at least two TSV columns at {}:{}"
                                 .format(source, line_no))
            tag = fields[0].split(None, 1)[0]
            if not tag.startswith("BX:Z:"):
                raise ValueError("missing BX tag at {}:{}".format(source, line_no))
            barcode = tag[5:]
            if barcode != current:
                current = barcode
                kept_in_group = 0
            if (len(fields[1]) >= min_read_length and
                    kept_in_group < max_reads_per_umi):
                out.write(line)
                kept_in_group += 1
                kept += 1
    os.replace(temp, destination)
    return kept


def write_report(path, metrics, input_path, input_format, min_read_length,
                 max_reads_per_umi, mode="report", action="report_only",
                 depth_baseline=None, length_baseline=None):
    with open(path, "w") as out:
        out.write("setting\tvalue\n")
        out.write("mode\t{}\n".format(mode))
        out.write("action\t{}\n".format(action))
        out.write("input\t{}\n".format(os.path.abspath(input_path)))
        out.write("input_format\t{}\n".format(input_format))
        if depth_baseline:
            out.write("depth_baseline\t{}\n".format(
                os.path.abspath(depth_baseline)))
        if length_baseline:
            out.write("read_length_baseline\t{}\n".format(
                os.path.abspath(length_baseline)))
        out.write("min_read_length\t{}\n".format(min_read_length))
        out.write("max_reads_per_umi\t{}\n".format(max_reads_per_umi))
        out.write("threshold_psi_large_drift\t{:.6f}\n".format(
            PSI_LARGE_DRIFT))
        out.write("threshold_top1pct_adverse_delta\t{:.6f}\n".format(
            TOP1PCT_ADVERSE_DELTA))
        out.write("threshold_short_read_adverse_delta\t{:.6f}\n".format(
            SHORT_READ_ADVERSE_DELTA))
        out.write("threshold_safe_post_median\t{}\n".format(SAFE_POST_MEDIAN))
        out.write("threshold_safe_post_below25_fraction\t{:.6f}\n".format(
            SAFE_POST_BELOW25_FRACTION))
        for key, value in metrics.items():
            if isinstance(value, bool):
                rendered = str(value)
            elif isinstance(value, float):
                rendered = "{:.6f}".format(value)
            else:
                rendered = str(value)
            out.write("{}\t{}\n".format(key, rendered))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--r2", required=True,
                        help="barcode-grouped R2 TSV or FASTQ")
    parser.add_argument("--r2-format", choices=("tsv", "fastq"), default="tsv")
    parser.add_argument("--min-read-length", type=int, default=300)
    parser.add_argument("--max-reads-per-umi", type=int, default=300)
    parser.add_argument("--depth-baseline")
    parser.add_argument("--read-length-baseline")
    parser.add_argument("--mode", choices=("off", "report", "auto"),
                        default="report")
    parser.add_argument("--reads-out")
    parser.add_argument("--out")
    parser.add_argument("--build-baseline-dir")
    parser.add_argument("--apply-decision",
                        help="apply an existing decision report without rescanning")
    args = parser.parse_args()

    if args.min_read_length < 1 or args.max_reads_per_umi < 1:
        parser.error("--min-read-length and --max-reads-per-umi must be positive")
    iterator = iter_tsv(args.r2) if args.r2_format == "tsv" else iter_fastq(args.r2)

    if args.apply_decision:
        if not args.reads_out or args.out or args.build_baseline_dir:
            parser.error("--apply-decision requires --reads-out and no --out")
        with open(args.apply_decision) as report:
            decision = dict(line.rstrip("\n").split("\t", 1)
                            for line in report if "\t" in line)
        action = decision.get("action")
        if action == "salvage":
            if args.r2_format != "tsv":
                parser.error("automatic salvage currently requires TSV input")
            kept = filter_tsv(
                args.r2, args.reads_out,
                int(decision["min_read_length"]),
                int(decision["max_reads_per_umi"]))
            print("applied noisy preprocess decision: salvage ({} reads)".format(kept))
        elif action in ("pass_through", "report_only", "disabled"):
            _copy_passthrough(args.r2, args.reads_out)
            print("applied noisy preprocess decision: pass_through")
        else:
            parser.error("unrecognized action in {}: {!r}".format(
                args.apply_decision, action))
        return

    if args.build_baseline_dir:
        if args.out or args.reads_out or args.depth_baseline or args.read_length_baseline:
            parser.error("--build-baseline-dir cannot be combined with report/output options")
        metrics = build_baseline(iterator, args.build_baseline_dir, args.r2,
                                 args.min_read_length, args.max_reads_per_umi)
        print("wrote Zymo baselines: {} reads, {} UMIs".format(
            metrics["total_reads"], metrics["n_umi"]))
        return

    if not args.out:
        parser.error("--out is required outside --build-baseline-dir mode")
    if args.mode != "off" and (not args.depth_baseline or
                               not args.read_length_baseline):
        parser.error("--depth-baseline and --read-length-baseline are required")
    if args.mode == "off":
        metrics = {"candidate_verdict": "not_evaluated"}
        action = "disabled"
    else:
        baseline = load_baseline(args.depth_baseline, args.read_length_baseline)
        metrics = summarize(iterator, args.min_read_length,
                            args.max_reads_per_umi, baseline)
        action = "report_only"
        if args.mode == "auto":
            if metrics["candidate_verdict"] == "salvage_candidate":
                action = "salvage"
                if args.reads_out:
                    if args.r2_format != "tsv":
                        parser.error("automatic salvage currently requires TSV input")
                    metrics["actual_postfilter_reads"] = filter_tsv(
                        args.r2, args.reads_out, args.min_read_length,
                        args.max_reads_per_umi)
            else:
                action = "pass_through"
                if args.reads_out:
                    _copy_passthrough(args.r2, args.reads_out)
    if args.mode in ("off", "report") and args.reads_out:
        _copy_passthrough(args.r2, args.reads_out)

    write_report(args.out, metrics, args.r2, args.r2_format,
                 args.min_read_length, args.max_reads_per_umi,
                 mode=args.mode, action=action,
                 depth_baseline=args.depth_baseline,
                 length_baseline=args.read_length_baseline)
    print("noisy preprocess QC: {} (action={})".format(
        metrics["candidate_verdict"], action))


if __name__ == "__main__":
    main()
