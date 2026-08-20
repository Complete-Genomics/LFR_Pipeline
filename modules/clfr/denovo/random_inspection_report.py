#!/usr/bin/env python3
"""Write a denominator-explicit QA report for random_inspection.smk.

The local-reference identity comparison is intentionally separate from the
optional strict chimera labels.  The former works on an arbitrary community;
the latter needs a small, known-composition control reference whose headers
encode a biological source (for Zymo: ``Bacillus_subtilis_16S_1``).
"""
import argparse
import csv
import math
import os
from collections import Counter, defaultdict


ARMS = ("nofilter", "plain", "ml")
IDENTITY_BINS = [(0, 80), (80, 85), (85, 90), (90, 95), (95, 97), (97, 99), (99, 100.0001)]
LENGTH_BINS = [(0, 500), (500, 750), (750, 1000), (1000, 1250), (1250, 1500), (1500, 2000), (2000, math.inf)]
READ_BINS = [(0, 1), (1, 3), (3, 6), (6, 11), (11, 21), (21, 51), (51, 101), (101, math.inf)]


def barcode_from_contig(value):
    return value.split(">", 1)[0]


def label_for_target(target, separator):
    if separator and separator in target:
        return target.split(separator, 1)[0]
    return target


def read_primary_fasta(path):
    sequences = {}
    name = None
    chunks = []

    def save():
        if name is not None:
            barcode = barcode_from_contig(name)
            sequence = "".join(chunks)
            if len(sequence) > len(sequences.get(barcode, "")):
                sequences[barcode] = sequence

    with open(path) as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith(">"):
                save()
                name = line[1:].split()[0]
                chunks = []
            elif line:
                chunks.append(line)
    save()
    return sequences


def barcode_from_read_id(read_id):
    if not read_id.startswith("BX:Z:") or len(read_id) < 20:
        raise ValueError("unexpected sgrep read ID: {!r}".format(read_id))
    return read_id[5:20]


def read_counts(path):
    counts = Counter()
    with open(path) as handle:
        for line_number, line in enumerate(handle, 1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2:
                continue
            try:
                counts[barcode_from_read_id(fields[0])] += 1
            except ValueError as exc:
                raise ValueError("{}:{}: {}".format(path, line_number, exc))
    return counts


def load_paired_identity(path):
    values = defaultdict(list)
    with open(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"plain", "ml"}
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise ValueError("{} must contain plain and ml identity columns".format(path))
        for row in reader:
            for arm in ("plain", "ml"):
                try:
                    values[arm].append(float(row[arm]))
                except (TypeError, ValueError):
                    continue
    return values


def strict_chimera_labels(path, separator):
    """Return barcode -> 0/1 for strict, reference-defined labels only.

    A quarter is decisive at best identity >=97% and a >=3-point advantage
    over its second biological source.  A barcode is labelled only when at
    least two quarters are decisive; it is a chimera when those sources differ.
    """
    hits = defaultdict(lambda: defaultdict(dict))
    with open(path) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3 or "||Q" not in fields[0]:
                continue
            contig, quarter = fields[0].rsplit("||Q", 1)
            try:
                identity = float(fields[2])
            except ValueError:
                continue
            barcode = barcode_from_contig(contig)
            label = label_for_target(fields[1], separator)
            hits[barcode][quarter][label] = max(identity, hits[barcode][quarter].get(label, -1.0))

    labels = {}
    for barcode, quarters in hits.items():
        decisive = []
        for source_hits in quarters.values():
            ranked = sorted(source_hits.items(), key=lambda item: (-item[1], item[0]))
            if not ranked:
                continue
            source, best = ranked[0]
            second = ranked[1][1] if len(ranked) > 1 else 0.0
            if best >= 97.0 and best - second >= 3.0:
                decisive.append(source)
        if len(decisive) >= 2:
            labels[barcode] = int(len(set(decisive)) >= 2)
    return labels


def percentile(values, probability):
    if not values:
        return float("nan")
    values = sorted(values)
    index = (len(values) - 1) * probability
    lo, hi = int(math.floor(index)), int(math.ceil(index))
    if lo == hi:
        return values[lo]
    return values[lo] + (values[hi] - values[lo]) * (index - lo)


def distribution(values):
    values = list(values)
    if not values:
        return {key: float("nan") for key in ("n", "mean", "min", "p10", "p25", "median", "p75", "p90", "max")}
    return {
        "n": len(values), "mean": sum(values) / len(values), "min": min(values),
        "p10": percentile(values, 0.10), "p25": percentile(values, 0.25),
        "median": percentile(values, 0.50), "p75": percentile(values, 0.75),
        "p90": percentile(values, 0.90), "max": max(values),
    }


def histogram(rows, arm, metric, values, bins):
    counts = [0] * len(bins)
    for value in values:
        for index, (lower, upper) in enumerate(bins):
            if lower <= value < upper:
                counts[index] += 1
                break
    total = len(values)
    for (lower, upper), count in zip(bins, counts):
        upper_text = "inf" if math.isinf(upper) else "{:g}".format(upper)
        rows.append({"arm": arm, "metric": metric, "bin": "[{:g}, {})".format(lower, upper_text),
                     "count": count, "fraction": count / total if total else "", "denominator": total})


def format_value(value):
    if isinstance(value, float):
        if math.isnan(value):
            return "NA"
        return "{:.2f}".format(value)
    return str(value)


def markdown_table(rows, fields):
    header = "| " + " | ".join(fields) + " |"
    separator = "|" + "|".join("---:" if field not in ("arm", "metric") else "---" for field in fields) + "|"
    body = ["| " + " | ".join(format_value(row.get(field, "")) for field in fields) + " |" for row in rows]
    return "\n".join([header, separator] + body)


def write_tsv(path, rows, fields):
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for arm in ARMS:
        parser.add_argument("--{}-reads".format(arm), required=True)
        parser.add_argument("--{}-contigs".format(arm), required=True)
    parser.add_argument("--paired-identity", required=True)
    parser.add_argument("--quarter-hits-nofilter")
    parser.add_argument("--quarter-hits-plain")
    parser.add_argument("--quarter-hits-ml")
    parser.add_argument("--chimera-label-separator", default="_16S")
    parser.add_argument("--out-md", required=True)
    parser.add_argument("--out-metrics", required=True)
    parser.add_argument("--out-identity-hist", required=True)
    parser.add_argument("--out-contig-len-hist", required=True)
    parser.add_argument("--out-reads-hist", required=True)
    args = parser.parse_args()

    quarter_paths = [args.quarter_hits_nofilter, args.quarter_hits_plain, args.quarter_hits_ml]
    if any(quarter_paths) and not all(quarter_paths):
        parser.error("provide quarter-hit files for all three arms, or for none")

    reads = {arm: read_counts(getattr(args, "{}_reads".format(arm))) for arm in ARMS}
    contigs = {arm: read_primary_fasta(getattr(args, "{}_contigs".format(arm))) for arm in ARMS}
    input_umis = set(reads["nofilter"])
    if not input_umis:
        raise ValueError("nofilter read pool contains no UMI")
    identities = load_paired_identity(args.paired_identity)

    metric_rows = []
    length_rows, read_rows, identity_rows = [], [], []
    for arm in ARMS:
        read_values = [reads[arm].get(barcode, 0) for barcode in input_umis]
        contig_values = [len(sequence) for sequence in contigs[arm].values()]
        read_summary = distribution(read_values)
        length_summary = distribution(contig_values)
        assembled = len(set(contigs[arm]) & input_umis)
        metric_rows.extend([
            {"section": "assembly", "arm": arm, "metric": "input_umis", "value": len(input_umis), "denominator": len(input_umis)},
            {"section": "assembly", "arm": arm, "metric": "read_bearing_umis", "value": sum(value > 0 for value in read_values), "denominator": len(input_umis)},
            {"section": "assembly", "arm": arm, "metric": "assembled_umis", "value": assembled, "denominator": len(input_umis)},
            {"section": "assembly", "arm": arm, "metric": "assembly_rate", "value": assembled / len(input_umis), "denominator": len(input_umis)},
        ])
        for metric, value in read_summary.items():
            metric_rows.append({"section": "reads_per_umi", "arm": arm, "metric": metric, "value": value, "denominator": len(input_umis)})
        for metric, value in length_summary.items():
            metric_rows.append({"section": "primary_contig_length", "arm": arm, "metric": metric, "value": value, "denominator": len(contig_values)})
        histogram(read_rows, arm, "reads_per_input_umi", read_values, READ_BINS)
        histogram(length_rows, arm, "primary_contig_bp", contig_values, LENGTH_BINS)

    for arm in ("plain", "ml"):
        values = identities[arm]
        summary = distribution(values)
        for metric, value in summary.items():
            metric_rows.append({"section": "identity_to_fixed_local_reference", "arm": arm,
                                "metric": metric, "value": value, "denominator": len(values)})
        for threshold, name, compare in ((97.0, "at_least_97_pct", lambda x, t: x >= t),
                                         (90.0, "below_90_pct", lambda x, t: x < t)):
            metric_rows.append({"section": "identity_to_fixed_local_reference", "arm": arm,
                                "metric": name, "value": sum(compare(x, threshold) for x in values) / len(values) if values else float("nan"),
                                "denominator": len(values)})
        histogram(identity_rows, arm, "identity_pct", values, IDENTITY_BINS)

    chimera_rows = []
    if all(quarter_paths):
        for arm, path in zip(ARMS, quarter_paths):
            labels = strict_chimera_labels(path, args.chimera_label_separator)
            labels = {barcode: label for barcode, label in labels.items() if barcode in contigs[arm]}
            total = len(labels)
            chimeras = sum(labels.values())
            chimera_rows.append({"arm": arm, "assembled_umis": len(set(contigs[arm]) & input_umis),
                                 "strictly_classifiable": total, "strict_clean": total - chimeras,
                                 "strict_chimera": chimeras,
                                 "strict_chimera_rate": chimeras / total if total else float("nan")})
            metric_rows.extend([
                {"section": "strict_reference_chimera", "arm": arm, "metric": "strictly_classifiable", "value": total, "denominator": len(set(contigs[arm]) & input_umis)},
                {"section": "strict_reference_chimera", "arm": arm, "metric": "strict_chimera", "value": chimeras, "denominator": total},
                {"section": "strict_reference_chimera", "arm": arm, "metric": "strict_chimera_rate", "value": chimeras / total if total else float("nan"), "denominator": total},
            ])

    write_tsv(args.out_metrics, metric_rows, ["section", "arm", "metric", "value", "denominator"])
    write_tsv(args.out_identity_hist, identity_rows, ["arm", "metric", "bin", "count", "fraction", "denominator"])
    write_tsv(args.out_contig_len_hist, length_rows, ["arm", "metric", "bin", "count", "fraction", "denominator"])
    write_tsv(args.out_reads_hist, read_rows, ["arm", "metric", "bin", "count", "fraction", "denominator"])

    assembly = [row for row in metric_rows if row["section"] == "assembly" and row["metric"] in ("assembled_umis", "assembly_rate")]
    read_summary = [row for row in metric_rows if row["section"] == "reads_per_umi" and row["metric"] in ("n", "median", "p10", "p90", "max")]
    length_summary = [row for row in metric_rows if row["section"] == "primary_contig_length" and row["metric"] in ("n", "median", "p10", "p90", "max")]
    identity_summary = [row for row in metric_rows if row["section"] == "identity_to_fixed_local_reference" and row["metric"] in ("n", "mean", "median", "p10", "p90", "at_least_97_pct", "below_90_pct")]
    report = [
        "# Random-inspection assembly QA report", "",
        "All assembly rates use the no-filter input UMI set as denominator. Length distributions contain assembled primary (`k41_0`) contigs only. Read distributions include zero for an input UMI that a filter removed completely.",
        "", "## Assembly rate", "", markdown_table(assembly, ["arm", "metric", "value", "denominator"]),
        "", "## Reads per input UMI", "", markdown_table(read_summary, ["arm", "metric", "value", "denominator"]),
        "", "## Primary contig length (bp)", "", markdown_table(length_summary, ["arm", "metric", "value", "denominator"]),
        "", "## Identity to fixed local reference", "",
        "`plain` and `ml` are compared to the same per-UMI reference assigned by the no-filter arm. This is a paired consistency score, not absolute taxonomic truth; the no-filter top-hit identity is intentionally not included because it selected that reference.",
        "", markdown_table(identity_summary, ["arm", "metric", "value", "denominator"]),
    ]
    if chimera_rows:
        report.extend([
            "", "## Strict reference-defined chimera rate", "",
            "A quarter is decisive only at identity >=97% and a >=3-point margin over the next biological source. The chimera-rate denominator is only contigs with at least two decisive quarters; unclassifiable contigs are not counted as clean.",
            "", markdown_table(chimera_rows, ["arm", "assembled_umis", "strictly_classifiable", "strict_clean", "strict_chimera", "strict_chimera_rate"]),
        ])
    else:
        report.extend([
            "", "## Chimera rate", "",
            "Not computed: strict chimera truth requires `random_inspection.chimera_reference`, a small known-composition control reference with biological-source labels in its FASTA headers. Do not substitute the broad Greengenes database for this rate.",
        ])
    report.extend([
        "", "## Distribution files", "",
        "- `{}`: all machine-readable summary metrics".format(os.path.basename(args.out_metrics)),
        "- `{}`: fixed bins for paired identity scores".format(os.path.basename(args.out_identity_hist)),
        "- `{}`: fixed bins for primary-contig lengths".format(os.path.basename(args.out_contig_len_hist)),
        "- `{}`: fixed bins for reads per input UMI".format(os.path.basename(args.out_reads_hist)), "",
    ])
    with open(args.out_md, "w") as handle:
        handle.write("\n".join(report))


if __name__ == "__main__":
    main()
