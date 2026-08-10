#!/usr/bin/env python3
"""Summarize ERCC consensus FASTA length recovery and variants from PAF.

Use an existing minimap2 PAF when available. For exact SNP/indel positions,
generate PAF with minimap2 --cs=long. Without cs tags, the script still reports
alignment coverage and NM mismatch/edit distance when present, but cannot list
per-base variants.
"""

from __future__ import print_function

import argparse
import csv
import math
import os
import re
import sys


ERCC_RE = re.compile(r"(ERCC-\d+)")
CS_RE = re.compile(r"(:[0-9]+|\*[a-z][a-z]|\+[a-z]+|-[a-z]+|=[a-z]+|~[a-z]{2}[0-9]+[a-z]{2})")


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--fasta", required=True, help="Consensus FASTA file.")
    parser.add_argument(
        "--paf",
        default=None,
        help="Optional minimap2 PAF for length and variant evaluation. Not needed when FASTA headers contain ERCC IDs.",
    )
    parser.add_argument("--ercc_ref", default=None, help="ERCC truth table with 'ERCC ID' and 'Sequence'.")
    parser.add_argument("--summary", default="consensus/consensus_ercc_summary.tsv")
    parser.add_argument("--variants", default="consensus/consensus_ercc_variants.tsv")
    parser.add_argument("--length_stats", default="consensus/consensus_ercc_length_stats.tsv")
    parser.add_argument(
        "--fasta_length_stats",
        default=None,
        help="Optional TSV with consensus FASTA sequence-length count, mean, median, min, and max.",
    )
    parser.add_argument(
        "--barcode_vs_concentration",
        default="consensus/consensus_ercc_barcode_vs_concentration",
        help="Output prefix for assembled-barcode count per ERCC ID vs. known Mix 1 concentration.",
    )
    parser.add_argument(
        "--barcode_vs_concentration_plot",
        default=None,
        help="PNG output path. Defaults to <barcode_vs_concentration>.png.",
    )
    parser.add_argument(
        "--include_unobserved_ercc",
        action="store_true",
        help="Include all ERCC truth entries with zero barcode count; by default plot only ERCC IDs present in the consensus FASTA.",
    )
    return parser.parse_args()


def default_ercc_ref():
    here = os.path.dirname(os.path.abspath(__file__))
    return os.path.abspath(os.path.join(here, "..", "calc_frag_len", "ercc_truth.txt"))


def read_fasta_lengths(path):
    lengths = {}
    name = None
    length = 0
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    lengths[name] = length
                name = line[1:].split()[0]
                length = 0
            else:
                length += len(line)
    if name is not None:
        lengths[name] = length
    return lengths


def fasta_length_stats(fasta_lengths):
    values = sorted(fasta_lengths.values())
    count = len(values)
    if not values:
        return {
            "metric": "consensus_fasta_length_bp",
            "count": 0,
            "mean": "",
            "median": "",
            "min": "",
            "max": "",
        }
    midpoint = count // 2
    median = values[midpoint] if count % 2 else (values[midpoint - 1] + values[midpoint]) / 2.0
    return {
        "metric": "consensus_fasta_length_bp",
        "count": count,
        "mean": sum(values) / float(count),
        "median": median,
        "min": values[0],
        "max": values[-1],
    }


def read_ercc_truth(path):
    refs = {}
    with open(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if "ERCC ID" not in reader.fieldnames or "Sequence" not in reader.fieldnames:
            raise ValueError("ERCC reference must contain 'ERCC ID' and 'Sequence'")
        for row in reader:
            ercc_id = row["ERCC ID"].strip()
            seq = row["Sequence"].strip().upper()
            if ercc_id and seq:
                refs[ercc_id] = seq
    return refs


def ercc_from_text(text):
    match = ERCC_RE.search(text)
    return match.group(1) if match else None


def parse_tags(fields):
    tags = {}
    for field in fields[12:]:
        parts = field.split(":", 2)
        if len(parts) == 3:
            tags[parts[0]] = parts[2]
    return tags


def read_best_ercc_paf(path, ercc_refs):
    best = {}
    with open(path) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                continue
            query = fields[0]
            target = fields[5]
            ercc_id = ercc_from_text(target) or ercc_from_text(query)
            if ercc_id not in ercc_refs:
                continue
            record = {
                "query": query,
                "query_len": int(fields[1]),
                "query_start": int(fields[2]),
                "query_end": int(fields[3]),
                "strand": fields[4],
                "target": target,
                "target_len": int(fields[6]),
                "target_start": int(fields[7]),
                "target_end": int(fields[8]),
                "matches": int(fields[9]),
                "block_len": int(fields[10]),
                "mapq": int(fields[11]),
                "tags": parse_tags(fields),
                "ercc_id": ercc_id,
            }
            old = best.get(query)
            if old is None or (record["matches"], record["block_len"], record["mapq"]) > (
                old["matches"], old["block_len"], old["mapq"]
            ):
                best[query] = record
    return best


def count_variants_from_cs(query, ercc_id, cs, target_start, strand):
    variants = []
    snp_count = 0
    ins_count = 0
    del_count = 0
    ref_pos = target_start + 1

    for token in CS_RE.findall(cs):
        op = token[0]
        payload = token[1:]
        if op == ":":
            ref_pos += int(payload)
        elif op == "=":
            ref_pos += len(payload)
        elif op == "*":
            ref_base = payload[0].upper()
            query_base = payload[1].upper()
            snp_count += 1
            variants.append({
                "consensus_id": query,
                "ercc_id": ercc_id,
                "variant_type": "snp",
                "ref_pos_1based": ref_pos,
                "ref_base": ref_base,
                "consensus_base": query_base,
                "strand": strand,
            })
            ref_pos += 1
        elif op == "+":
            ins_count += len(payload)
            variants.append({
                "consensus_id": query,
                "ercc_id": ercc_id,
                "variant_type": "insertion",
                "ref_pos_1based": ref_pos,
                "ref_base": "-",
                "consensus_base": payload.upper(),
                "strand": strand,
            })
        elif op == "-":
            del_count += len(payload)
            variants.append({
                "consensus_id": query,
                "ercc_id": ercc_id,
                "variant_type": "deletion",
                "ref_pos_1based": ref_pos,
                "ref_base": payload.upper(),
                "consensus_base": "-",
                "strand": strand,
            })
            ref_pos += len(payload)
        elif op == "~":
            match = re.match(r"[a-z]{2}([0-9]+)[a-z]{2}", payload)
            if match:
                ref_pos += int(match.group(1))

    return snp_count, ins_count, del_count, variants


def read_ercc_concentrations(path):
    """ERCC ID -> known Mix 1 concentration (attomoles/ul) from the truth table."""
    mix1_col = "concentration in Mix 1 (attomoles/ul)"
    concentrations = {}
    with open(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if "ERCC ID" not in reader.fieldnames or mix1_col not in reader.fieldnames:
            raise ValueError("ERCC reference must contain 'ERCC ID' and '%s'" % mix1_col)
        for row in reader:
            ercc_id = row["ERCC ID"].strip()
            if not ercc_id:
                continue
            try:
                concentrations[ercc_id] = float(row[mix1_col])
            except (TypeError, ValueError):
                continue
    return concentrations


def barcode_from_consensus_id(consensus_id, ercc_id):
    marker = "_" + ercc_id
    if marker in consensus_id:
        barcode = consensus_id.split(marker, 1)[0]
        if barcode:
            return barcode
    return consensus_id


def assembled_barcodes_from_fasta(fasta_lengths):
    """ERCC ID -> unique barcodes encoded in consensus FASTA headers."""
    barcodes = {}
    for consensus_id in fasta_lengths:
        ercc_id = ercc_from_text(consensus_id)
        if not ercc_id:
            continue
        barcodes.setdefault(ercc_id, set()).add(
            barcode_from_consensus_id(consensus_id, ercc_id)
        )
    return barcodes


def assembled_barcodes_from_paf(best_paf):
    """ERCC ID -> unique consensus barcodes when ERCC IDs are absent from FASTA headers."""
    barcodes = {}
    for record in best_paf.values():
        ercc_id = record["ercc_id"]
        barcodes.setdefault(ercc_id, set()).add(
            barcode_from_consensus_id(record["query"], ercc_id)
        )
    return barcodes


def log2p1(value):
    return math.log(value + 1, 2)


def pearson_corr(xs, ys):
    n = len(xs)
    if n < 2:
        return float("nan")
    mean_x = sum(xs) / n
    mean_y = sum(ys) / n
    var_x = sum((x - mean_x) ** 2 for x in xs)
    var_y = sum((y - mean_y) ** 2 for y in ys)
    if var_x == 0 or var_y == 0:
        return float("nan")
    cov = sum((x - mean_x) * (y - mean_y) for x, y in zip(xs, ys))
    return cov / math.sqrt(var_x * var_y)


def spearman_corr(xs, ys):
    def rank(values):
        order = sorted(range(len(values)), key=lambda i: values[i])
        ranks = [0.0] * len(values)
        i = 0
        while i < len(order):
            j = i
            while j + 1 < len(order) and values[order[j + 1]] == values[order[i]]:
                j += 1
            avg_rank = (i + j) / 2.0 + 1
            for k in range(i, j + 1):
                ranks[order[k]] = avg_rank
            i = j + 1
        return ranks

    return pearson_corr(rank(xs), rank(ys))


def linear_fit(xs, ys):
    n = len(xs)
    if n < 2:
        return float("nan"), float("nan")
    mean_x = sum(xs) / n
    mean_y = sum(ys) / n
    var_x = sum((x - mean_x) ** 2 for x in xs)
    if var_x == 0:
        return float("nan"), float("nan")
    cov = sum((x - mean_x) * (y - mean_y) for x, y in zip(xs, ys))
    slope = cov / var_x
    intercept = mean_y - slope * mean_x
    return slope, intercept


def barcode_vs_concentration(barcodes_by_ercc, ercc_ref_path, output_prefix, include_unobserved):
    """Compare assembled-transcript (barcode) yield per ERCC ID against known Mix 1 concentration.

    Each consensus transcript in the PAF comes from one barcode's assembly, so the count of
    distinct consensus transcripts that best-map to an ERCC ID is a proxy for how many barcodes
    successfully assembled that transcript. ERCC concentration is known ground truth, so a clean
    log2-log2 dose-response relationship (like the read/BX-tag check in ercc_count.py) indicates
    assembly yield tracks input abundance rather than dropping out unevenly across the spike-in
    concentration range.
    """
    concentrations = read_ercc_concentrations(ercc_ref_path)
    ercc_ids = set(concentrations) if include_unobserved else set(barcodes_by_ercc)

    rows = []
    for ercc_id in sorted(ercc_ids):
        if ercc_id not in concentrations:
            continue
        rows.append({
            "ercc_id": ercc_id,
            "concentration_mix1_attomoles_ul": concentrations[ercc_id],
            "assembled_barcode_count": len(barcodes_by_ercc.get(ercc_id, set())),
        })

    xs = [log2p1(row["concentration_mix1_attomoles_ul"]) for row in rows]
    ys = [log2p1(row["assembled_barcode_count"]) for row in rows]
    pearson = pearson_corr(xs, ys)
    spearman = spearman_corr(xs, ys)
    slope, intercept = linear_fit(xs, ys)
    r_squared = pearson ** 2 if not math.isnan(pearson) else float("nan")

    stats_row = {
        "n": len(rows),
        "pearson_log2_mix1": pearson,
        "spearman_log2_mix1": spearman,
        "r_squared_log2_mix1": r_squared,
        "slope_log2_mix1": slope,
        "intercept_log2_mix1": intercept,
    }

    write_table(
        "%s.tsv" % output_prefix, rows,
        ["ercc_id", "concentration_mix1_attomoles_ul", "assembled_barcode_count"],
    )
    write_table(
        "%s_correlation.tsv" % output_prefix, [stats_row],
        ["n", "pearson_log2_mix1", "spearman_log2_mix1", "r_squared_log2_mix1", "slope_log2_mix1", "intercept_log2_mix1"],
    )
    return rows, stats_row


def plot_barcode_vs_concentration(rows, stats_row, output_path, fasta_stats):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    directory = os.path.dirname(output_path)
    if directory:
        os.makedirs(directory, exist_ok=True)

    xs = [log2p1(row["concentration_mix1_attomoles_ul"]) for row in rows]
    ys = [log2p1(row["assembled_barcode_count"]) for row in rows]
    figure, axis = plt.subplots(figsize=(8, 6))
    axis.scatter(xs, ys, color="#0072B2", edgecolors="white", linewidths=0.5, alpha=0.85)

    slope = stats_row["slope_log2_mix1"]
    intercept = stats_row["intercept_log2_mix1"]
    if xs and not math.isnan(slope) and not math.isnan(intercept):
        x_min, x_max = min(xs), max(xs)
        axis.plot(
            [x_min, x_max],
            [slope * x_min + intercept, slope * x_max + intercept],
            color="#D55E00",
            linewidth=1.5,
        )

    axis.set_title("ERCC Consensus Barcode Yield")
    axis.set_xlabel("Mix 1 concentration: log2(attomoles/uL + 1)")
    axis.set_ylabel("Consensus barcodes: log2(count + 1)")
    axis.grid(axis="both", alpha=0.25)
    axis.text(
        0.03,
        0.97,
        "ERCC n={n}\nPearson r={pearson_log2_mix1:.3f}\nSpearman rho={spearman_log2_mix1:.3f}\nR2={r_squared_log2_mix1:.3f}\n\n"
        "FASTA n={fasta_count}\nLength bp: mean={fasta_mean:.1f}, median={fasta_median:.1f}\n"
        "min={fasta_min}, max={fasta_max}".format(
            n=stats_row["n"],
            pearson_log2_mix1=stats_row["pearson_log2_mix1"],
            spearman_log2_mix1=stats_row["spearman_log2_mix1"],
            r_squared_log2_mix1=stats_row["r_squared_log2_mix1"],
            fasta_count=fasta_stats["count"],
            fasta_mean=fasta_stats["mean"],
            fasta_median=fasta_stats["median"],
            fasta_min=fasta_stats["min"],
            fasta_max=fasta_stats["max"],
        ),
        transform=axis.transAxes,
        verticalalignment="top",
    )
    figure.tight_layout()
    figure.savefig(output_path, dpi=180)
    plt.close(figure)


def write_table(path, rows, fieldnames):
    directory = os.path.dirname(path)
    if directory:
        os.makedirs(directory, exist_ok=True)
    with open(path, "w") as out:
        writer = csv.DictWriter(out, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def length_pct_stats(summaries):
    values = sorted(float(row["consensus_len_pct_ref"]) for row in summaries)
    if not values:
        return {
            "metric": "consensus_len_pct_ref",
            "count": 0,
            "mean": "",
            "median": "",
            "min": "",
            "max": "",
        }

    count = len(values)
    middle = count // 2
    if count % 2:
        median = values[middle]
    else:
        median = (values[middle - 1] + values[middle]) / 2.0

    return {
        "metric": "consensus_len_pct_ref",
        "count": count,
        "mean": sum(values) / count,
        "median": median,
        "min": values[0],
        "max": values[-1],
    }


def main():
    args = parse_args()
    ercc_ref = args.ercc_ref or default_ercc_ref()
    ercc_refs = read_ercc_truth(ercc_ref)
    fasta_lengths = read_fasta_lengths(args.fasta)
    best_paf = read_best_ercc_paf(args.paf, ercc_refs) if args.paf else {}
    fasta_stats = fasta_length_stats(fasta_lengths)

    if args.fasta_length_stats:
        write_table(
            args.fasta_length_stats,
            [fasta_stats],
            ["metric", "count", "mean", "median", "min", "max"],
        )

    summaries = []
    variants = []
    missing_cs = 0
    for query, record in sorted(best_paf.items()):
        ercc_id = record["ercc_id"]
        ref_len = len(ercc_refs[ercc_id])
        consensus_len = fasta_lengths.get(query, record["query_len"])
        aligned_ref_bases = record["target_end"] - record["target_start"]
        aligned_query_bases = record["query_end"] - record["query_start"]
        nm = record["tags"].get("NM", "")
        cs = record["tags"].get("cs")

        snp_count = ""
        ins_count = ""
        del_count = ""
        if cs:
            snp_count, ins_count, del_count, record_variants = count_variants_from_cs(
                query, ercc_id, cs, record["target_start"], record["strand"]
            )
            variants.extend(record_variants)
        else:
            missing_cs += 1

        summaries.append({
            "consensus_id": query,
            "ercc_id": ercc_id,
            "ref_len": ref_len,
            "consensus_len": consensus_len,
            "consensus_len_pct_ref": 100.0 * consensus_len / ref_len if ref_len else 0,
            "aligned_ref_bases": aligned_ref_bases,
            "aligned_ref_pct": 100.0 * aligned_ref_bases / ref_len if ref_len else 0,
            "aligned_query_bases": aligned_query_bases,
            "query_aligned_pct": 100.0 * aligned_query_bases / consensus_len if consensus_len else 0,
            "matches": record["matches"],
            "block_len": record["block_len"],
            "mapq": record["mapq"],
            "strand": record["strand"],
            "target_start_1based": record["target_start"] + 1,
            "target_end_1based": record["target_end"],
            "nm": nm,
            "snp_count": snp_count,
            "insertion_count": ins_count,
            "deletion_count": del_count,
            "has_snp": "yes" if snp_count != "" and snp_count > 0 else ("unknown" if snp_count == "" else "no"),
            "has_indel": "yes" if ins_count != "" and (ins_count + del_count) > 0 else ("unknown" if ins_count == "" else "no"),
        })

    if args.paf:
        summary_fields = [
            "consensus_id", "ercc_id", "ref_len", "consensus_len", "consensus_len_pct_ref",
            "aligned_ref_bases", "aligned_ref_pct", "aligned_query_bases", "query_aligned_pct",
            "matches", "block_len", "mapq", "strand", "target_start_1based", "target_end_1based",
            "nm", "snp_count", "insertion_count", "deletion_count", "has_snp", "has_indel",
        ]
        variant_fields = [
            "consensus_id", "ercc_id", "variant_type", "ref_pos_1based",
            "ref_base", "consensus_base", "strand",
        ]
        write_table(args.summary, summaries, summary_fields)
        write_table(args.variants, variants, variant_fields)
        write_table(
            args.length_stats,
            [length_pct_stats(summaries)],
            ["metric", "count", "mean", "median", "min", "max"],
        )

        if missing_cs:
            sys.stderr.write(
                "WARNING: %s PAF records do not have cs tags; exact SNP/indel positions are unknown. "
                "Regenerate PAF with minimap2 --cs=long to enable variant calls.\n" % missing_cs
            )
        sys.stderr.write("Wrote %s ERCC consensus summaries, %s variants, and length stats.\n" % (
            len(summaries), len(variants)
        ))
    else:
        sys.stderr.write("No PAF supplied; skipped length and variant evaluation.\n")

    barcodes_by_ercc = assembled_barcodes_from_fasta(fasta_lengths)
    if not barcodes_by_ercc:
        barcodes_by_ercc = assembled_barcodes_from_paf(best_paf)
    barcode_rows, corr_stats = barcode_vs_concentration(
        barcodes_by_ercc,
        ercc_ref,
        args.barcode_vs_concentration,
        args.include_unobserved_ercc,
    )
    barcode_plot = args.barcode_vs_concentration_plot or "%s.png" % args.barcode_vs_concentration
    plot_barcode_vs_concentration(barcode_rows, corr_stats, barcode_plot, fasta_stats)
    sys.stderr.write(
        "Wrote %s barcode-vs-concentration rows to %s.tsv and %s "
        "(Pearson r=%.3f, R2=%.3f, n=%s).\n" % (
            len(barcode_rows), args.barcode_vs_concentration,
            barcode_plot,
            corr_stats["pearson_log2_mix1"], corr_stats["r_squared_log2_mix1"], corr_stats["n"],
        )
    )


if __name__ == "__main__":
    main()
