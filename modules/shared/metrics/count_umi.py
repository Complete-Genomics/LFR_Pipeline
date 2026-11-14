#!/usr/bin/env python3
# import argparse

import pandas as pd


# parser = argparse.ArgumentParser()
# parser.add_argument("-i", "--input", required=True, help="3-column TSV: index, read_count, umi")
# parser.add_argument("-o", "--prefix", required=True)
# args = parser.parse_args()

import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument("logs", nargs="*", default=['split_stat_read1.log'],
                    help="1 log (split_stat_read1.log) or 2 logs when a sample is in 2 batches "
                         "(duplicate (index, umi) are summed)")
parser.add_argument("--cutoff", type=int, default=25)
parser.add_argument("--outdir", default=".")
args = parser.parse_args()
missing = [f for f in args.logs if not os.path.exists(f)]
if missing:
    parser.error("log file(s) not found: " + ", ".join(missing))
if len(args.logs) not in (1, 2):
    parser.error("expected 1 or 2 log files, got %d" % len(args.logs))
# prefix = 'split_log'
inputs = args.logs
cutoff = args.cutoff


def bin_umi(per_umi):
    """Summarize reads and UMIs by reads per UMI."""
    bins = [
        ("N=1", per_umi["read_count"] == 1),
        ("N=2", per_umi["read_count"] == 2),
        ("N=3", per_umi["read_count"] == 3),
        ("N=4", per_umi["read_count"] == 4),
        ("5<=N<10", per_umi["read_count"].between(5, 9)),
        ("10<=N<15", per_umi["read_count"].between(10, 14)),
        ("15<=N<20", per_umi["read_count"].between(15, 19)),
        ("20<=N<25", per_umi["read_count"].between(20, 24)),
        ("N>=25", per_umi["read_count"] >= 25),
    ]
    summary = pd.DataFrame([
        {
            "bin": label,
            "n_reads": per_umi.loc[mask, "read_count"].sum(),
            "n_umi": mask.sum(),
        }
        for label, mask in bins
    ])
    summary["%reads"] = (summary["n_reads"] / per_umi["read_count"].sum() * 100).round(2)
    summary["%umi"] = (summary["n_umi"] / len(per_umi) * 100).round(2)
    return summary


def read_log(f):
    """Read one split_stat_read1.log: 4 header lines, then '<index>\t<reads>\t<UMI>'."""
    df = pd.read_csv(
        f,
        sep="\t",
        header=None, skiprows=4,
        names=["index", "read_count", "umi"],
        dtype={"index": str, "umi": str},
    )
    df["read_count"] = pd.to_numeric(df["read_count"], errors="raise")
    return df


def write_summary(per_umi, cutoff, outdir):
    counts = per_umi.loc[per_umi["read_count"] >= cutoff, "read_count"]
    stats = pd.DataFrame([{
        "n_umi": counts.count(),
        "mean": round(counts.mean(), 2),
        "median": counts.median(),
        "min": counts.min(),
        "max": counts.max(),
        "std": round(counts.std(), 2),
        "cutoff": cutoff,
        # "n_umi_gt_1000": (counts > 1000).sum(),
    }])
    stats.to_csv(f"{outdir}/N{cutoff}.per_umi_read_count.stats.tsv", sep="\t", index=False)
    bin_umi(per_umi).to_csv(
        f"{outdir}/N{cutoff}.per_umi_read_count.bins.tsv", sep="\t", index=False
    )


def count_umi_1log(log, cutoff, outdir):
    """Single batch: split_stat_read1.log."""
    write_summary(read_log(log), cutoff, outdir)


def count_umi_2logs(log1, log2, cutoff, outdir):
    """Sample split in 2 batches: duplicate (index, umi) across the two logs get their read counts summed."""
    df = pd.concat([read_log(log1), read_log(log2)], ignore_index=True)
    per_umi = df.groupby(["index", "umi"], as_index=False)["read_count"].sum()
    write_summary(per_umi, cutoff, outdir)


if len(inputs) == 1:
    count_umi_1log(inputs[0], cutoff, args.outdir)
else:
    count_umi_2logs(inputs[0], inputs[1], cutoff, args.outdir)
# per_umi.to_csv(f"{prefix}.per_umi_read_count.tsv", sep="\t", index=False)

# plot_counts = per_umi.loc[
#     per_umi["read_count"].between(1, 1000)
# ].copy()

# sns.set_theme(style="whitegrid")
# fig, ax = plt.subplots(figsize=(7, 5))

# sns.histplot(
#     data=per_umi,
#     x="read_count",
#     bins=100,
#     stat="density",
#     color="#4C78A8",
#     edgecolor="white",
#     ax=ax,
# )

# ax.set_xlim(1, 1000)
# ax.set_xlabel("Read count per (index, UMI)")
# ax.set_ylabel("Density")
# ax.set_title("Per-UMI read-count distribution (1–1000)")
# fig.tight_layout()
# fig.savefig(f"{prefix}.per_umi_read_count.density_1_1000.png", dpi=300)

# top_counts = (
#     per_umi.nlargest(100_000, "read_count")["read_count"]
# )

# stats = pd.DataFrame([{
#     "n_umi": top_counts.count(),
#     # "n_umi_gt_1000": (top_counts > 1000).sum(),
#     "mean": round(top_counts.mean(), 2),
#     "median": top_counts.median(),
#     "std": round(top_counts.std(), 2),
#     "min": top_counts.min(),
#     "max": top_counts.max(),
# }])

# stats.to_csv(
#     f"{prefix}.top100k_per_umi_read_count.stats.tsv",
#     sep="\t",
#     index=False,
# )
