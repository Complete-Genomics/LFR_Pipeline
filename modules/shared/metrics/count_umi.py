#!/usr/bin/env python3
# import argparse

import pandas as pd


# parser = argparse.ArgumentParser()
# parser.add_argument("-i", "--input", required=True, help="3-column TSV: index, read_count, umi")
# parser.add_argument("-o", "--prefix", required=True)
# args = parser.parse_args()

input = 'split_stat_read1.log'
# prefix = 'split_log'
cutoff = 25


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
    return pd.DataFrame([
        {
            "bin": label,
            "n_reads": per_umi.loc[mask, "read_count"].sum(),
            "n_umi": mask.sum(),
        }
        for label, mask in bins
    ])


df = pd.read_csv(
    input,
    sep="\t",
    header=None, skiprows=4, 
    names=["index", "read_count", "umi"],
    dtype={"index": str, "umi": str},
)

df["read_count"] = pd.to_numeric(df["read_count"], errors="raise")
per_umi = df
# 同一 (index, umi) 多行时，合并其 read count
# per_umi = (
#     df.groupby(["index", "umi"], as_index=False)["read_count"]
#       .sum()
# )
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

stats.to_csv(f"N{cutoff}.per_umi_read_count.stats.tsv", sep="\t", index=False)
bin_umi(per_umi).to_csv(
    f"N{cutoff}.per_umi_read_count.bins.tsv", sep="\t", index=False
)
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
