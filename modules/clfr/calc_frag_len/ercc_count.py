"""ERCC concentration correlation against mapped reads and mapped barcodes."""
import pandas as pd
import argparse
import numpy as np
import os
from collections import defaultdict

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
except ImportError:
    plt = None

parser = argparse.ArgumentParser()
parser.add_argument("--ercc_count", type=str)
parser.add_argument("--ercc_bam", type=str)
parser.add_argument("--output", type=str, required=True)
parser.add_argument("--ercc_ref", type=str, required=True, help="Path to ERCC truth/reference table")
parser.add_argument("--mapq", type=int, default=0)

args = parser.parse_args()
infile_ercc_count = args.ercc_count
output = args.output
infile_ref = args.ercc_ref


def read_ercc_truth(path):
    df = pd.read_csv(path, index_col=None, sep='\t')
    if 'ERCC ID' not in df.columns:
        raise ValueError("ERCC reference must contain column: ERCC ID")
    if 'concentration in Mix 1 (attomoles/ul)' not in df.columns:
        raise ValueError("ERCC reference must contain Mix 1 concentration column")
    if 'concentration in Mix 2 (attomoles/ul)' not in df.columns:
        raise ValueError("ERCC reference must contain Mix 2 concentration column")
    return df


def read_featurecounts(path):
    df = pd.read_csv(path, index_col=None, sep='\t', comment='#')
    if 'Geneid' not in df.columns:
        raise ValueError("featureCounts output must contain column: Geneid")
    if 'Length' not in df.columns:
        raise ValueError("featureCounts output must contain column: Length")

    count_col = df.columns[-1]
    out = df[['Geneid', 'Length', count_col]].copy()
    out.columns = ['ERCC ID', 'Length', 'mapped_reads_count']
    out['mapped_reads_count'] = pd.to_numeric(out['mapped_reads_count'], errors='coerce').fillna(0)
    out['Length'] = pd.to_numeric(out['Length'], errors='coerce')
    return out


def barcode_from_read(read):
    try:
        return read.get_tag('BX')
    except KeyError:
        pass

    name = read.query_name
    if '#' not in name:
        return None
    barcode = name.split('#', 1)[1]
    barcode = barcode.split('#', 1)[0]
    barcode = barcode.split('/', 1)[0]
    return barcode or None


def read_ercc_bam_counts(path, mapq):
    try:
        import pysam
    except ImportError:
        raise ImportError("pysam is required when using --ercc_bam")

    read_counts = defaultdict(int)
    barcode_sets = defaultdict(set)

    with pysam.AlignmentFile(path, 'rb') as bam:
        for read in bam.fetch(until_eof=True):
            if read.is_unmapped:
                continue
            if read.is_secondary or read.is_supplementary:
                continue
            if read.mapping_quality < mapq:
                continue
            ercc_id = bam.get_reference_name(read.reference_id)
            if not ercc_id or not ercc_id.startswith('ERCC-'):
                continue
            read_counts[ercc_id] += 1
            barcode = barcode_from_read(read)
            if barcode:
                barcode_sets[ercc_id].add(barcode)

    ercc_ids = sorted(set(read_counts) | set(barcode_sets))
    return pd.DataFrame({
        'ERCC ID': ercc_ids,
        'mapped_reads_count': [read_counts[ercc_id] for ercc_id in ercc_ids],
        'mapped_barcode_count': [len(barcode_sets[ercc_id]) for ercc_id in ercc_ids],
    })


def log2p(series):
    return np.log2(pd.to_numeric(series, errors='coerce').fillna(0) + 1)


def correlation_stats(x, y):
    stats = {
        'n': int(len(x)),
        'pearson': np.nan,
        'spearman': np.nan,
        'r_squared': np.nan,
        'slope': np.nan,
        'intercept': np.nan,
    }
    if len(x) >= 2 and x.nunique() > 1 and y.nunique() > 1:
        stats['pearson'] = x.corr(y, method='pearson')
        stats['spearman'] = x.corr(y, method='spearman')
        stats['r_squared'] = stats['pearson'] ** 2
        stats['slope'], stats['intercept'] = np.polyfit(x, y, deg=1)
    return stats


def correlation_rows(df, feature_names):
    rows = []
    x = log2p(df['concentration in Mix 1 (attomoles/ul)'])
    for feature_name in feature_names:
        if feature_name not in df.columns:
            continue
        y = log2p(df[feature_name])
        stats = correlation_stats(x, y)
        rows.append({
            'feature': feature_name,
            'n': stats['n'],
            'pearson_log2_mix1': stats['pearson'],
            'spearman_log2_mix1': stats['spearman'],
            'r_squared_log2_mix1': stats['r_squared'],
            'slope_log2_mix1': stats['slope'],
            'intercept_log2_mix1': stats['intercept'],
        })
    return pd.DataFrame(rows)


def output_base(output):
    base, ext = os.path.splitext(output)
    return base if ext else output


def plot_path(output, feature_name):
    return "%s_%s.png" % (output_base(output), feature_name)


def plot_grouped_points(ax, plot_df, x, y):
    if 'subgroup' not in plot_df.columns:
        ax.scatter(x, y, s=46, alpha=0.78, edgecolors="k", linewidths=0.4)
        return

    colors = {
        'A': '#1f77b4',
        'B': '#2ca02c',
        'C': '#ff7f0e',
        'D': '#d62728',
    }
    for subgroup in sorted(plot_df['subgroup'].dropna().unique()):
        mask = plot_df['subgroup'] == subgroup
        ax.scatter(
            x[mask],
            y[mask],
            label="subgroup %s" % subgroup,
            color=colors.get(subgroup, '#7f7f7f'),
            s=46,
            alpha=0.78,
            edgecolors="k",
            linewidths=0.4,
        )


def format_feature_label(feature_name):
    labels = {
        'mapped_reads_count': 'mapped reads count',
        'mapped_barcode_count': 'mapped barcode count',
        'count_normalized': 'length-normalized reads count',
    }
    return labels.get(feature_name, feature_name.replace('_', ' '))


def draw_ercc_panel(ax, df, feature_name):
    keep_cols = ['concentration in Mix 1 (attomoles/ul)', feature_name]
    if 'subgroup' in df.columns:
        keep_cols.append('subgroup')
    plot_df = df[keep_cols].copy()
    plot_df = plot_df.dropna(subset=['concentration in Mix 1 (attomoles/ul)', feature_name])
    if plot_df.empty:
        ax.set_axis_off()
        return

    x = log2p(plot_df['concentration in Mix 1 (attomoles/ul)'])
    y = log2p(plot_df[feature_name])
    stats = correlation_stats(x, y)

    plot_grouped_points(ax, plot_df, x, y)
    if not np.isnan(stats['slope']):
        order = np.argsort(x.values)
        x_sorted = x.values[order]
        y_fit = stats['intercept'] + stats['slope'] * x_sorted
        ax.plot(x_sorted, y_fit, color="black", lw=2.0)

    label = format_feature_label(feature_name)
    ax.set_ylabel("log2(%s + 1)" % label)
    ax.set_xlabel("log2(ERCC Mix 1 concentration + 1)")
    ax.set_title("%s vs ERCC concentration" % label)
    ax.grid(True, color="#dddddd", linewidth=0.6, alpha=0.8)

    stat_label = "n=%s\nPearson r=%.3f\nSpearman rho=%.3f\nR2=%.3f" % (
        stats['n'],
        stats['pearson'] if not np.isnan(stats['pearson']) else float('nan'),
        stats['spearman'] if not np.isnan(stats['spearman']) else float('nan'),
        stats['r_squared'] if not np.isnan(stats['r_squared']) else float('nan'),
    )
    ax.text(
        0.04,
        0.96,
        stat_label,
        transform=ax.transAxes,
        ha='left',
        va='top',
        fontsize=9,
        bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='#cccccc', alpha=0.9),
    )


def plot_ercc(df, feature_name, output):
    if plt is None:
        return

    fig, ax = plt.subplots(figsize=(8, 7))
    draw_ercc_panel(ax, df, feature_name)
    if 'subgroup' in df.columns:
        ax.legend(frameon=False, loc='lower right', fontsize=9)
    fig.tight_layout()
    fig.savefig(plot_path(output, feature_name), dpi=180)
    plt.close(fig)


def plot_ercc_summary(df, feature_names, output):
    if plt is None:
        return
    if not feature_names:
        return

    ncols = min(2, len(feature_names))
    nrows = int(np.ceil(float(len(feature_names)) / ncols))
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(8 * ncols, 6.5 * nrows))
    axes = np.array(axes).reshape(-1)
    for ax, feature_name in zip(axes, feature_names):
        draw_ercc_panel(ax, df, feature_name)
    for ax in axes[len(feature_names):]:
        ax.set_axis_off()
    if 'subgroup' in df.columns:
        handles, labels = axes[0].get_legend_handles_labels()
        if handles:
            fig.legend(handles, labels, frameon=False, loc='upper center', ncol=len(handles))
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig("%s_summary.png" % output_base(output), dpi=180)
    plt.close(fig)


df_truth = read_ercc_truth(infile_ref)
if args.ercc_bam:
    df_counts = read_ercc_bam_counts(args.ercc_bam, args.mapq)
elif infile_ercc_count:
    df_counts = read_featurecounts(infile_ercc_count)
else:
    raise ValueError("Either --ercc_bam or --ercc_count is required")
df_out = df_truth.merge(df_counts, on='ERCC ID', how='left')

if 'Length' not in df_out.columns:
    df_out['Length'] = np.nan
if 'size' in df_out.columns:
    df_out['Length'] = df_out['Length'].fillna(pd.to_numeric(df_out['size'], errors='coerce'))
if 'mapped_reads_count' not in df_out.columns:
    df_out['mapped_reads_count'] = 0
df_out['mapped_reads_count'] = df_out['mapped_reads_count'].fillna(0)
df_out['Length'] = pd.to_numeric(df_out['Length'], errors='coerce')
df_out['count_normalized'] = df_out['mapped_reads_count'] / df_out['Length'].replace(0, np.nan)
df_out['ratio_mix1'] = df_out['count_normalized'] / df_out['concentration in Mix 1 (attomoles/ul)']
df_out['ratio_mix2'] = df_out['count_normalized'] / df_out['concentration in Mix 2 (attomoles/ul)']

features = ['mapped_reads_count', 'count_normalized']
if 'mapped_barcode_count' in df_out.columns:
    df_out['mapped_barcode_count'] = df_out['mapped_barcode_count'].fillna(0)
    features.append('mapped_barcode_count')

df_out.to_csv(output, index=False, sep='\t')
correlation_rows(df_out, features).to_csv(
    output.replace('.txt', '_correlation.txt'),
    index=False,
    sep='\t',
)

for feature in features:
    plot_ercc(df_out, feature, output)
plot_ercc_summary(df_out, features, output)
