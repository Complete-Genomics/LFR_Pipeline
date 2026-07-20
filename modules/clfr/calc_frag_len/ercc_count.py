"""ERCC concentration correlation against mapped reads and mapped barcodes."""
import pandas as pd
import argparse
import numpy as np

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
except ImportError:
    plt = None

parser = argparse.ArgumentParser()
parser.add_argument("--ercc_count", type=str, required=True)
parser.add_argument("--output", type=str, required=True)
parser.add_argument("--ercc_ref", type=str, required=True, help="Path to ERCC truth/reference table")
parser.add_argument("--frag_bc_df", type=str, help="Optional frag_and_bc_dataframe.tsv with Chrom/Barcode/N_Reads")

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


def read_frag_barcode_counts(path):
    df = pd.read_csv(path, index_col=None, sep='\t')
    required = {'Chrom', 'Barcode', 'N_Reads'}
    missing = required - set(df.columns)
    if missing:
        raise ValueError("frag dataframe missing columns: {}".format(', '.join(sorted(missing))))

    df = df[df['Chrom'].astype(str).str.startswith('ERCC-')].copy()
    df['N_Reads'] = pd.to_numeric(df['N_Reads'], errors='coerce').fillna(0)
    return df.groupby('Chrom').agg(
        mapped_barcode_count=('Barcode', 'nunique'),
        mapped_fragment_count=('Barcode', 'size'),
        fragment_reads_count=('N_Reads', 'sum'),
    ).reset_index().rename(columns={'Chrom': 'ERCC ID'})


def correlation_rows(df, feature_names):
    rows = []
    x = np.log2(df['concentration in Mix 1 (attomoles/ul)'] + 1)
    for feature_name in feature_names:
        if feature_name not in df.columns:
            continue
        y = np.log2(df[feature_name].fillna(0) + 1)
        pearson = np.nan
        spearman = np.nan
        if x.nunique() > 1 and y.nunique() > 1:
            pearson = x.corr(y, method='pearson')
            spearman = x.corr(y, method='spearman')
        rows.append({
            'feature': feature_name,
            'pearson_log2_mix1': pearson,
            'spearman_log2_mix1': spearman,
        })
    return pd.DataFrame(rows)


def plot_ercc(df, feature_name, output):
    if plt is None:
        return
    plot_df = df[['concentration in Mix 1 (attomoles/ul)', feature_name]].dropna()
    if plot_df.empty:
        return

    fig, ax = plt.subplots(figsize=(8, 8))
    x = np.log2(plot_df['concentration in Mix 1 (attomoles/ul)'] + 1)
    y = np.log2(plot_df[feature_name] + 1)
    ax.scatter(x, y, s=60, alpha=0.7, edgecolors="k")
    if len(plot_df) >= 2 and x.nunique() > 1 and y.nunique() > 1:
        b, a = np.polyfit(x, y, deg=1)
        ax.plot(x, a + b * x, color="k", lw=2.5)
    ax.set_ylabel(feature_name+'_log2')
    ax.set_xlabel('mix1_concentration_log2')
    ax.set_title(feature_name+'_vs_mix1_concentration')
    plt.savefig("_".join([output.split(".txt")[0],feature_name,"png"]))
    plt.clf()


df_truth = read_ercc_truth(infile_ref)
df_counts = read_featurecounts(infile_ercc_count)
df_out = df_truth.merge(df_counts, on='ERCC ID', how='left')

if 'size' in df_out.columns:
    df_out['Length'] = df_out['Length'].fillna(pd.to_numeric(df_out['size'], errors='coerce'))
df_out['mapped_reads_count'] = df_out['mapped_reads_count'].fillna(0)
df_out['count_normalized'] = df_out['mapped_reads_count'] / df_out['Length']
df_out['ratio_mix1'] = df_out['count_normalized'] / df_out['concentration in Mix 1 (attomoles/ul)']
df_out['ratio_mix2'] = df_out['count_normalized'] / df_out['concentration in Mix 2 (attomoles/ul)']

features = ['mapped_reads_count', 'count_normalized']
if args.frag_bc_df:
    df_barcode = read_frag_barcode_counts(args.frag_bc_df)
    df_out = df_out.merge(df_barcode, on='ERCC ID', how='left')
    for col in ('mapped_barcode_count', 'mapped_fragment_count', 'fragment_reads_count'):
        df_out[col] = df_out[col].fillna(0)
    features.extend(['mapped_barcode_count', 'mapped_fragment_count', 'fragment_reads_count'])

df_out.to_csv(output, index=False, sep='\t')
correlation_rows(df_out, features).to_csv(
    output.replace('.txt', '_correlation.txt'),
    index=False,
    sep='\t',
)

for feature in features:
    plot_ercc(df_out, feature, output)
