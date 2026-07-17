"""
bin fragment and calculate mean coverage, from 5' to 3' to check if higher coverage at 3' end
flip start/end if SE, PE R1 aligned reverse
Feature 1 (bin_num10): split each frag into configured number of bins, count reads per bin as depth, boxplot
Feature 2 (bin_size100): bin by 100bp, count reads per bin as depth, boxplot
usage:
python frag_mean_cov_boxplot.py --read_len 100 --num_bin_in_frag 100 --name_sample V350096874_L01 --dirname Calc_Frag_Length_{split} --library_type gdna --ylim 100 --minreads 50 --processNline 50000
output: mean_cov_bin pickle files and png plots
"""
import numpy as np
import pandas as pd
from ast import literal_eval
import scipy.stats as st
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import subprocess
import seaborn as sns
import pickle
import argparse
import os
import traceback

parser = argparse.ArgumentParser()
parser.add_argument("--name_sample", type=str, help="name of sample1")
parser.add_argument("--library_type", type=str, help="gdna or mrna", default="gdna")
parser.add_argument("--read_len", type=int, required=True, help="100")
parser.add_argument("--num_bin_in_frag", type=int, default=100, help="100")
parser.add_argument("--dirname", type=str, help="Calc_Frag_Length_{split}")
parser.add_argument("--ylim", type=int, help="100")
parser.add_argument("--minreads", type=int, default=50, help="50")
parser.add_argument("--processNline", type=int, default=None, help="Number of lines to process from dataframe; default: all lines")

# processNline = 50000


def preprocess(dirname, read_len, num_bin_in_frag, minreads, processNline=None):
    infile = 'frag_and_bc_dataframe.tsv'
    infile_bc_df = dirname + f"/{infile}"
    read_csv_kwargs = dict(index_col=None, sep='\t')
    if processNline is not None:
        read_csv_kwargs['nrows'] = processNline
    bc_df = pd.read_csv(infile_bc_df, **read_csv_kwargs)
    bc_df.Positions = bc_df.Positions.apply(literal_eval)
    bc_df.Cigar_match = bc_df.Cigar_match.apply(literal_eval)
    bc_df = bc_df[bc_df['N_Reads'] >= minreads]
    
    # Calculate reads_covered_fragment_percent and frag_mean_depth from Positions
    def calc_gap_length_total(positions):
        total_gap = 0
        for i in range(1, len(positions)):
            gap = positions[i] - (positions[i-1] + read_len)
            if gap > 0:
                total_gap += gap
        return total_gap
    
    bc_df['gap_length_total'] = bc_df['Positions'].apply(calc_gap_length_total)
    bc_df['reads_covered_fragment_percent'] = 1 - bc_df['gap_length_total'] / bc_df['Frag_Length']
    bc_df['frag_mean_depth'] = bc_df['N_Reads'] * read_len / bc_df['Frag_Length']
    
    return bc_df


def frag_mean_cov_100bins(bc_df, read_len, num_bin_in_frag, minreads, nk):
    mean_cov_bin, top_nk_mean, tail_nk_mean = [], [], []
    slopes = []
    slope_x = [i + 1 for i in range(num_bin_in_frag)]

    for i in range(bc_df.shape[0]):
        strand = bc_df.iloc[i]['Strand']
        pos = bc_df.iloc[i]['Positions']
        frag_len = bc_df.iloc[i]['Frag_Length']
        cigar_match = bc_df.iloc[i]['Cigar_match']
        frag_start = pos[0]
        frag_end = pos[-1] + read_len

        bin_size = frag_len // num_bin_in_frag
        hist, edges = np.histogram(pos, bins=num_bin_in_frag, range=(frag_start, frag_end), density=False)

        cov = []
        for j in range(len(hist)):
            idx_start_include = sum(hist[:j])
            idx_end_notinclude = sum(hist[:j + 1])
            cov.append(sum(cigar_match[idx_start_include:idx_end_notinclude]))

        for j in range(len(hist) - 1):
            idx = sum(hist[:j + 1]) - 1
            if pos[idx] + cigar_match[idx] >= edges[j + 1]:
                cov[j + 1] += pos[idx] + cigar_match[idx] - edges[j + 1]
                cov[j] -= edges[j + 1] - pos[idx]

        mean_cov = np.array(cov) / bin_size
        if strand == 1:
            mean_cov_bin.append(mean_cov)
            slope, intercept, r_value, p_value, std_err = st.linregress(slope_x, mean_cov)
        elif strand == 0:
            mean_cov_bin.append(mean_cov[::-1])
            slope, intercept, r_value, p_value, std_err = st.linregress(slope_x, mean_cov[::-1])

        slopes.append(slope)

    return mean_cov_bin, slopes, top_nk_mean, tail_nk_mean


def frag_mean_cov_binSize(bc_df, read_len, frag_bin, bin_size, minreads, nk):
    mean_cov_bin = []

    for i in range(bc_df.shape[0]):
        strand = bc_df.iloc[i]['Strand']
        pos = bc_df.iloc[i]['Positions']
        frag_len = bc_df.iloc[i]['Frag_Length']
        cigar_match = bc_df.iloc[i]['Cigar_match']
        frag_start = pos[0]
        frag_end = pos[-1] + read_len
        num_bin_in_frag = frag_bin // bin_size
        x = [j + 1 for j in range(num_bin_in_frag)]
        hist, edges = np.histogram(pos, bins=num_bin_in_frag, range=(frag_start + 1, frag_start + bin_size * num_bin_in_frag), density=False)

        cov = []
        for j in range(len(hist)):
            idx_start_include = sum(hist[:j])
            idx_end_notinclude = sum(hist[:j + 1])
            cov.append(sum(cigar_match[idx_start_include:idx_end_notinclude]))

        for j in range(len(hist) - 1):
            idx = sum(hist[:j + 1]) - 1
            if pos[idx] + cigar_match[idx] >= edges[j + 1]:
                cov[j + 1] += pos[idx] + cigar_match[idx] - edges[j + 1]
                cov[j] -= edges[j + 1] - pos[idx]

        mean_cov = np.array(cov) / bin_size
        if strand == 1:
            mean_cov_bin.append(mean_cov)
        elif strand == 0:
            mean_cov_bin.append(mean_cov[::-1])

    return mean_cov_bin


class frag_mean_cov_boxplot(object):
    def __init__(self):
        self.df = pd.DataFrame()
        return

    def set_param(self, name_sample, dirname, ylim):
        self.name_sample = name_sample
        self.dirname = dirname
        self.ylim = ylim

    def parse_mean_cov_bin(self, frag_bin, bin_type):
        batch = self.name_sample
        _dir = self.dirname
        with open(_dir + '/{}/mean_cov_bin_{}'.format(bin_type, str(frag_bin)), 'rb') as fp:
            mean_cov_bin = pickle.load(fp)
        if len(mean_cov_bin) > 0:
            num_bin_in_frag = len(mean_cov_bin[0])

        if len(mean_cov_bin) == 0:
            bin_size = bin_frag_by
            num_bin_in_frag = frag_bin // bin_size
            mean_cov_bin = np.array([np.zeros(num_bin_in_frag)]) + 1

        _df = pd.DataFrame(mean_cov_bin)
        _df.columns = ['bin' + str(i) for i in range(1, num_bin_in_frag + 1)]
        df = _df.stack().reset_index()
        df.columns = ['idx', 'fragment_length', 'mean_coverage']
        df['sample'] = self.name_sample

        return df

    def box_plot_bin(self, df_combined_log, frag_bin, bin_type):
        sns.set(rc={'figure.figsize': (10, 12)})
        sns.set_style("white")
        sns.set(context='notebook', style='ticks', font_scale=1.8)
        ax = sns.boxplot(x='fragment_length', y='mean_coverage', hue='sample', palette='tab10', showfliers=False, data=df_combined_log)

        ax.set_xlabel('Fragment length: ' + str(frag_bin), {'weight': 'bold'})
        ax.set_ylabel('Mean coverage', {'weight': 'bold'})
        num_bin = df_combined_log['fragment_length'].nunique()
        num_frag = df_combined_log.shape[0] // max(num_bin, 1)
        ax.set_title('num_frag={}'.format(num_frag))
        ax.legend(loc="upper right")
        ax.set(ylim=(0, self.ylim))
        fig = ax.get_figure()
        fig.savefig(self.dirname + "/{}/mean_cov_bin_boxplot_{}.png".format(bin_type, frag_bin))
        fig.clear()

    def line_plot_bin10(self, df, frag_bin, bin_type):
        sns.set(rc={'figure.figsize': (10, 12)})
        sns.set_style("white")
        sns.set(context='notebook', style='ticks', font_scale=1.8)

        mean_by_bin = df.groupby('fragment_length')['mean_coverage'].mean().reset_index()
        mean_by_bin = mean_by_bin.sort_values('fragment_length')

        ax = plt.figure().add_subplot(111)
        ax.plot(mean_by_bin['fragment_length'], mean_by_bin['mean_coverage'],
                marker='o', linewidth=2, markersize=4, label=self.name_sample)

        ax.set_xlabel('Position in fragment (bp)', {'weight': 'bold'})
        ax.set_ylabel('Mean coverage', {'weight': 'bold'})
        ax.set_title('Fragment length: {}bp, bin size: 10bp'.format(frag_bin))
        ax.legend(loc="upper right")
        ax.set(ylim=(0, self.ylim))
        ax.grid(True, alpha=0.3)

        fig = ax.get_figure()
        fig.savefig(self.dirname + "/{}/mean_cov_line_bin10_{}.png".format(bin_type, frag_bin))
        fig.clear()


def _boxplot(_list, dirname, data_bin_type):
    c = pd.DataFrame()
    for i in range(len(_list)):
        data = _list[i]
        if len(data) == 0:
            continue
        tmp = pd.DataFrame(data)
        tmp.columns = [data_bin_type]
        tmp['frag_len'] = str(1000 + i * 100)
        c = pd.concat([c, tmp], ignore_index=True)

    if c.empty:
        return

    sns.set(rc={'figure.figsize': (10, 12)})
    sns.set_style("white")
    sns.set(context='notebook', style='ticks', font_scale=1.8)
    ax = sns.boxplot(x='frag_len', y=data_bin_type, palette='tab10', showfliers=False, data=c)

    ax.set_xlabel('fragment_length', {'weight': 'bold'})
    ax.set_ylabel(data_bin_type, {'weight': 'bold'})
    ax.set_title('{}_by_fragment_length'.format(data_bin_type))
    fig = ax.get_figure()
    plt.xticks(rotation=45)
    fig.savefig(dirname + "/{}.png".format(data_bin_type))
    fig.clear()


def main():
    args = parser.parse_args()
    name_sample = args.name_sample
    library_type = args.library_type
    dirname = args.dirname
    ylim = args.ylim
    minreads = args.minreads
    read_len = args.read_len
    num_bin_in_frag = args.num_bin_in_frag
    processNline = args.processNline

    if library_type == 'gdna':
        frag_len_max = 10000
        frag_len_min = 1000
    elif library_type == 'mrna':
        frag_len_max = 3000
        frag_len_min = 200

    global bin_frag_by
    bin_frag_by = 100

    # Create output directories
    os.makedirs(dirname + '/bin_num10', exist_ok=True)
    os.makedirs(dirname + '/bin_size100', exist_ok=True)

    # Calculate mean coverage
    bc_df = preprocess(dirname, read_len, num_bin_in_frag, minreads, processNline)
    print(f"bc_df shape: {bc_df.shape}")
    print(f"bc_df columns: {list(bc_df.columns)}")
    print(f"Sample Strand values: {bc_df['Strand'].head() if 'Strand' in bc_df.columns else 'NO Strand column'}")
    _reads_covered_fragment_percent, _frag_mean_depth = [], []

    n = int(1 + (frag_len_max - frag_len_min) / bin_frag_by)
    nk = 1000

    for i in range(n):
        try:
            frag_bin = frag_len_min + bin_frag_by * i
            _df = bc_df[bc_df['Frag_Length'].between(frag_len_min + bin_frag_by * i, frag_len_min + bin_frag_by * (i + 1))]
            print(f"frag_bin={frag_bin}, num_frags={_df.shape[0]}")
            if _df.shape[0] == 0:
                continue
            mean_cov_bin_100bins, slopes, top_nk_mean, tail_nk_mean = frag_mean_cov_100bins(_df, read_len, num_bin_in_frag, minreads, nk)
            _reads_covered_fragment_percent.append(_df['reads_covered_fragment_percent'])
            _frag_mean_depth.append((_df['frag_mean_depth']))

            with open(dirname + '/bin_num10/mean_cov_bin_' + str(frag_bin), 'wb') as fp:
                pickle.dump(mean_cov_bin_100bins, fp)

            # bin_size100
            bin_size = 100
            mean_cov_binSize = frag_mean_cov_binSize(_df, read_len, frag_bin, bin_size, minreads, nk)

            with open(dirname + '/bin_size100/mean_cov_bin_' + str(frag_bin), 'wb') as fp:
                pickle.dump(mean_cov_binSize, fp)
        except Exception as e:
            print(f"Error at frag_bin={frag_bin}: {e}")
            traceback.print_exc()

    type1 = 'bin_num10'
    with open(dirname + '/' + type1 + '/mean_cov_bin_' + str(minreads) + '_done', 'w') as f:
        f.write("done\n")
    type2 = 'bin_size100'
    with open(dirname + '/' + type2 + '/mean_cov_bin_' + str(minreads) + '_done', 'w') as f:
        f.write("done\n")

    with open(dirname + '_reads_covered_fragment_percent', 'wb') as fp:
        pickle.dump(_reads_covered_fragment_percent, fp)
    with open(dirname + '_frag_mean_depth', 'wb') as fp:
        pickle.dump(_frag_mean_depth, fp)

    # Plotting
    s = frag_mean_cov_boxplot()
    s.set_param(name_sample, dirname, ylim)

    for i in range(n):
        frag_bin = frag_len_min + bin_frag_by * i
        try:
            bin_type = 'bin_num10'
            mean_cov_bin_combined_log = s.parse_mean_cov_bin(frag_bin, bin_type)
            s.box_plot_bin(mean_cov_bin_combined_log, frag_bin, bin_type)
            s.line_plot_bin10(mean_cov_bin_combined_log, frag_bin, bin_type)
        except Exception as e:
            print(e)
        try:
            bin_type = 'bin_size100'
            mean_cov_bin_combined_log = s.parse_mean_cov_bin(frag_bin, bin_type)
            s.box_plot_bin(mean_cov_bin_combined_log, frag_bin, bin_type)
            s.line_plot_bin10(mean_cov_bin_combined_log, frag_bin, bin_type)
        except Exception as e:
            print(e)

    with open(dirname + '_'.join(['mean_cov_bin_boxplot', str(minreads), 'done']), 'w') as f:
        f.write("done\n")

    # Fragment metrics plots
    with open(dirname + '_reads_covered_fragment_percent', 'rb') as fp:
        _reads_covered_fragment_percent = pickle.load(fp)
    _boxplot(_reads_covered_fragment_percent, dirname, data_bin_type='reads_covered_fragment_percent')

    with open(dirname + '_frag_mean_depth', 'rb') as fp:
        _frag_mean_depth = pickle.load(fp)
    _boxplot(_frag_mean_depth, dirname, data_bin_type='frag_mean_depth')


if __name__ == "__main__":
    main()
