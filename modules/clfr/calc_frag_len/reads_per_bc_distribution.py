"""
a patch to generate reads_per_bc_distribution.tsv from frags_reads_per_bc.tsv
usage:
python reads_per_bc_distribution.py --dirname <analysis_dir>
"""
import pandas as pd
import seaborn as sns
import matplotlib
# matplotlib.use('Agg')
from collections import Counter
import argparse
pd.set_option('display.precision', 4)

# def reads_per_bc_distribution_mapped(dirname):
#     in_file = dirname + f'/Calc_Frag_Length_{split_distance}/frags_reads_per_bc.tsv'
#     barcode_summary = pd.read_csv(in_file, index_col=None, sep='\t')

#     # df = pd.DataFrame.from_dict(Counter(barcode_summary.Reads), orient='index').reset_index()
#     df = barcode_summary.groupby('Reads')['Barcode'].apply(list).reset_index()
#     df['bc_count'] = df['Barcode'].apply(lambda x: len(x))
#     df.columns = ['reads_per_bc', 'bc_name', 'bc_count']
#     df_sorted = df.sort_values('reads_per_bc')
#     df_sorted['reads_count'] = df_sorted.bc_count*df_sorted.reads_per_bc
#     total_reads = sum(df_sorted.reads_count)
#     df_sorted['reads_count_pct'] =df_sorted['reads_count']/total_reads
#     with open(dirname+f'/Calc_Frag_Length_{split_distance}/additional_summary.txt', 'a') as f:
#         print(f"total_reads_mapped:\t{total_reads}",file = f)
#     df_sorted.to_csv(dirname + f"/Calc_Frag_Length_{split_distance}/reads_per_bc_distribution.tsv", sep='\t', index=False, float_format='%.7f')
#     return

def reads_per_bc_distribution_regardless_mapping(dirname):
    split_log = dirname +'/split_stat_read1.log'
    summary = pd.read_csv(split_log, header=None, index_col=None, sep='\t', skiprows=4)        
    summary.columns = ['index', 'Reads', 'BC']

    df = pd.DataFrame.from_dict(Counter(summary.Reads), orient='index').reset_index()
    df.columns = ['reads_per_bc','bc_count']
    df_sorted = df.sort_values('reads_per_bc')
    df_sorted['reads_count'] = df_sorted.bc_count*df_sorted.reads_per_bc
    total_reads = sum(df_sorted.reads_count)
    df_sorted['reads_count_pct'] =df_sorted['reads_count']/total_reads
    with open(dirname+f'/Calc_Frag_Length_{split_distance}/additional_summary.txt', 'a') as f:
        print(f"total_reads_regardless_mapping:\t{total_reads}",file = f)
    df_sorted.to_csv(dirname + f"/Calc_Frag_Length_{split_distance}/reads_per_bc_distribution_regardless_mapping.tsv", sep='\t', index=False, float_format='%.7f')
    return

def frag_length_distribution(dirname):
    infile_bc_df = dirname + f'/Calc_Frag_Length_{split_distance}/frag_and_bc_dataframe.tsv'
    bc_df = pd.read_csv(infile_bc_df, index_col=None, sep='\t')

    sns.set(rc={'figure.figsize':(15,8)})
    frag_plot = sns.distplot([i for i in bc_df['Frag_Length'] if i<=cutoff])
    plt.savefig(dirname + "/frag_length_distribution.pdf")
    plt.clf()
    return


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--dirname", type=str)
    args = parser.parse_args() 
    dirname = args.dirname

    # reads_per_bc_distribution_mapped(dirname)
    reads_per_bc_distribution_regardless_mapping(dirname)
    # frag_length_distribution(dirname)

if __name__ == "__main__":
    split_distance=50000
    cutoff=7000
    main()  