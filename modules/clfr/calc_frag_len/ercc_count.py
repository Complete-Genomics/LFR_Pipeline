"""
raw ercc count only vs mix1
"""
import pandas as pd
from Bio import SeqIO
import pickle
import argparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("--ercc_count", type=str)
# parser.add_argument("--ercc_frag", type=str)
parser.add_argument("--output", type=str)
parser.add_argument("--ercc_ref", type=str, help="Path to ERCC reference GC file (ercc_ref_gc.txt)")

args = parser.parse_args() 
infile_ercc_count = args.ercc_count
# infile_ercc_frag = args.ercc_frag
output = args.output
infile_ref = args.ercc_ref
### sample ercc count
# infile_ercc_count='ercc_count.txt'
df =pd.read_csv(infile_ercc_count, index_col=None, sep='\t')
# df3 = pd.read_csv(infile_ercc_frag, index_col=None, sep='\t')
# df4 = df.groupby('Chrom')['N_Reads'].agg(['sum', 'count']).reset_index().rename(columns={'sum':'Reads', 'count':'Frags'})

### ref ercc
df2 =pd.read_csv(infile_ref, index_col=None, sep='\t')
df2['raw_count']=df.iloc[:,-1]
df2['Length'] = df['Length']
df2['count_normalized'] = df2['raw_count']/df2['Length']
df2['ratio_mix1']= df2['count_normalized']/df2['concentration in Mix 1 (attomoles/ul)']
df2['ratio_mix2']= df2['count_normalized']/df2['concentration in Mix 2 (attomoles/ul)']
# df2['frag_cnt'] = df4['Frags']
df2.to_csv(output, index=False, sep='\t')

### plot mix1 vs. count_normalized
def plot_ercc(df2, feature_name, output):
    fig, ax = plt.subplots(figsize = (8, 8))
    x=np.log2(df2['concentration in Mix 1 (attomoles/ul)']+1)
    y=np.log2(df2[feature_name]+1)
    ax.scatter(x, y, s=60, alpha=0.7, edgecolors="k")
    b, a = np.polyfit(x,y , deg=1)
    ax.plot(x, a + b * x, color="k", lw=2.5)
    ax.set_ylabel(feature_name+'_log2')
    ax.set_xlabel('mix1_concentration_log2')
    ax.set_title(feature_name+'_vs_mix1_concentration')
    plt.savefig("_".join([output.split(".txt")[0],feature_name,"png"]))
    plt.clf()

## correlation, against read counts
plot_ercc(df2, 'count_normalized', output)

## correlation, against frag counts
# plot_ercc(df2, 'ercc_frag', output)