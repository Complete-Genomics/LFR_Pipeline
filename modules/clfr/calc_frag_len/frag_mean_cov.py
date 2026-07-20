"""
for .smk
bin fragment and calculate mean coverage, from 5' to 3' to check if higher coverage at 3' end
flip start/end if SE, PE R1 aligned reverse
split each frag into 10 bins
usage:
python frag_mean_cov.py --read_len 100 --num_bin_in_frag 10 --name_sample V350096874_L01 --bin_frag_by 100
output: mean_cov_bin, slopes
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
parser = argparse.ArgumentParser()
# parser.add_argument("--frag_bc_df", type=str, help="frag_and_bc_dataframe.tsv")
parser.add_argument("--library_type", type=str, help="gdna or mrna", default="gdna")
parser.add_argument("--read_len", type=int, help="100")
parser.add_argument("--num_bin_in_frag", type=int, help="10")
parser.add_argument("--dirname", type=str, help="Calc_Frag_Length_{split}")
parser.add_argument("--minreads", type=int, default=4,
                    help="Minimum number of reads to keep a fragment [default: 4].")
# parser.add_argument("--frag_len_min", type=int, default=1000,help="select the min length of bin fragment length")
# parser.add_argument("--frag_len_max", type=int, default=5100,help="select the max length of bin fragment length")
# parser.add_argument("--bin_frag_by", type=int, default=100,help="the step from frag_len_min to frag_len_max")


def preprocess(dirname, read_len, num_bin_in_frag, minreads):
    # batch, lane = name_sample.split('_')
    infile='frag_and_bc_dataframe.tsv'
    # infile='test.tsv'
    infile_bc_df = dirname+f"/{infile}"
    bc_df = pd.read_csv(infile_bc_df, index_col=None, sep='\t')
    bc_df.Positions=bc_df.Positions.apply(literal_eval)
    # bc_df.Strand=bc_df.Strand.apply(literal_eval)
    bc_df.Cigar_match = bc_df.Cigar_match.apply(literal_eval)
    if 'Strand' not in bc_df.columns:
        bc_df['Strand'] = 1
    bc_df = bc_df[bc_df['N_Reads']>=minreads]
    return bc_df


def frag_mean_cov_10bins(bc_df, read_len, num_bin_in_frag, minreads, nk):
    
    mean_cov_bin, top_nk_mean, tail_nk_mean = [], [], []
    slopes = []
    slope_x = [i+1 for i in range(num_bin_in_frag)]

    for i in range(bc_df.shape[0]):       
        strand = bc_df.iloc[i]['Strand']
        pos = bc_df.iloc[i]['Positions']
        frag_len = bc_df.iloc[i]['Frag_Length']
        cigar_match = bc_df.iloc[i]['Cigar_match']
        # for SE reads, should be all in same strand
        # if len(strand)==1: [0,1] either 0,1
        frag_start = pos[0]
        frag_end = pos[-1]+read_len

        # top_nk_mean.append(len([i for i in pos if i<=pos[0]+nk])*read_len/nk)
        # tail_nk_mean.append(len([i for i in pos if i>=pos[-1]-nk])*read_len/nk)

        bin_size = frag_len//num_bin_in_frag
        ## hist[5,8,3,4,2], 5 frag in bin1
        hist, edges = np.histogram(pos,bins=num_bin_in_frag,range=(frag_start, frag_end),density=False)

        ## cov per bin
        cov = []
        for i in range(len(hist)):
            idx_start_include = sum(hist[:i])
            idx_end_notinclude = sum(hist[:i+1])
            cov.append(sum(cigar_match[idx_start_include: idx_end_notinclude]))

        ## check read spaning 2 bins
        for i in range(len(hist)-1):
            idx=sum(hist[:i+1])-1

            if pos[idx] +cigar_match[idx]>=edges[i+1]:
                cov[i+1]+=pos[idx] +cigar_match[idx]-edges[i+1]
                cov[i]-=edges[i+1]-pos[idx]

        mean_cov = np.array(cov)/bin_size
        ## flip frag start/end if in - strand, 0 means - strand, see calc_frag_len_optimized.py
        if strand ==1:
            mean_cov_bin.append(mean_cov)
            slope, intercept, r_value, p_value, std_err = st.linregress(slope_x, mean_cov)
        elif strand==0:
            mean_cov_bin.append(mean_cov[::-1])
            slope, intercept, r_value, p_value, std_err = st.linregress(slope_x, mean_cov[::-1])

        slopes.append(slope)
        ## if both + and - strand present in a frag
        # else:
        #     print(i)
    ## mean_cov_bin=for all frag with len x, [frag1_list, frag2_list,..., fragn], frag1_list=[bin1_mean_cov_bin, bin2_mean_cov_bin, ..., bin10_mean_cov_bin]

    return mean_cov_bin, slopes, top_nk_mean, tail_nk_mean

def frag_mean_cov_binSize100(bc_df, read_len, frag_bin, bin_size, minreads, nk):
    
    mean_cov_bin = []
    
    
    for i in range(bc_df.shape[0]):       
        strand = bc_df.iloc[i]['Strand']
        pos = bc_df.iloc[i]['Positions']
        frag_len = bc_df.iloc[i]['Frag_Length']
        cigar_match = bc_df.iloc[i]['Cigar_match']
        frag_start = pos[0]
        frag_end = pos[-1]+read_len
        num_bin_in_frag  = frag_bin//bin_size
        x = [i+1 for i in range(num_bin_in_frag)]
        hist, edges = np.histogram(pos,bins=num_bin_in_frag,range=(frag_start+1, frag_start+bin_size*num_bin_in_frag),density=False)

        cov = []
        for i in range(len(hist)):
            idx_start_include = sum(hist[:i])
            idx_end_notinclude = sum(hist[:i+1])
            cov.append(sum(cigar_match[idx_start_include: idx_end_notinclude]))

        for i in range(len(hist)-1):
            idx=sum(hist[:i+1])-1

            if pos[idx] +cigar_match[idx]>=edges[i+1]:
                cov[i+1]+=pos[idx] +cigar_match[idx]-edges[i+1]
                cov[i]-=edges[i+1]-pos[idx]

        mean_cov = np.array(cov)/bin_size
        if strand ==1:
            mean_cov_bin.append(mean_cov)
            slope, intercept, r_value, p_value, std_err = st.linregress(x, mean_cov)
        elif strand==0:
            mean_cov_bin.append(mean_cov[::-1])
            slope, intercept, r_value, p_value, std_err = st.linregress(x, mean_cov[::-1])

        # slopes.append(slope)

    return mean_cov_bin

def plot_mean_frag_bin(mean_cov_bin, frag_bin):
    xlab ='frag_bin'
    ylab='mean_cov_bin'
    save_title='mean_cov_frag_bin_'+frag_bin
    title='mean_cov_frag_bin'
    ## just plot 100 frag
    for i in range(100):
        plt.plot(mean_cov_bin[i], label=str(i))

    plt.ylim((-5, 50))
    plt.title(title)
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    plt.grid(linestyle='dotted')
    plt.savefig(save_title)


def main():
    args = parser.parse_args()
    # infile_bc_df = args.frag_bc_df
    num_bin_in_frag = args.num_bin_in_frag
    read_len = args.read_len
    dirname = args.dirname
    library_type = args.library_type
    minreads = args.minreads
    # bin_frag_by = args.bin_frag_by
    # frag_len_min = args.frag_len_min
    # frag_len_max = args.frag_len_max
    bc_df = preprocess(dirname, read_len, num_bin_in_frag, minreads)

    ### bin fragment length, expect different bias between short and long fragments
    ## TODO: robust, put frag_len_min, frag_len_max, write_binary to params
    _reads_covered_fragment_percent, _frag_mean_depth =[], []
    _top_nk_mean, _tail_nk_mean = [],[]
    
    if library_type=='gdna':
        frag_len_max = 10000
        frag_len_min = 1000
    elif library_type=='mrna':
        frag_len_max = 3000
        frag_len_min = 200
    write_binary = True
    bin_frag_by = 100
    bin_size =100
    # n = 42
    n = int(1+(frag_len_max-frag_len_min)/bin_frag_by)
    ## top or tail nk base pair seq
    nk=1000
    # nk_tag = str(nk//1000)+'k'
    
    for i in range(n):
        try:
            ## 10bins per frag
            frag_bin= frag_len_min+bin_frag_by*i
            _df = bc_df[bc_df['Frag_Length'].between(frag_len_min+bin_frag_by*i,frag_len_min+bin_frag_by*(i+1))]
            mean_cov_bin_10bins, slopes, top_nk_mean, tail_nk_mean = frag_mean_cov_10bins(_df, read_len, num_bin_in_frag, minreads, nk)
            _reads_covered_fragment_percent.append(_df['reads_covered_fragment_percent'])
            _frag_mean_depth.append((_df['frag_mean_depth']))
            # _top_nk_mean.append(top_nk_mean)
            # _tail_nk_mean.append(tail_nk_mean)

            if write_binary:
                with open(dirname+'_'.join(['/bin_num10/mean_cov_bin', str(frag_bin)]), 'wb') as fp:
                    pickle.dump(mean_cov_bin_10bins, fp)
                # with open(dirname+'_'.join(['/bin_num10/slopes', str(frag_bin)]), 'wb') as fp:
                #     pickle.dump(slopes, fp)

            ## bin_size100
            bin_size=100
            mean_cov_binSize100 = frag_mean_cov_binSize100(_df, read_len, frag_bin, bin_size, minreads, nk)

            if write_binary:
                with open(dirname+'_'.join(['/bin_size100/mean_cov_bin', str(frag_bin)]), 'wb') as fp:
                    pickle.dump(mean_cov_binSize100, fp)
        except Exception as e: 
            print(e)

    type ='bin_num10'
    with open(dirname+'_'.join(['/{}/mean_cov_bin'.format(type),str(minreads),'done']), 'w') as f:
        f.write("done\n")
    type2 ='bin_size100'
    with open(dirname+'_'.join(['/{}/mean_cov_bin'.format(type2),str(minreads),'done']), 'w') as f:
        f.write("done\n")

    with open(dirname+'_reads_covered_fragment_percent', 'wb') as fp:
        pickle.dump(_reads_covered_fragment_percent, fp)
    # _boxplot(_reads_covered_fragment_percent, dirname, data_type='reads_covered_fragment_percent')
        
    with open(dirname+'_frag_mean_depth', 'wb') as fp:
        pickle.dump(_frag_mean_depth, fp)
    # _boxplot(_frag_mean_depth, dirname, data_type='frag_mean_depth')

    # with open(dirname+'_top_nk_mean', 'wb') as fp:
    #     pickle.dump(_top_nk_mean, fp)
    # with open(dirname+'_tail_nk_mean', 'wb') as fp:
    #     pickle.dump(_tail_nk_mean, fp)
    # print(tail_nk_mean)

    ### frag depth by gc_ratio
    # try:
    #     _frag_mean_depth_gc=[]
    #     for i in range(10):
    #         _df = bc_df[bc_df['frag_GC_ratio'].between(0.1*i,0.1*(i+1))]
    #         _frag_mean_depth_gc.append((_df['frag_mean_depth']))
    #     with open(dirname+'_frag_mean_depth_gc', 'wb') as fp:
    #         pickle.dump(_frag_mean_depth_gc, fp)
    # except Exception as e: print(e)
if __name__ == "__main__":
    main()



    
