"""
utility functions for calc_frag_len.py, used for both stLFR and stcLFR
"""
from statistics import mean, stdev, median
import matplotlib.pyplot as plt
import multiprocessing as mp
import pandas as pd
import numpy as np
import seaborn as sns
import pysam, argparse, os, sys
from collections import Counter, defaultdict
import subprocess
import pickle
from ast import literal_eval
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime
import concurrent.futures

def reads_notIn_frag(dirname, barcode_collection, min_reads, min_frag, max_frag):
    """
    filter raw frag tb by min_reads4, min_frag 5000, get pct of reads not in frag
    output: barcode_collection_filtered
    """

    def writeout(test_barcodes, name):
        barcode_summary = (test_barcodes.groupby('Barcode')['N_Reads']
                           .agg(['sum', 'count'])
                           .reset_index()
                           .rename(columns={'sum':'Reads', 'count':'Frags'})
                       )

        reads_per_bc_bins = []                 
        reads_per_bc_bins.append(barcode_summary[barcode_summary['Reads']==1].Reads.sum())
        reads_per_bc_bins.append(barcode_summary[(barcode_summary['Reads'] >= 2) & (barcode_summary['Reads'] < 5)].Reads.sum())
        reads_per_bc_bins.append(barcode_summary[(barcode_summary['Reads'] >= 5) & (barcode_summary['Reads'] < 10)].Reads.sum())
        reads_per_bc_bins.append(barcode_summary[(barcode_summary['Reads'] >= 10) & (barcode_summary['Reads'] < 25)].Reads.sum())
        reads_per_bc_bins.append(barcode_summary[(barcode_summary['Reads'] >= 25) & (barcode_summary['Reads'] < 100)].Reads.sum())
        reads_per_bc_bins.append(barcode_summary[barcode_summary['Reads'] >= 100].Reads.sum())
        reads_per_bc_bins.append(len(barcode_summary[barcode_summary['Reads'] == 1].index))
        reads_per_bc_bins.append(len(barcode_summary[(barcode_summary['Reads'] >= 2) & (barcode_summary['Reads'] < 5)].index))
        reads_per_bc_bins.append(len(barcode_summary[(barcode_summary['Reads'] >= 5) & (barcode_summary['Reads'] < 10)].index))
        reads_per_bc_bins.append(len(barcode_summary[(barcode_summary['Reads'] >= 10) & (barcode_summary['Reads'] < 25)].index))
        reads_per_bc_bins.append(len(barcode_summary[(barcode_summary['Reads'] >= 25) & (barcode_summary['Reads'] < 100)].index))
        reads_per_bc_bins.append((barcode_summary[barcode_summary['Reads'] >= 100].shape[0]))
        reads_per_bc_bins.append((barcode_summary.shape[0]))

        try:
            reads_one = reads_per_bc_bins[0]
            reads_two = reads_per_bc_bins[1]
            reads_three = reads_per_bc_bins[2]
            reads_four = reads_per_bc_bins[3]
            reads_five = reads_per_bc_bins[4]
            reads_six = reads_per_bc_bins[5]

            bcs_one = reads_per_bc_bins[6]
            bcs_two = reads_per_bc_bins[7]
            bcs_three = reads_per_bc_bins[8]
            bcs_four = reads_per_bc_bins[9]
            bcs_five = reads_per_bc_bins[10]
            bcs_six = reads_per_bc_bins[11]

            bcs_total = reads_per_bc_bins[12]
            reads_total = sum(reads_per_bc_bins[0:6])
            dirname = './'
            out_stats =   f"{dirname}/{name}.txt"
            with open(out_stats, "w") as frag_stats:
                print(f"Reads/BC (1):\t{reads_one}\n"
                    f"Reads/BC (2,5]:\t{reads_two}\n"
                    f"Reads/BC (5,10]:\t{reads_three}\n"
                    f"Reads/BC (10,25]:\t{reads_four}\n"
                    f"Reads/BC (25,100]:\t{reads_five}\n"
                    f"Reads/BC (100+):\t{reads_six}\n"
                    f"Total Reads Mapped:\t{reads_total}\n"

                    f"Total BCs (1):\t{bcs_one}\n"
                    f"Total BCs (2,5]:\t{bcs_two}\n"
                    f"Total BCs (5,10]:\t{bcs_three}\n"
                    f"Total BCs (10,25]:\t{bcs_four}\n"
                    f"Total BCs (25,100]:\t{bcs_five}\n"
                    f"Total BCs (100+):\t{bcs_six}\n"
                    f"Total Barcodes Mapped:\t{bcs_total}",
                    file = frag_stats)
        except Exception as e: print('reads not in frag failed') 
    
    ## reads not in frag table
    df = barcode_collection
    test_barcodes = df[(df.N_Reads < min_reads) | (df.Frag_Length < min_frag)]
    writeout(test_barcodes, f'{dirname}/reads_notIn_frag')
    
    ## 
    f = open(dirname+'/additional_summary_minreads'+str(min_reads)+'.txt', "a")
    f.write('frag_cnt_all={}\n'.format(barcode_collection.shape[0]))
    barcode_collection_filtered = barcode_collection[barcode_collection['N_Reads'] >= min_reads]
    f.write('frag_cnt_minreadsX={}\n'.format(barcode_collection_filtered.shape[0]))
    # print('frag_cnt_minreadsX_reads={}'.format(barcode_collection_filtered['N_Reads'].sum()))
    frag_less_minreads = barcode_collection[barcode_collection['N_Reads'] < min_reads]
    f.write('frag_less_minreadsX_cnt={}\n'.format(frag_less_minreads.shape[0]))
    f.write('frag_less_minreadsX_reads={}\n'.format(frag_less_minreads['N_Reads'].sum()))

    # enough reads, but short
    frag_shorter = barcode_collection_filtered[barcode_collection_filtered['Frag_Length'] < min_frag]
    frag_longer = barcode_collection_filtered[barcode_collection_filtered['Frag_Length'] > max_frag]
    frag_shorter_total_reads= frag_shorter['N_Reads'].sum()
    frag_longer_total_reads= frag_longer['N_Reads'].sum()
    frag_shorter_cnt = frag_shorter.shape[0]
    frag_longer_cnt = frag_longer.shape[0]
    f.write('frag_shorter_total_reads={}, frag_count={}\n'.format(frag_shorter_total_reads, frag_shorter_cnt))
    f.write('frag_longer_total_reads={}, frag_count={}\n'.format(frag_longer_total_reads, frag_longer_cnt))
    f.close()
    barcode_collection_filtered = barcode_collection_filtered[barcode_collection_filtered['Frag_Length']>=min_frag]
    
    # with open(dirname+'/frag_longer', 'wb') as fp:
    #         pickle.dump(frag_longer.drop('ReadID', axis=1), fp, protocol=4)
    
    try:
        secondary_frag_distribution(barcode_collection,min_reads,min_frag)
    except Exception as e: print('secondary_frag_distribution failed')
    if barcode_collection_filtered.shape[0]==0:
        print('error: no frag pass min_reads+min_frag filter')
    print(f'barcode_collection_filtered_shape {barcode_collection_filtered.shape[0]}')
    return barcode_collection_filtered

def frag_dup_rate(dirname, barcode_collection_filtered, min_reads):
    ## per frag pct of dup reads removed
    df = barcode_collection_filtered
    if 'dup' not in df:
        df['dup'] = df['Positions_dup'].apply(len)-df['N_Reads']
    df['N_Reads_dup'] = df['Positions_dup'].apply(len)

    df['dup_per_frag']=df['dup']/df['N_Reads_dup']
    frag_dup_mean = df.loc[:, 'dup_per_frag'].mean()
    dup_sum  = sum(df.dup)
    frag_dup_total = dup_sum/(sum(df.N_Reads)+dup_sum)

    with open(dirname+'/additional_summary_minreads'+str(min_reads)+'.txt', 'a') as f:
        f.write("frag_dup_mean:\t{}\n".format(round(frag_dup_mean),2))
        f.write("frag_dup_total:\t{}\n".format(frag_dup_total))
        f.write("dup_sum:\t{}\n".format(dup_sum)) 
    return
    
def get_chroms(bam_path):
    '''
    scrape bam file to get valid chromosomes
    '''
    chroms = []
    try:
        bamfile = pysam.AlignmentFile(bam_path, "rb")
        for i in range(0, bamfile.nreferences):
            chroms.append(bamfile.get_reference_name(i))    
    finally:
        bamfile.close()
    
    return chroms

def mismatch_byN(s1, s2, num_mismatch=1):
    ## chech if 2 BC mismatch only by N
    cnt_mismatch=0
    
    for c1, c2 in zip(s1, s2):
        if c1 != c2:
            cnt_mismatch+=1
    if cnt_mismatch==num_mismatch:
        return True
    else:
        return False

def bed_sum(bed_file_merged):
    total_len_merged = 0
    with open(bed_file_merged, 'r') as f:
        for line in f:
            info = line.strip().split('\t')
            if len(info) < 3:
                continue
            chrom = info[0]
            start = int(info[1])
            end = int(info[2])
            total_len_merged +=end-start+1
    # print(total_len_merged)
    return total_len_merged

def write_merged_bed(df_bed, bed_file_merged):
    """Write a sorted, merged BED file without depending on bedtools."""
    if df_bed.empty:
        open(bed_file_merged, "w").close()
        return 0

    bed_rows = df_bed[['Chrom', 'start', 'end']].dropna().copy()
    if bed_rows.empty:
        open(bed_file_merged, "w").close()
        return 0

    bed_rows['start'] = bed_rows['start'].astype(int)
    bed_rows['end'] = bed_rows['end'].astype(int)
    bed_rows = bed_rows.sort_values(['Chrom', 'start', 'end'])

    merged_rows = []
    current_chrom = None
    current_start = None
    current_end = None

    for row in bed_rows.itertuples(index=False):
        chrom = row.Chrom
        start = int(row.start)
        end = int(row.end)
        if current_chrom is None or chrom != current_chrom or start > current_end:
            if current_chrom is not None:
                merged_rows.append((current_chrom, current_start, current_end))
            current_chrom = chrom
            current_start = start
            current_end = end
        else:
            current_end = max(current_end, end)

    if current_chrom is not None:
        merged_rows.append((current_chrom, current_start, current_end))

    with open(bed_file_merged, "w") as out:
        for chrom, start, end in merged_rows:
            out.write(f"{chrom}\t{start}\t{end}\n")

    return sum(end - start + 1 for _, start, end in merged_rows)

def manipulate_df_gaps(test_barcodes, min_frag, max_frag, min_reads, read_len):
    '''
    get distance between reads, filter dataframe and add new columns
    '''
    def get_read_gap(read_positions, read_len):
        '''
        get distances between starting positions
        '''
        gap_vals = [0]
        for i in range(1, len(read_positions)):
            _gap = read_positions[i] - (read_positions[i-1]+read_len)
            if _gap>0:
                gap_vals.append(_gap)

        return gap_vals

    # filter by minimum reads and fragment length
    # test_barcodes = test_barcodes[test_barcodes['N_Reads'] >= min_reads]
    # test_barcodes = test_barcodes[test_barcodes['Frag_Length'].between(min_frag, max_frag)]
    test_barcodes['Frag_Gaps'] = test_barcodes['Positions'].apply(get_read_gap, read_len=read_len)
    # create columns for min, max and mean distance between reads
    test_barcodes['Min_Frag_Gap'] = test_barcodes['Frag_Gaps'].apply(min)
    test_barcodes['Max_Frag_Gap'] = test_barcodes['Frag_Gaps'].apply(max)
    test_barcodes['Mean_Frag_Gap'] = round(test_barcodes['Frag_Gaps'].apply(mean),2)
    test_barcodes['gap_length_total'] = test_barcodes['Frag_Gaps'].apply(sum)
    test_barcodes['reads_covered_fragment_percent'] = 1-test_barcodes['gap_length_total']/test_barcodes['Frag_Length']
    test_barcodes['frag_mean_depth'] = round(test_barcodes['N_Reads']*read_len/test_barcodes['Frag_Length'],2)

    # inputSeqFile = open(ref, "r")
    # SeqDict = SeqIO.to_dict(SeqIO.parse(inputSeqFile, "fasta"))
    # inputSeqFile.close()
    gc_ratio_lst=[]

    try:
        if test_barcodes['Chrom'][0].startswith('chr'):
            ref_dict='/research/rv-02/home/eanderson/CGI_WGS_Pipeline/Data_and_Tools/data/hg38/SeqDict_hg38'
        else:
            ref_dict='/research/rv-02/home/eanderson/CGI_WGS_Pipeline/Data_and_Tools/data/hg38/ERCC/SeqDict_hg38cdnaERCC'
    
        with open (ref_dict, 'rb') as fp:
            _SeqDict = pickle.load(fp)
        for k,v in test_barcodes[['Chrom','Min_Pos', 'Max_Pos']].iterrows():
            chrom, start_idx, end_idx=v[0], v[1]-1, v[2]-1
            gc_ratio_lst.append(frag_gc_ratio(chrom, start_idx, end_idx, _SeqDict))
        test_barcodes['frag_GC_ratio']=gc_ratio_lst
    except Exception as e: print('frag_GC_ratio failed')   
    # print(test_barcodes.columns)
    return test_barcodes


def write_out_tsv_and_summary1(test_barcodes, dirname, write_out_tsv, min_reads, plot_cutoff):
    # Write a small preview artifact for Snakemake/downstream checks without
    # materializing the full fragment dataframe.
    columns = ['Chrom', 'Positions', 'Frag_Length', 'Cigar_match', 'Barcode', 'N_Reads']
    df = test_barcodes.reindex(columns=columns).head(100)
    df.to_csv(f"{dirname}/frag_and_bc_dataframe.tsv", sep='\t', index=False)

    unique_barcode_count = test_barcodes['Barcode'].nunique() if 'Barcode' in test_barcodes else 0
    fragment_count = test_barcodes.shape[0]
    total_frag_read_count = test_barcodes['N_Reads'].sum() if 'N_Reads' in test_barcodes else 0
    frags_per_bc = fragment_count/unique_barcode_count if unique_barcode_count else 0
    avg_frag_len = test_barcodes['Frag_Length'].sum()/fragment_count if fragment_count else 0
    med_frag_len = test_barcodes['Frag_Length'].median() if fragment_count else 0
    avg_frag_read_count = total_frag_read_count/fragment_count if fragment_count else 0
    avg_bc_read_count = total_frag_read_count/unique_barcode_count if unique_barcode_count else 0
    total_len_merged = 0
    total_len_merged_pct = 0

    try:
        ## create frag.bed file
        df_bed = test_barcodes.copy()
        df_bed['end'] = df_bed['Min_Pos']+df_bed['Frag_Length']
        df_bed['start'] = df_bed['Min_Pos']
        df_bed = df_bed[['Chrom', 'start', 'end', 'Barcode', 'N_Reads']]
        bed_file = f"{dirname}/frag_minreads"+str(min_reads)+".bed"
        df_bed.to_csv(bed_file, index=False, header=False, sep='\t')
        ## frag.merged.bed
        bed_file_merged =   f"{dirname}/frag_minreads"+str(min_reads)+"_merged.bed"
        total_len_merged = write_merged_bed(df_bed, bed_file_merged)
        GCA_000001405_15_GRCh38_len= 2934876545
        total_len_merged_pct = total_len_merged/GCA_000001405_15_GRCh38_len
    except Exception as e:
        bed_file_merged = f"{dirname}/frag_minreads"+str(min_reads)+"_merged.bed"
        open(bed_file_merged, "w").close()
        print(f'step1 failed: {e}', file=sys.stderr)
    try:
        out_stats =  f"{dirname}/frag_summary_minreads{min_reads}.txt"
        with open(out_stats, "w") as frag_stats:
            print(f"Unique Barcode Count:\t{unique_barcode_count}\n"
                f"Avg Frag Per BC:\t{frags_per_bc}\n"
                f"Avg BC Read Count:\t{avg_bc_read_count}\n"
                f"Fragment Count:\t{fragment_count}\n"
                f"Avg Frag Length:\t{avg_frag_len}\n"
                f"Median Frag Length:\t{med_frag_len}\n"
                f"Avg Frag Read Count:\t{avg_frag_read_count}\n"
                f"Total Reads in Frag:\t{total_frag_read_count}\n"
                f"Fragments_length_total_pct_genome:\t{total_len_merged_pct}\n"
                f"Fragments_length_total:\t{total_len_merged}\n"
                f"Ave_read_per_kbp:\t{1000*avg_frag_read_count/med_frag_len if med_frag_len else 0}",
                file = frag_stats)

        # mapped_bc_notin_fragment = bcs_total-unique_barcode_count 
        # mapped_bc_notin_fragment_pct = mapped_bc_notin_fragment/bcs_total
        # mapped_reads_notin_fragment = reads_total-total_frag_read_count 
        # mapped_reads_notin_fragment_pct = mapped_reads_notin_fragment/reads_total        
        # with open(dirname+'/additional_summary'+str(min_reads)+'.txt', 'a') as f:
        #     print(f"mapped_reads_notin_fragment:\t{mapped_reads_notin_fragment}\n"
        #         f"mapped_reads_notin_fragment_pct:\t{mapped_reads_notin_fragment_pct}",
        #         file = f)

    except Exception as e: print('write to frag_summary_minreads.txt failed')

    # create plots of frag length distribution and n reads
    frag_lengths = [i for i in test_barcodes['Frag_Length'] if i<=plot_cutoff]
    if frag_lengths:
        sns.distplot(frag_lengths)
    plt.savefig( f"{dirname}/frag_length_distribution.pdf")
    plt.clf()
    if not test_barcodes.empty:
        sns.distplot(test_barcodes['N_Reads'], kde=False)
    plt.savefig(f"{dirname}/n_read_distribution.pdf")
    plt.clf()

    try:
        per_frag_coverage(test_barcodes)
    except Exception as e: print('per_frag_coverage failed') 

    ## frag with sameBC that are 1MB apart
    try:
        apart_dist = 1000000
        frag_sameBC_apart(apart_dist, test_barcodes)
    except Exception as e: print('frag_sameBC_apart failed') 

def write_out_tsv_and_summary2(test_barcodes, dirname, reads_per_bc_bins, write_out_tsv, min_reads):
    try:

        reads_one = reads_per_bc_bins[0]
        reads_two = reads_per_bc_bins[1]
        reads_three = reads_per_bc_bins[2]
        reads_four = reads_per_bc_bins[3]
        reads_five_nine = reads_per_bc_bins[4]
        reads_ten_twenties = reads_per_bc_bins[5]
        reads_twenties_fifty = reads_per_bc_bins[6]
        reads_fifty_hundo = reads_per_bc_bins[7]
        reads_hundo_plus = reads_per_bc_bins[8]
        bcs_one = reads_per_bc_bins[9]
        bcs_two = reads_per_bc_bins[10]
        bcs_three = reads_per_bc_bins[11]
        bcs_four = reads_per_bc_bins[12]
        bcs_five_nine = reads_per_bc_bins[13]
        bcs_ten_twenties = reads_per_bc_bins[14]
        bcs_twenties_fifty = reads_per_bc_bins[15]
        bcs_fifty_hundo = reads_per_bc_bins[16]
        bcs_hundo_plus = reads_per_bc_bins[17]
        bcs_total = reads_per_bc_bins[18]
        reads_total = sum(reads_per_bc_bins[0:9])

        out_stats =  f"{dirname}/reads_per_bc_bins.txt"
        with open(out_stats, "w") as frag_stats:
            print(f"Reads/BC (1):\t{reads_one}\n"
                f"Reads/BC (2):\t{reads_two}\n"
                f"Reads/BC (3):\t{reads_three}\n"
                f"Reads/BC (4):\t{reads_four}\n"
                f"Reads/BC (5, 10]:\t{reads_five_nine}\n"
                f"Reads/BC (10, 25]:\t{reads_ten_twenties}\n"
                f"Reads/BC (25, 50]:\t{reads_twenties_fifty}\n"
                f"Reads/BC (50, 100]:\t{reads_fifty_hundo}\n"
                f"Reads/BC (100+):\t{reads_hundo_plus}\n"
                f"Total Reads Mapped:\t{reads_total}\n"
                f"Total BCs (1):\t{bcs_one}\n"
                f"Total BCs (2):\t{bcs_two}\n"
                f"Total BCs (3):\t{bcs_three}\n"
                f"Total BCs (4):\t{bcs_four}\n"
                f"Total BCs (5, 10]:\t{bcs_five_nine}\n"
                f"Total BCs (10, 25]:\t{bcs_ten_twenties}\n"
                f"Total BCs (25, 50]:\t{bcs_twenties_fifty}\n"
                f"Total BCs (50, 100]:\t{bcs_fifty_hundo}\n"
                f"Total BCs (100+):\t{bcs_hundo_plus}\n"
                f"Total Barcodes Mapped:\t{bcs_total}",
                file = frag_stats)
    except Exception as e: print('write to reads_per_bc_bins.txt failed')  

def write_out_barcode_summary(test_barcodes, dirname, write_out_tsv):
    '''
    get reads per bc bins and write out full tsv of fragment data.
    the full tsv is unfiltered whatsoever. if use barcode_summary[barcode_summary['Reads'].between(2,5) 会有精度问题 5包括在2个bin里
    '''
    try:
        # Create target Directory
        os.mkdir(dirname)
        print("Directory " , dirname,  " Created ") 
    except FileExistsError:
        print("Directory " , dirname,  " already exists")

    print(f"Writing out summary files to {dirname}", file=sys.stderr)
    reads_per_bc_bins=[]
    # group dataframe by N_reads, count and sum groups, rename columns
    # then subset into bins
    barcode_summary = (test_barcodes.groupby('Barcode')['N_Reads']
                           .agg(['sum', 'count'])
                           .reset_index()
                           .rename(columns={'sum':'Reads', 'count':'Frags'})
                       )
    reads_per_bc_bins.append(barcode_summary[barcode_summary['Reads'] == 1].Reads.sum())
    reads_per_bc_bins.append(barcode_summary[barcode_summary['Reads'] == 2].Reads.sum())
    reads_per_bc_bins.append(barcode_summary[barcode_summary['Reads'] == 3].Reads.sum())
    reads_per_bc_bins.append(barcode_summary[barcode_summary['Reads'] == 4].Reads.sum())
    reads_per_bc_bins.append(barcode_summary[(barcode_summary['Reads'] >= 5) & (barcode_summary['Reads'] < 10)].Reads.sum())
    reads_per_bc_bins.append(barcode_summary[(barcode_summary['Reads'] >= 10) & (barcode_summary['Reads'] < 25)].Reads.sum())
    reads_per_bc_bins.append(barcode_summary[(barcode_summary['Reads'] >= 25) & (barcode_summary['Reads'] < 50)].Reads.sum())
    reads_per_bc_bins.append(barcode_summary[(barcode_summary['Reads'] >= 50) & (barcode_summary['Reads'] < 100)].Reads.sum())
    reads_per_bc_bins.append(barcode_summary[barcode_summary['Reads'] >= 100].Reads.sum())
    reads_per_bc_bins.append(len(barcode_summary[barcode_summary['Reads'] == 1].index))
    reads_per_bc_bins.append(len(barcode_summary[barcode_summary['Reads'] == 2].index))
    reads_per_bc_bins.append(len(barcode_summary[barcode_summary['Reads'] == 3].index))
    reads_per_bc_bins.append(len(barcode_summary[barcode_summary['Reads'] == 4].index))
    reads_per_bc_bins.append(len(barcode_summary[(barcode_summary['Reads'] >= 5) & (barcode_summary['Reads'] < 10)].index))
    reads_per_bc_bins.append(len(barcode_summary[(barcode_summary['Reads'] >= 10) & (barcode_summary['Reads'] < 25)].index))
    reads_per_bc_bins.append(len(barcode_summary[(barcode_summary['Reads'] >= 25) & (barcode_summary['Reads'] < 50)].index))
    reads_per_bc_bins.append(len(barcode_summary[(barcode_summary['Reads'] >= 50) & (barcode_summary['Reads'] < 100)].index))
    reads_per_bc_bins.append(len(barcode_summary[barcode_summary['Reads'] >= 100].index))
    reads_per_bc_bins.append(len(barcode_summary.index))
    
    # output frags_per_bc.pdf
    # frag_per_bc = sns.distplot(barcode_summary['Frags'], kde=False, rug=True)
    # plt.savefig(dirname + "/frags_per_bc.png")
    # plt.clf()
    # if write_out_tsv:
    #     barcode_summary.to_csv(dirname + "/frags_reads_per_bc.tsv", sep='\t', index=False)


    ### detailed reads distribution, mapped, unmapped, moved to reads_per_bc_distribution.py
    # try:
    #     # output reads per BC distribution table, not binned table as Reads/BC, Total BCs 
    #     df = pd.DataFrame.from_dict(Counter(barcode_summary.Reads), orient='index').reset_index()
    #     df.columns = ['reads_per_bc','bc_count']
    #     df_sorted = df.sort_values('reads_per_bc')
    #     df_sorted['reads_count'] = df_sorted.bc_count*df_sorted.reads_per_bc
    #     total_reads = sum(df_sorted.reads_count)
    #     df_sorted['reads_count_pct'] =df_sorted['reads_count']/total_reads
    #     print('total_reads_mapped={}'.format(total_reads))
    #     df_sorted.to_csv(dirname + "/reads_per_bc_distribution.tsv", sep='\t', index=False, float_format='%.7f')

    #     def reads_per_bc_distribution_unmapped(dirname):
    #         split_log = 'split_stat_read1.log'
    #         summary = pd.read_csv(split_log, header=None, index_col=None, sep='\t', skiprows=4)        
    #         summary.columns = ['index', 'Reads', 'BC']

    #         df = pd.DataFrame.from_dict(Counter(summary.Reads), orient='index').reset_index()
    #         df.columns = ['reads_per_bc','bc_count']
    #         df_sorted = df.sort_values('reads_per_bc')
    #         df_sorted['reads_count'] = df_sorted.bc_count*df_sorted.reads_per_bc
    #         total_reads = sum(df_sorted.reads_count)
    #         df_sorted['reads_count_pct'] =df_sorted['reads_count']/total_reads
    #         print('total_reads_unmapped={}'.format(total_reads))
    #         df_sorted.to_csv(dirname + "/reads_per_bc_distribution_nomapping.tsv", sep='\t', index=False, float_format='%.7f')
    #         return

    #     reads_per_bc_distribution_unmapped(dirname)
    return reads_per_bc_bins, barcode_summary
    # except Exception as e: print(e)
        
def write_out_barcode_summary_nodup(test_barcodes, dirname, write_out_tsv):
    '''
    get reads per bc bins and write out full tsv of fragment data.
    the full tsv is unfiltered whatsoever.
    '''
    try:
        # Create target Directory
        os.mkdir(dirname)
        print("Directory " , dirname,  " Created ") 
    except FileExistsError:
        print("Directory " , dirname,  " already exists")

    print(f"Writing out summary files to {dirname}", file=sys.stderr)
    reads_per_bc_bins=[]
    # group dataframe by N_reads, count and sum groups, rename columns
    # then subset into bins
    barcode_summary = (test_barcodes.groupby('Barcode')['N_Reads']
                           .agg(['sum', 'count'])
                           .reset_index()
                           .rename(columns={'sum':'Reads', 'count':'Frags'})
                       )
    reads_per_bc_bins.append(barcode_summary[barcode_summary['Reads'] == 1].Reads.sum())
    reads_per_bc_bins.append(barcode_summary[barcode_summary['Reads'] == 2].Reads.sum())
    reads_per_bc_bins.append(barcode_summary[barcode_summary['Reads'] == 3].Reads.sum())
    reads_per_bc_bins.append(barcode_summary[barcode_summary['Reads'] == 4].Reads.sum())
    reads_per_bc_bins.append(barcode_summary[(barcode_summary['Reads'] >= 5) & (barcode_summary['Reads'] < 10)].Reads.sum())
    reads_per_bc_bins.append(barcode_summary[(barcode_summary['Reads'] >= 10) & (barcode_summary['Reads'] < 25)].Reads.sum())
    reads_per_bc_bins.append(barcode_summary[(barcode_summary['Reads'] >= 25) & (barcode_summary['Reads'] < 50)].Reads.sum())
    reads_per_bc_bins.append(barcode_summary[(barcode_summary['Reads'] >= 50) & (barcode_summary['Reads'] < 100)].Reads.sum())
    reads_per_bc_bins.append(barcode_summary[barcode_summary['Reads'] >= 100].Reads.sum())
    reads_per_bc_bins.append(len(barcode_summary[barcode_summary['Reads'] == 1].index))
    reads_per_bc_bins.append(len(barcode_summary[barcode_summary['Reads'] == 2].index))
    reads_per_bc_bins.append(len(barcode_summary[barcode_summary['Reads'] == 3].index))
    reads_per_bc_bins.append(len(barcode_summary[barcode_summary['Reads'] == 4].index))
    reads_per_bc_bins.append(len(barcode_summary[(barcode_summary['Reads'] >= 5) & (barcode_summary['Reads'] < 10)].index))
    reads_per_bc_bins.append(len(barcode_summary[(barcode_summary['Reads'] >= 10) & (barcode_summary['Reads'] < 25)].index))
    reads_per_bc_bins.append(len(barcode_summary[(barcode_summary['Reads'] >= 25) & (barcode_summary['Reads'] < 50)].index))
    reads_per_bc_bins.append(len(barcode_summary[(barcode_summary['Reads'] >= 50) & (barcode_summary['Reads'] < 100)].index))
    reads_per_bc_bins.append(len(barcode_summary[barcode_summary['Reads'] >= 100].index))
    reads_per_bc_bins.append(len(barcode_summary.index))
    
    # output frags_per_bc.pdf
    frag_per_bc = sns.distplot(barcode_summary['Frags'], kde=False, rug=True)
    plt.savefig( f"{dirname}/frags_per_bc.pdf")
    plt.clf()
    if write_out_tsv:
        barcode_summary.to_csv(f"{dirname}/frags_reads_per_bc.tsv", sep='\t', index=False)


    return(reads_per_bc_bins)

def frag_gc_ratio(chrom, start_idx, end_idx, SeqDict):
    _seq=str(SeqDict[chrom][start_idx:end_idx+1].seq)
    gc=_seq.count('G')+_seq.count('C')
    gc_ratio = gc/(len(_seq)+1)
    return gc_ratio

def reads_per_bc_distribution_mapped(dirname, barcode_summary):
    # in_file = dirname + f'/Calc_Frag_Length_{split_distance}/frags_reads_per_bc.tsv'
    # barcode_summary = pd.read_csv(in_file, index_col=None, sep='\t')

    # df = pd.DataFrame.from_dict(Counter(barcode_summary.Reads), orient='index').reset_index()
    df = barcode_summary.groupby('Reads')['Barcode'].apply(list).reset_index()
    df['bc_count'] = df['Barcode'].apply(lambda x: len(x))
    df.columns = ['reads_per_bc', 'bc_name', 'bc_count']
    df_sorted = df.sort_values('reads_per_bc')
    df_sorted['reads_count'] = df_sorted.bc_count*df_sorted.reads_per_bc
    total_reads = sum(df_sorted.reads_count)
    df_sorted['reads_count_pct'] =df_sorted['reads_count']/total_reads
    with open(f'{dirname}/additional_summary.txt', 'a') as f:
        print(f"total_reads_mapped:\t{total_reads}",file = f)
    # df_sorted.to_csv(dirname + f"/reads_per_bc_distribution.tsv", sep='\t', index=False, float_format='%.7f')
    return df_sorted

def safe_literal_list(value):
    if isinstance(value, list):
        return value
    if isinstance(value, tuple):
        return list(value)
    if pd.isna(value):
        return []
    try:
        parsed = literal_eval(value)
    except (ValueError, SyntaxError):
        return []
    if isinstance(parsed, list):
        return parsed
    if isinstance(parsed, tuple):
        return list(parsed)
    return [parsed]

def parse_bc100(barcode_collection, df_sorted, min_reads, reads_range, min_frag, dirname):
    ## bc100+ reads list

    df_bc100= df_sorted[df_sorted['reads_per_bc']>=reads_range].copy()
    if df_bc100.empty:
        return
    df_bc100['bc_name']=df_bc100['bc_name'].apply(safe_literal_list)
    lst_bc100 = []
    for i in df_bc100['bc_name']:
        lst_bc100.extend(i)
    if not lst_bc100:
        return
    print('Total BCs (100+): {}'.format(len(lst_bc100)))

    ## frag with bc100+
    frag_bc100= barcode_collection[barcode_collection['Barcode'].isin(lst_bc100)]
    frag_bc100.to_csv(f'{dirname}/bc100_frag_minreads{min_reads}_reads.csv', index=False)
    len_frag_all= frag_bc100.shape[0]
    print('frag_bc100_cnt={}'.format(len_frag_all))
    frag_bc100_minreadsX_pass = frag_bc100[frag_bc100.N_Reads>=min_reads]
    frag_bc100_minreadsX = frag_bc100[frag_bc100.N_Reads<min_reads]
    len_frag_minreadsX=frag_bc100_minreadsX.shape[0]
    print(f'failed_frag_minreads{min_reads}_cnt={len_frag_minreadsX}')
    _sum=frag_bc100_minreadsX['N_Reads'].sum()
    print(f'failed_frag_minreads{min_reads}_reads={_sum}')
    # frag_bc100_minreadsX.to_csv(dirname+f'failed_frag_minreads{min_reads}_reads.csv', index=False)

    frag_bc100_minlenXk = frag_bc100_minreadsX_pass[frag_bc100_minreadsX_pass.Frag_Length<min_frag]
    len_frag_minreadsX_minlenXk = frag_bc100_minlenXk.shape[0]
    print(f'failed_frag_minFragLen{min_frag}_cnt={len_frag_minreadsX_minlenXk}')
    _sum=frag_bc100_minlenXk['N_Reads'].sum()
    print(f'failed_frag_minFragLen{min_frag}_reads={_sum}')
    # frag_bc100_minlenXk.to_csv(dirname+f'failed_frag_minFragLen{min_frag}.csv', index=False)

def frag_sameBC_apart(apart_dist, test_barcodes):
    '''
    1/ percent filtered frag with same BC that are 5mb apart,
    '''
    dirname = './'
    split_distance = 50000
    if not test_barcodes:
        infile=dirname + f'/Calc_Frag_Length_{split_distance}/frag_and_bc_dataframe.tsv'
        barcode_collection2 = pd.read_csv(infile, index_col=None, sep='\t')
    else:
        barcode_collection2 = test_barcodes
    # with open (infile, 'rb') as fp:
    #     barcode_collection2 = pickle.load(fp)

    # barcode_collection2 = barcode_collection2[barcode_collection2['N_Reads'] >= min_reads]
    # barcode_collection2 = barcode_collection2[barcode_collection2['Frag_Length']>=min_frag]
    sameBC_gap=[]
    sameBC = defaultdict(list)
    sameBC_len_Nreads = defaultdict(list)
    sameBC_frag_Xmb = defaultdict(list)    
    sameBC_frag_Xmb2 = defaultdict(list)
    barcode_collection2.Positions = barcode_collection2['Positions'].apply(literal_eval)
    for i in range(barcode_collection2.shape[0]):
        bc = barcode_collection2['Barcode'][i]
        Min_Pos = int(barcode_collection2['Positions'][i][0])
        Max_Pos = int(barcode_collection2['Positions'][i][-1])
        n_reads, frag_len =barcode_collection2['N_Reads'][i], barcode_collection2['Frag_Length'][i]
        sameBC[bc].append((Min_Pos, Max_Pos))
        sameBC_len_Nreads[bc].append((n_reads, frag_len))

    for k,v in sameBC.items():
        for i in range(1,len(v)):
            gap = v[i][0]-v[i-1][1]
            # start of chr2 < end of chr1, thus 0<
            if 0<gap: 
                sameBC_gap.append(gap)
                if gap>= apart_dist:
                    sameBC_frag_Xmb[bc].append((v[i-1],v[i]))
                
    n_mb = apart_dist//1000000
    print(f'count of sameBC_frag_{n_mb}mb: {len(sameBC_frag_Xmb)}')
    print(f'percent of sameBC_frag_{n_mb}mb/all frag: {100*len(sameBC_frag_Xmb)/len(sameBC)}%')
    with open(dirname + f'/Calc_Frag_Length_{split_distance}/sameBC_frag_{n_mb}mb', 'wb') as fp:
        pickle.dump(sameBC_frag_Xmb, fp)
    
    sns.set(rc={'figure.figsize':(15,8)})
    sns.set_style("whitegrid")
    _fig = sns.distplot(sameBC_gap,bins=50,kde=True)
    plt.title(f"sameBC_gap_length_distribution cnt={len(sameBC_gap)}")
    # plt.savefig("sameBC_gap_length_distribution.pdf")
    plt.savefig( dirname + f'/Calc_Frag_Length_{split_distance}/sameBC_gap_length_distribution.downsample.pdf')
    plt.clf()

    
    return

def secondary_frag_distribution(barcode_collection2,min_reads,min_frag):
    '''
        2/ distribution of secondary frag (AAA: ([])) with N_reads, frag_len
    '''
    # if library_type.lower()=='stlfr':
    #     min_reads = 4
    #     min_frag = 5000 
    # elif library_type.lower()=='clfr':
    #     min_reads = 50
    #     min_frag = 500 

    sameBC_len_Nreads = defaultdict(list)
    dirname = './'
    split_distance = 50000
    # infile=dirname + f'/Calc_Frag_Length_{split_distance}/barcode_collection'
    # with open (infile, 'rb') as fp:
    #     barcode_collection2 = pickle.load(fp)

    for i in range(barcode_collection2.shape[0]):
        bc = barcode_collection2['Barcode'][i]
        n_reads, frag_len =barcode_collection2['N_Reads'][i], barcode_collection2['Frag_Length'][i]
        sameBC_len_Nreads[bc].append((n_reads, frag_len))

    distribution_len, distribution_Nreads, secondary_frag_cnt = [], [], []
    for k,v in sameBC_len_Nreads.items():
        if len(v)>1:
            
            _len = [i[1] for i in v]
            _Nreads = [i[0] for i in v]

            if max(_len)>=min_frag and max(_Nreads)>=min_reads:
                # print(v)
                distribution_len.extend([i for i in _len if i<min_frag])
                sec_Nreads=[i for i in _Nreads if i<min_reads]
                distribution_Nreads.extend(sec_Nreads)
                secondary_frag_cnt.append(len(sec_Nreads))

    plot_distribution_len = sns.distplot(distribution_len)
    plt.savefig(dirname + f"Calc_Frag_Length_{split_distance}/secondary_frag_distribution_len.pdf")
    plt.clf()

    plot_distribution_Nreads = sns.distplot(distribution_Nreads)
    plt.savefig(dirname + f"Calc_Frag_Length_{split_distance}/secondary_frag_distribution_Nreads.pdf")
    plt.clf()
            
    plot_secondary_frag_cnt = sns.distplot(secondary_frag_cnt)
    plt.savefig(dirname + f"Calc_Frag_Length_{split_distance}/secondary_frag_distribution_cnt.pdf")
    plt.clf()

    print(f'median secondary_frag_distribution_len {median(distribution_len)}')
    print(f'median secondary_frag_distribution_Nreads {median(distribution_Nreads)}')   
    print(f'median secondary_frag_distribution_cnt {median(secondary_frag_cnt)}')

    print(f'median secondary_frag_distribution_len {sum(distribution_len)}')
    print(f'median secondary_frag_distribution_Nreads {sum(distribution_Nreads)}')   
    print(f'median secondary_frag_distribution_cnt {sum(secondary_frag_cnt)}')

    return

def per_frag_coverage(module, min_reads):
    dirname = './'
    split_distance = 50000
    if module=='dataframe':
        subprocess.call(f"head -10000 Calc_Frag_Length_{split_distance}/frag_and_bc_dataframe.tsv > Calc_Frag_Length_{split_distance}/head.frag_and_bc_dataframe.tsv", shell=True)
        infile=dirname + f'/Calc_Frag_Length_{split_distance}/head.frag_and_bc_dataframe.tsv'
        df = pd.read_csv(infile, index_col=None, sep='\t')
        df.Cigar_match=df.Cigar_match.apply(literal_eval)
        df.Positions_dup=df.Positions_dup.apply(literal_eval)
    elif module=='pickle':
        df = test_barcodes.iloc[:10000,]

    # TODO
    df = df[df.N_Reads>=min_reads]
    print(df.shape)
    df['pct_coverage']=(df.Frag_Length -df.gap_length_total)/df.Frag_Length
    sns.set(rc={'figure.figsize':(15,8)})
    sns.set_style("whitegrid")
    _fig = sns.distplot(df.pct_coverage,bins=50,kde=True)
    plt.title("per_frag_coverage_distribution")
    # plt.savefig("per_frag_coverage_distribution.pdf")
    plt.savefig( dirname + f'/Calc_Frag_Length_{split_distance}/per_frag_coverage_distribution.downsample.pdf')
    plt.clf()
    return

def parallel_filterN100(barcode_collection):
    # with open (barcode_collection_path, 'rb') as fp:
    #     barcode_collection = pickle.load(fp)
    

    def write_file(barcode_chrom, chrom):
        try:
            # chrom = barcode_chrom.Chrom[0]
            start = barcode_chrom.iloc[0].Min_Pos
            end = barcode_chrom.iloc[-1].Max_Pos
            n = 20
            chunk_size = (end-start)//n

            tmp = barcode_chrom[barcode_chrom.N_Reads>=100]
            filter100 =tmp.ReadID
            read_idx_lst = np.array([])
            for i in filter100:
                read_idx_lst = np.append(read_idx_lst, i)
            read_idx_lst = set(read_idx_lst)
            print(f'read_idx_lst {len(read_idx_lst)}')

            def write_bam(read_idx_lst, start, end, chunk_size, i):
                bam_path = "Align/data.sort.removedup_rm000.bam"
                bamfile = pysam.AlignmentFile(bam_path, "rb")
                outfile = pysam.AlignmentFile(f"Calc_Frag_Length_50000/chunk_{chrom}.{i}.N100.bam", "wb", template=bamfile)
                s = start+i*chunk_size
                e = s+chunk_size
                for read in bamfile.fetch(chrom, s, e):
                    if read.query_name in read_idx_lst:
                        outfile.write(read)
                outfile.close()
                # bamfile.close()
                return
            
            with concurrent.futures.ThreadPoolExecutor(max_workers=20) as executor:
                futures = [executor.submit(write_bam, read_idx_lst, start, end, chunk_size, i) for i in range(n)]
                for future in concurrent.futures.as_completed(futures):
                    try:
                        res =future.result()
                    except Exception as e:
                        print(f"Error occurred: {e}")
            
            # print(f'{chrom} done')
        except Exception as e: print(e)
        return
            
    chrom_lst = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"]

    for chrom in chrom_lst:
        
        barcode_chrom = barcode_collection[barcode_collection.Chrom==chrom]
        print(f'{chrom},{barcode_chrom.shape}')
        write_file(barcode_chrom, chrom)
    return 

def parallel_filterN100_2(barcode_collection):
    '''
    chunk barcode_collection, filter Nread>100, output filterN100.bam
    '''

    def write_file(barcode_chrom, chrom):
        try:

            n = 20
            nrow = barcode_chrom.shape[0]
            chunk_size = nrow//n
            # print(chunk_size)
            def write_bam(barcode_chrom, chrom, chunk_size, i):
                df = barcode_chrom.iloc[i*chunk_size:(i+1)*chunk_size]

                start = df.iloc[0].Min_Pos
                end = df.iloc[-1].Max_Pos

                tmp = df[df.N_Reads>=100]
                filter100 =tmp.ReadID
                read_idx_lst = np.array([])
                for j in filter100:
                    read_idx_lst = np.append(read_idx_lst, j)
                read_idx_lst = set(read_idx_lst)
                print(f'read_idx_lst {len(read_idx_lst)}')

                bam_path = "Align/data.sort.removedup_rm000.bam"
                bamfile = pysam.AlignmentFile(bam_path, "rb")
                outfile = pysam.AlignmentFile(f"Calc_Frag_Length_50000/chunk_{chrom}.{i}.N100.bam", "wb", template=bamfile)
                s = start+i*chunk_size
                e = s+chunk_size
                for read in bamfile.fetch(chrom, s, e):
                    if read.query_name in read_idx_lst:
                        outfile.write(read)
                outfile.close()
                # bamfile.close()
                return
            
            with concurrent.futures.ThreadPoolExecutor(max_workers=20) as executor:
                futures = [executor.submit(write_bam, barcode_chrom, chrom, chunk_size, i) for i in range(n)]
                for future in concurrent.futures.as_completed(futures):
                    try:
                        res =future.result()
                    except Exception as e:
                        print(f"Error occurred: {e}")
            
            # print(f'{chrom} done')
        except Exception as e: print(e)
        return
            
    chrom_lst = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"]

    for chrom in chrom_lst:
        barcode_chrom = barcode_collection[barcode_collection.Chrom==chrom].reset_index(drop=True)
        print(f'{chrom},{barcode_chrom.shape}')
        write_file(barcode_chrom, chrom)
    return
