#!/usr/bin/env python3

'''
calculate frag 
output .bed, each line as 1 read mapped region, not a frag
'''

from statistics import mean, stdev, median
import matplotlib.pyplot as plt
import multiprocessing as mp
import pandas as pd
import numpy as np
import seaborn as sns
import pysam, argparse, os, sys
import ray
from collections import Counter, defaultdict
from utility3 import mismatch_byN, manipulate_df_gaps, write_out_barcode_summary, write_out_tsv_and_summary1, write_out_tsv_and_summary2, get_chroms
from utility3 import reads_notIn_frag, frag_dup_rate, reads_per_bc_distribution_mapped, parse_bc100
# from utility3 import parallel_filterN100_2
import pickle
from subprocess import call 
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime
import concurrent.futures


def main():
    '''
    run fragment calculations using ray
    '''
    args = get_arguments()
    bam_path = args.bampath
    min_frag = args.minfrag
    max_frag = 15000
    plot_cutoff = 7000
    split_dist = args.splitdist
    min_reads = args.minreads
    read_len = args.readlen
    include_dups = args.includedups
    write_out_tsv = args.writeouttsvs
    n_threads = args.threads
    chroms = args.chroms
    dirname = args.outdir
    n_tolerance = args.n_tolerance
    cbc_len = args.cbc_len
    sequence_type = args.sequence_type
    BC_type = args.BC_type
    mapping_quality = args.mapping_quality
    # nearest_bc = args.nearest_bc
    nearest_bc=5
    reads_range=100
    # print(f"Calculating Fragment Lengths for Chroms {chroms}", file=sys.stderr)
    # shut down ray, sometimes it's already initialized and causes an error
    ray.shutdown()
    # initialize ray
    ray.init(num_cpus=n_threads, object_store_memory=60*1024*1024*1024, memory=30*1024*1024*1024)
    print(ray.available_resources(), file=sys.stderr)

    # process bam file
    barcode_collection = get_reads(bam_path, split_dist, chroms, read_len, n_tolerance)
    # print(barcode_collection.columns)
    
    # try:                      
    #     parallel_filterN100_2(barcode_collection)
    # except Exception as e: print(e)

    try:
        for chrom in [chroms]:
            # with open(f'Align/barcode_collection_{chrom}', 'wb') as fp:
            #     pickle.dump(barcode_collection, fp)
            # barcode_collection.to_csv(f'Align/{chrom}.barcode_collection.csv', index=False)
            # print(median(barcode_collection['N_Reads']))
            df = barcode_collection[barcode_collection['N_Reads']>=min_reads]
            df['N_line'] = df['line'].apply(len)
            print(f'filter N_Reads {df.shape}')
            df = df[df.Frag_Length>=min_frag]
            print(f'filter Frag_Length {df.shape}')
            # df['line'].to_csv(f'Align/{chrom}.line.csv', index=False, header=False)

            # df[['Barcode','Positions','N_Reads','Min_Pos','Max_Pos','Frag_Length', 'N_line']].to_csv(f'Align/{chrom}.filter100.csv', index=False)

            with open(f'Align/{chrom}.chr2bx', 'w') as f:
                # Barcode,Chrom,Positions,line,N_Reads,Min_Pos,Max_Pos,Frag_Length
                # for i in range(df.shape[0]):
                #     for j in df['line'][i]:
                #         f.write(f"{df[Barcode][i]}+'\t'+{j}+\n")
                for k,v in df.iterrows():
                    bx = v[0]
                    new_line = v[3]
                    for line in new_line:
                        f.write(f"{bx}"+'\t'+f"{line}"+"\n")
        
    except Exception as e: print(e)

    

    # with open(dirname+'/additional_summary_minreads'+str(min_reads)+'.txt', 'a') as f:
    #     print(f"read_flag_failed:\t{read_flag_failed}\n"
    #         f"poorQaulity:\t{poorQaulity}\n"
    #         f"mapped_flag_failed:\t{mapped_flag_failed}\n"
    #         f"bc_flag_failed:\t{bc_flag_failed}",
    #         file = f)
    # # write out barcodes summary prior to any filtering
    # raw_reads_per_bc_bins, barcode_summary = write_out_barcode_summary(barcode_collection, dirname, write_out_tsv)
    # # with open(dirname+'/raw_reads_per_bc_bins', 'wb') as fp:
    # #     pickle.dump(raw_reads_per_bc_bins, fp)
    # # filter to keep barcodes greater than min_frag and min_reads, process

    # try:
    #     barcode_collection_filtered = reads_notIn_frag(dirname, barcode_collection, min_reads, min_frag, max_frag)
    # except Exception as e: print(e)

    # barcode_collection_filtered = manipulate_df_gaps(barcode_collection_filtered, min_frag, max_frag, min_reads, read_len)

    # ### frag-dependent dup rate
    # try:
    #     frag_dup_rate(dirname, barcode_collection_filtered, min_reads)
    # except Exception as e: print(e)

    # # output fragment calculations
    # write_out_tsv_and_summary2(barcode_collection_filtered, dirname, raw_reads_per_bc_bins, write_out_tsv, min_reads)
    # write_out_tsv_and_summary1(barcode_collection_filtered, dirname, write_out_tsv, min_reads, plot_cutoff)

    # try:
    #     script = "cat {dirname}/frag_summary_minreads{min_reads}.txt > {dirname}/frag_and_bc_summary.txt && head -2 {dirname}/additional_summary_minreads{min_reads}.txt >> {dirname}/frag_and_bc_summary.txt && cat {dirname}/raw_reads_per_bc_bins.txt >> {dirname}/frag_and_bc_summary.txt".format(dirname=dirname, min_reads=str(min_reads))    
    #     call(script, shell=True)
    #     # mean_dup_per_frag=barcode_collection_filtered.loc[:, 'dup_per_frag'].mean()
    #     # mean_dup_all_frag = sum(barcode_collection_filtered.dup)/sum(barcode_collection_filtered.N_Reads)
    #     # print("mean_dup_per_frag={}, mean_dup_all_frag={}".format(mean_dup_per_frag, mean_dup_all_frag))
    # except Exception as e: print(e)

    # try:
    #     df_sorted = reads_per_bc_distribution_mapped(dirname, barcode_summary)
    #     parse_bc100(barcode_collection, df_sorted, min_reads, reads_range, min_frag, dirname)
    # except Exception as e: print(e)
    
    # return

    
def get_arguments():
    '''
    Use arg parse to get arguments
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument("--bampath", type=str, help="Path to bamfile for fragment length calculations.")
    parser.add_argument("-m", "--minfrag", type=int, default=750,
                        help="Minimum fragment size [default: 750].")
    # parser.add_argument("-x", "--maxfrag", type=int, default=10000,
                        # help="Max fragment size [default: 10000] to eliminate randomBC inflation.")
    parser.add_argument("-s", "--splitdist", type=int, default=300000,
                        help="Distance between reads to call a new fragment [default: 300000].")
    parser.add_argument("--minreads", type=int, default=4,
                        help="Minimum number of reads to keep a fragment [default: 4].")
    parser.add_argument("--readlen", type=int, default=100,
                        help="Read length of dataset [default: 100].")
    parser.add_argument("-n", "--threads", type=int, default=1,
                        help="Number of threads to run in parallel [default: 1].")
    parser.add_argument("--chroms", type=str,
                        help="Comma seperated list of chromosomes to use [default: all]")
    parser.add_argument("--includedups", action="store_true",
                        help="Include reads marked as duplicates.")
    parser.add_argument("--writeouttsvs", action="store_true",
                        help="Write out large, detailed TSVs of various data manipulations.")
    parser.add_argument("--mapping_quality", type=int, default=30, help="mapping quality cutoff")                    
    parser.add_argument("--outdir", type=str, default="Calc_Frag_Length",
                        help="Specify output directory [default ./Calc_Frag_Length]")
    parser.add_argument("--n_tolerance", type=int, default=2,
                        help="Number of Ns to tolerate")
    parser.add_argument("--cbc_len", type=int, default=15,
                        help="select BC length to test, 15 to 18")                    
    parser.add_argument("--BC_type", type=str, help="random_bc or combined")                    

    parser.add_argument("--sequence_type", type=str, help="pe or se")                    
    args = parser.parse_args()
    
    # # scrape valid chromosomes from bam
    # valid_chroms = get_chroms(args.bampath)
    
    # # if chroms supplied, validate them against the bam file
    # if args.chroms:
    #     chroms = args.chroms.split(",") # validate these real quick...
    #     chrom_test = all(chrom in valid_chroms for chrom in chroms)
    #     if not chrom_test:
    #         parser.error("""Chroms supplied: {}
    #                         Chromosomes supplied must be valid chromosomes.
    #                         Valid Chromosomes: {}""".format(chroms, valid_chroms))
    #     args.chroms = chroms
    # # else use all chroms from bam file
    # else:
    #     args.chroms = valid_chroms
    
    return args


@ray.remote(num_return_vals=2)
def p_processing_chroms(bam_path, split_dist, chrom, read_len, n_tolerance):
    '''
    remote process to collect barcode and fragment information
    '''
    n = 0
    barcode_coll = {} # Collection of barcodes to be turned into a df
    barcode_subs = {} # Collection of barcodes tracking sub-fragments
    read_flag_failed = 0
    mapped_flag_failed =0
    poorQaulity = 0
    bc_flag_failed =0
    read_idx_lst = np.array([])

    # bc_pos_dict = defaultdict(list)
    # barcode_orderedlist = []

    def check_bc_tag(bc, chrom):
        '''
        check to see if barcode_chrom combination exists.
        if it does, create new sub-fragment
        '''
        bc_identifier = bc + "_" + chrom
        # add barcode with subfragment 0
        if bc_identifier not in barcode_subs:
            barcode_subs[bc_identifier] = 0

        # check if barcode subs last position is within split distance
        else:
            last_pos = barcode_coll[bc_identifier + "_" + str(barcode_subs[bc_identifier])][2][-1]
            # add new sub
            if reference_start - last_pos >= split_dist:
                barcode_subs[bc_identifier] += 1

        # return tag for the fragment
        bc_tag = bc_identifier + "_" + str(barcode_subs[bc_identifier])
        return bc_tag
    
    # create alignment file object
    # bamfile = pysam.AlignmentFile(bam_path, "rb")
    # print(f"Processing Chrom: {chrom}", file=sys.stderr)
    last_id = 'x'
    last_bc_tag = 'x'
    strand_dict = {}
    convert_strand = {'+':'-', '-':'+'}
    with open(f'Align/{chrom}_bed', 'r') as f:
        # chr21	5348879	5348889	E200008241_L01_C034R007_2463745#AACTTATGAACTTGA/1	0	+
        # process alignment file, BC as 15bp in readID, stLFR+transposonBC: #0_0_0#ATCGATCGAT#/1
        for line in f:
            n+=1
            info = line.strip().split('\t')
            bc = info[3].split('#')[1].split('/')[0]
            n_count = bc.count("N")
            readid= info[3]
            reference_start = int(info[1])
            if n_count<=n_tolerance:
                if readid!=last_id:
                    bc_tag = check_bc_tag(bc, chrom)
                    last_bc_tag = bc_tag
                    last_id = readid
                    if bc_tag not in barcode_coll:
                        tmp_strand = info[-1]
                        r1r2 = readid[-1]
                        if r1r2=='2':
                            strand =  convert_strand[tmp_strand] 
                        else:
                            strand =  tmp_strand
                        strand_dict[bc_tag] = strand

                        new_line = '\t'.join(info[1:-1]+[strand])
                        barcode_coll[bc_tag] = [bc, chrom, [reference_start], [new_line]]
                    else:
                        new_line = '\t'.join(info[1:-1]+[strand_dict[bc_tag]])
                        barcode_coll[bc_tag][2].append(reference_start)
                        barcode_coll[bc_tag][3].append(new_line)
                else:
                    # bed from same read
                    new_line = '\t'.join(info[1:-1]+[strand_dict[last_bc_tag]])
                    barcode_coll[last_bc_tag][3].append(new_line)

    # else:
    #     print("error: check sequence_type input")

    barcode_df = manipulate_df_prelim(barcode_coll, read_len)
    # print(manipulate_df_prelim.shape)
    n = 0
    return barcode_df, n

def get_reads(bam_path, split_dist, chroms, read_len, n_tolerance):
    '''
    create tasks for processing bam for each chromosome in parallel.
    then aggregate results.
    '''
    n_list = []
    barcode_collection_list, read_flag_failed_list, poorQaulity_list, bc_flag_failed_list, mapped_flag_failed_list = [], [], [], [], []
    # create jobs for ray for each chromosome
    for chrom in [chroms]:
        object_id = p_processing_chroms.remote(bam_path, split_dist, chrom, read_len, n_tolerance)
        barcode_collection_list.append(object_id[0])
        n_list.append(object_id[1]) 
    
    # get all results from ray
    _barcode_collection = ray.get(barcode_collection_list)
    _n = ray.get(n_list)
    barcode_collection = pd.concat(_barcode_collection, ignore_index=True)
    n = sum(_n)

    # _bc_pos_dict = ray.get(bc_pos_dict_list)
    # insure that all chromosomes were processed
    # assert len(_barcode_collection) == len(chroms), f"length of done_tasks is not equal to chroms {len(done_tasks)} {len(chroms)}"
            
    # aggregate results into one df

    # print(len(barcode_collection))
    # read_flag_failed = sum(_read_flag_failed)
    # poorQaulity = sum(_poorQaulity)
    # bc_flag_failed = sum(_bc_flag_failed)
    # mapped_flag_failed = sum(_mapped_flag_failed)

    # with open('bc_pos_dict', 'wb') as fp:
    #     pickle.dump(bc_pos_dict, fp)

    return barcode_collection

def manipulate_df_prelim(barcode_collection, read_len):
    '''
    perform preliminary dataframe manipulations on barcode fragment data
    '''
    # create df from barcode_collection dict
    test_barcodes = pd.DataFrame.from_dict(barcode_collection, orient="index", columns=['Barcode','Chrom', 'Positions','line'])
    # print(test_barcodes.shape)
    test_barcodes['Positions'] = test_barcodes['Positions'].apply(sorted)
    test_barcodes['N_Reads'] = test_barcodes['Positions'].apply(len)
    test_barcodes['Min_Pos'] = test_barcodes['Positions'].apply(min)
    test_barcodes['Max_Pos'] = test_barcodes['Positions'].apply(max) + read_len
    test_barcodes['Frag_Length'] = test_barcodes['Max_Pos'] - test_barcodes['Min_Pos']
    

    return test_barcodes

# def parse_bc(bc_path):
#         bc_lst = []
#         with open(bc_path) as f:
#             for line in f:
#                 info = line.strip().split('\t')
#                 bc = info[0]
#                 bc_lst.append(bc)
#         return bc_lst

# def convert_stLFR(bc_lst):
#     bc_dict = {}
#     replacements = ['A', 'T', 'C', 'G']
#     for bc in bc_lst:
#         for i in range(len(bc)):
#             for replacement in replacements:                    
#                 new_bc = bc[:i] + replacement + bc[i+1:]
#                 bc_dict[new_bc]= bc
#     return bc_dict

if __name__ == "__main__":
    main()

