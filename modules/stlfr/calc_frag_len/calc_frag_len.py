#!/usr/bin/env python3
"""calc_frag_len.py

Usage:
  calc_frag_len.py [--threads <int>] [--splitdist <int>] 
                   [--minfrag <int>] [--minreads <int>] 
                   [--readlen <int>] [--chroms <comma,separated,chroms>] 
                   [--includedups] <bampath>

Options:
  -h --help                    Show this screen.
  -m <int> --minfrag <int>     Minimum fragment size [default: 750].
  -s <int> --splitdist <int>   Distance between reads to call a new fragment [default: 300000].
  --minreads <int>             Minimum number of reads to keep a fragment [default: 4].
  --readlen <int>              Read length of dataset [default: 100].
  -n <int> --threads <int>     Number of threads to run in parallel [default: 1].
  --chroms                     Comma seperated list of chromosomes to use [default: all].
  --includedups                Include reads marked as duplicates.
  --writeouttsvs               Write out large, detailed TSVs of various data manipulations. 
  --outdir                     Specify output directory [default ./Calc_Frag_Length]
"""

# from random import sample, seed
from statistics import mean, stdev, median
import matplotlib.pyplot as plt
import multiprocessing as mp
import pandas as pd
import seaborn as sns
import pysam, argparse, os, sys
import ray
from collections import Counter, defaultdict
import pickle
from utility import mismatch_byN, manipulate_df_gaps, write_out_barcode_summary, write_out_tsv_and_summary1, write_out_tsv_and_summary2, get_chroms
from utility import reads_notIn_frag, frag_dup_rate, reads_per_bc_distribution_mapped, parse_bc100 
import pickle
from subprocess import call 
# Also SHOULD validate the arguments passed

def main():
    '''
    run fragment calculations using ray
    '''
    args = get_arguments()
    bam_path = args.bampath
    min_frag = args.minfrag
    max_frag = 300000
    plot_cutoff = max_frag
    split_dist = args.splitdist
    min_reads = args.minreads
    read_len = args.readlen
    include_dups = args.includedups
    write_out_tsv = args.writeouttsvs
    mapping_quality = args.mapping_quality
    n_threads = args.threads
    chroms = args.chroms
    dirname = args.outdir
    ray_mem = args.ray_mem
    # print(ray_mem)
    ray_mem_lst = [int(i) for i in ray_mem.split(",")]
    reads_range=100
    print(f"Calculating Fragment Lengths for Chroms {chroms}", file=sys.stderr)
    # shut down ray, sometimes it's already initialized and causes an error
    ray.shutdown()
    # initialize ray, 500M=memory=500 * 1024 * 1024, 2.5G=memory=2500 * 1024 * 1024
    ray.init(num_cpus=n_threads, object_store_memory=ray_mem_lst[0]*ray_mem_lst[1]*1024*1024, memory=ray_mem_lst[2]*ray_mem_lst[3]*1024*1024)
    print(ray.available_resources(), file=sys.stderr)

    # process bam file
    barcode_collection, read_flag_failed, poorQaulity, bc_flag_failed, mapped_flag_failed = get_reads(bam_path, split_dist, include_dups, 
                                   chroms, read_len, mapping_quality)
    try:             
        n = 10000                 
        with open(dirname+'/barcode_collection', 'wb') as fp:
            pickle.dump(barcode_collection.iloc[n:(n+10000),].drop('ReadID', axis=1), fp, protocol=4)      
        with open(dirname+'/barcode_collection10k', 'wb') as fp:
            # skip ERCC lines
            
            pickle.dump(barcode_collection.iloc[n:(n+10000),].drop('ReadID', axis=1), fp, protocol=4)
    except Exception as e: print('write barcode_collection failed')

    with open(dirname+'/additional_summary_minreads'+str(min_reads)+'.txt', 'a') as f:
        print(f"read_flag_failed:\t{read_flag_failed}\n"
            f"poorQaulity:\t{poorQaulity}\n"
            f"mapped_flag_failed:\t{mapped_flag_failed}\n"
            f"bc_flag_failed:\t{bc_flag_failed}",
            file = f)                            
    # write out barcodes summary prior to any filtering
    raw_reads_per_bc_bins, barcode_summary = write_out_barcode_summary(barcode_collection, dirname, write_out_tsv)
    # filter to keep barcodes greater than min_frag and min_reads, process
    try:
        barcode_collection_filtered = reads_notIn_frag(dirname, barcode_collection, min_reads, min_frag, max_frag)
    except Exception as e: print('reads_notIn_frag failed')

    barcode_collection_filtered = manipulate_df_gaps(barcode_collection_filtered, min_frag, max_frag, min_reads, read_len)
    ### frag-dependent dup rate
    # try:
    #     frag_dup_rate(dirname, barcode_collection_filtered)
    # except Exception as e: print(e)
    # output fragment calculations
    # write_out_tsv_and_summary(barcode_collection, dirname, raw_reads_per_bc_bins, write_out_tsv, min_reads)
    write_out_tsv_and_summary2(barcode_collection, dirname, raw_reads_per_bc_bins, write_out_tsv, min_reads)
    write_out_tsv_and_summary1(barcode_collection_filtered, dirname, write_out_tsv, min_reads, plot_cutoff)

    try:
        script = f"cat {dirname}/frag_summary_minreads{min_reads}.txt > {dirname}/frag_and_bc_summary.txt && grep frag_shorter_total_reads {dirname}/additional_summary_minreads{min_reads}.txt -A 1 >> {dirname}/frag_and_bc_summary.txt && cat {dirname}/reads_per_bc_bins.txt >> {dirname}/frag_and_bc_summary.txt"
        call(script, shell=True)
    except Exception as e: print('write to frag_and_bc_summary.txt failed')

    try:
        df_sorted = reads_per_bc_distribution_mapped(dirname, barcode_summary)
        parse_bc100(barcode_collection, df_sorted, min_reads, reads_range, min_frag, dirname)
    except Exception as e: print('parse_bc100 failed')

    return
    

def get_arguments():
    '''
    Use arg parse to get arguments
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument("--bampath", type=str, help="Path to bamfile for fragment length calculations.")
    parser.add_argument("-m", "--minfrag", type=int, default=750,
                        help="Minimum fragment size [default: 750].")
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
    parser.add_argument("--outdir", type=str, default="Calc_Frag_Length",
                        help="Specify output directory [default ./Calc_Frag_Length]")
    parser.add_argument("--ray_mem", type=str, default="200, 1024, 100, 1024",
                        help="mem allocation for Ray.")
    parser.add_argument("--mapping_quality", type=int, default=30, help="mapping quality cutoff")                                 
    args = parser.parse_args()
    
    # scrape valid chromosomes from bam
    valid_chroms = get_chroms(args.bampath)
    
    # if chroms supplied, validate them against the bam file
    if args.chroms:
        chroms = args.chroms.split(",") # validate these real quick...
        chrom_test = all(chrom in valid_chroms for chrom in chroms)
        if not chrom_test:
            parser.error("""Chroms supplied: {}
                            Chromosomes supplied must be valid chromosomes.
                            Valid Chromosomes: {}""".format(chroms, valid_chroms))
        args.chroms = chroms
    # else use all chroms from bam file
    else:
        args.chroms = valid_chroms
    
    return args


@ray.remote(num_return_vals=5)
def p_processing_chroms(bam_path, include_dups, split_dist, 
            chrom, read_len, mapping_quality):
    '''
    remote process to collect barcode and fragment information
    '''
    barcode_coll = {} # Collection of barcodes to be turned into a df
    barcode_subs = {} # Collection of barcodes tracking sub-fragments
    read_flag_failed = 0
    mapped_flag_failed =0
    poorQaulity = 0
    bc_flag_failed =0
    n = 0
    # bc_pos_dict = defaultdict(list)

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
            if read.reference_start - last_pos >= split_dist:
                barcode_subs[bc_identifier] += 1

        # return tag for the fragment
        bc_tag = bc_identifier + "_" + str(barcode_subs[bc_identifier])
        return bc_tag

    # create alignment file object
    bamfile = pysam.AlignmentFile(bam_path, "rb")
    # print(f"Processing Chrom: {chrom}", file=sys.stderr)

    # process alignment file
    for read in bamfile.fetch(chrom):
        n+=1
        # bc = read.query_name.split("#")[1]
        bc = read.get_tag('BX')

        pos = (chrom, str(read.reference_start))
        n_count = bc.count("N")
        bc_flag, read_flag, mapped_flag, mapq_flag=True, True, True, True
        if bc == "0_0_0":
            bc_flag=False
            bc_flag_failed+=1

        if include_dups: # 2048(chimeric align) + 256(secondary align)
            if read.flag & 0x900 != 0:
                read_flag=False
                read_flag_failed+=1
        else: # + 1024 (pcr dup)
            if read.flag & 0xD00 != 0:
                read_flag=False
                read_flag_failed+=1
        if read.flag & 0x4!=0:
            mapped_flag = False
            mapped_flag_failed+=1
        if read.mapping_quality < mapping_quality:
            mapq_flag = False
            poorQaulity+=1
        
        if read.flag & 0x10==16:
            read_strand = 0
        else:
            read_strand = 1

        # check if we have a valid read
        if bc_flag and mapq_flag and read_flag and mapped_flag:
            chrom = bamfile.get_reference_name(read.reference_id)
            cigar = read.cigartuples
            match_lst = []
            for c in cigar:
                if c[0]==0:
                    match_lst.append(c[1])
            cigar_match = sum(match_lst)

            bc_tag = check_bc_tag(bc, chrom)
            # create new entry for unseen bc_tag
            if bc_tag not in barcode_coll:
                barcode_coll[bc_tag] = [bc, chrom, [read.reference_start], [read.query_name],[read_strand], [cigar_match]]
            # otherwise add start and read name
            else:
                barcode_coll[bc_tag][2].append(read.reference_start)
                barcode_coll[bc_tag][3].append(read.query_name)
                barcode_coll[bc_tag][4].append(read_strand)
                barcode_coll[bc_tag][5].append(cigar_match)
        ## % of mapped reads that do not have barcode
                

    # try:
    #     with open('additional_summary.txt', 'w') as f:
    #         f.write("read_flag_failed:\t{}\n".format(read_flag_failed))
    #         f.write("poorQaulity:\t{}\n".format(poorQaulity))
    # except:
    # print("read_flag_failed:\t{}\n".format(read_flag_failed))
    # print("poorQaulity:\t{}\n".format(poorQaulity))
    # print(len(bc_pos_dict))
    # perform preliminary manipulation of df
    print(f"f{chrom} reads count {n}\n")
    barcode_df = manipulate_df_prelim(barcode_coll, read_len)
    return barcode_df, read_flag_failed, poorQaulity, bc_flag_failed, mapped_flag_failed

def get_reads(bam_path, split_dist, include_dups, 
              chroms, read_len, mapping_quality):
    '''
    create tasks for processing bam for each chromosome in parallel.
    then aggregate results.
    '''
    barcode_collection_list, read_flag_failed_list, poorQaulity_list, bc_flag_failed_list, mapped_flag_failed_list = [], [], [], [], []
    # create jobs for ray for each chromosome
    for chrom in chroms:
        object_id = p_processing_chroms.remote(bam_path, include_dups, split_dist, 
                                               chrom, read_len, mapping_quality)
        barcode_collection_list.append(object_id[0])
        read_flag_failed_list.append(object_id[1])
        poorQaulity_list.append(object_id[2])
        bc_flag_failed_list.append(object_id[3])
        mapped_flag_failed_list.append(object_id[4])

    # get all results from ray
    _barcode_collection = ray.get(barcode_collection_list)
    _read_flag_failed = ray.get(read_flag_failed_list)
    _poorQaulity = ray.get(poorQaulity_list)
    _bc_flag_failed = ray.get(bc_flag_failed_list)
    _mapped_flag_failed = ray.get(mapped_flag_failed_list)

    # _bc_pos_dict = ray.get(bc_pos_dict_list)
    # insure that all chromosomes were processed
    assert len(_barcode_collection) == len(chroms), f"length of done_tasks is not equal to chroms {len(done_tasks)} {len(chroms)}"
    # aggregate results into one df
    barcode_collection = pd.concat(_barcode_collection, ignore_index=True)
    # print(len(barcode_collection))
    read_flag_failed = sum(_read_flag_failed)
    poorQaulity = sum(_poorQaulity)
    bc_flag_failed = sum(_bc_flag_failed)
    mapped_flag_failed = sum(_mapped_flag_failed)
   
    # with open('bc_pos_dict', 'wb') as fp:
    #     pickle.dump(bc_pos_dict, fp)
    return barcode_collection, read_flag_failed, poorQaulity, bc_flag_failed, mapped_flag_failed

def manipulate_df_prelim(barcode_collection, read_len):
    '''
    perform preliminary dataframe manipulations on barcode fragment data
    '''
    # create df from barcode_collection dict
    test_barcodes = pd.DataFrame.from_dict(barcode_collection, orient="index", columns=['Barcode', 'Chrom', 'Positions', 'ReadID','Strand', 'Cigar_match'])
    test_barcodes['Strand'] = test_barcodes['Strand'].apply(lambda x: Counter(x).most_common(1)[0][0])
    test_barcodes['Cigar_match'] = test_barcodes['Cigar_match'].apply(list)
    # output dup position, check if dup cluster
    test_barcodes['Positions_dup'] = test_barcodes['Positions'].apply(sorted)
    test_barcodes['Positions'] = test_barcodes['Positions_dup']
    # test_barcodes['Positions'] = test_barcodes['Positions'].apply(set).apply(list).apply(sorted)
    # dedup, get number of reads, min and max positions
    test_barcodes['N_Reads'] = test_barcodes['Positions'].apply(len)
    test_barcodes['Min_Pos'] = test_barcodes['Positions'].apply(min)
    test_barcodes['Max_Pos'] = test_barcodes['Positions'].apply(max)
    # get total length of fragment
    test_barcodes['Frag_Length'] = test_barcodes['Max_Pos'] - test_barcodes['Min_Pos'] + read_len

    return test_barcodes

if __name__ == "__main__":
    main()

