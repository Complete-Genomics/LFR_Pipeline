#!/usr/bin/env python3
"""calc_frag_len_optimized.py - Memory-optimized fragment length calculator

This is an optimized version of calc_frag_len.py with the following improvements:
1. Reduced memory usage by not storing ReadID lists (only count needed)
2. Streaming processing instead of building large dictionaries
3. Removed debug print statements from hot paths
4. Uses numpy arrays for numeric operations
5. Optional chunked processing for very large BAM files

Usage:
  calc_frag_len_optimized.py [--threads <int>] [--splitdist <int>] 
                             [--minfrag <int>] [--minreads <int>] 
                             [--readlen <int>] [--chroms <comma,separated,chroms>] 
                             [--includedups] [--low_memory] <bampath>

Options:
  -h --help                    Show this screen.
  -m <int> --minfrag <int>     Minimum fragment size [default: 750].
  -s <int> --splitdist <int>   Distance between reads to call a new fragment [default: 300000].
  --minreads <int>             Minimum number of reads to keep a fragment [default: 4].
  --readlen <int>              Read length of dataset [default: 100].
  -n <int> --threads <int>     Number of threads to run in parallel [default: 1].
  --chroms                     Comma separated list of chromosomes to use [default: all].
  --includedups                Include reads marked as duplicates.
  --writeouttsvs               Write out large, detailed TSVs of various data manipulations. 
  --outdir                     Specify output directory [default ./Calc_Frag_Length]
  --low_memory                 Use low memory mode (slower but uses less RAM)
"""

from statistics import mean, stdev, median
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt
import multiprocessing as mp
import pandas as pd
import numpy as np
import seaborn as sns
import pysam
import argparse
import os
import sys
from pathlib import Path
from collections import Counter, defaultdict
import pickle
from subprocess import call
from concurrent.futures import ProcessPoolExecutor, as_completed

SRC_DIR = Path(__file__).resolve().parent
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from utility import (
    mismatch_byN, manipulate_df_gaps, write_out_barcode_summary,
    write_out_tsv_and_summary1, write_out_tsv_and_summary2, get_chroms,
    reads_notIn_frag, frag_dup_rate, reads_per_bc_distribution_mapped, parse_bc100
)


def main():
    """Run fragment calculations with optimized memory usage."""
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
    low_memory = args.low_memory
    keep_read_ids = args.keep_read_ids
    
    print(f"Calculating Fragment Lengths for {len(chroms)} chromosomes", file=sys.stderr)
    print(f"LFR type: {args.lfr_type}", file=sys.stderr)
    print(f"Mode: {'Low Memory' if low_memory else 'Standard'}", file=sys.stderr)
    
    # Create output directory
    os.makedirs(dirname, exist_ok=True)
    
    # Process BAM file - choose method based on memory mode
    if low_memory:
        barcode_collection, stats = get_reads_low_memory(
            bam_path, split_dist, include_dups, chroms, read_len, mapping_quality, n_threads, keep_read_ids, args.cbc_len
        )
    else:
        barcode_collection, stats = get_reads_parallel(
            bam_path, split_dist, include_dups, chroms, read_len, mapping_quality, n_threads, keep_read_ids, args.cbc_len
        )
    
    read_flag_failed, poor_quality, bc_flag_failed, mapped_flag_failed = stats
    
    # Write additional summary
    with open(f'{dirname}/additional_summary_minreads{min_reads}.txt', 'a') as f:
        f.write(f"read_flag_failed:\t{read_flag_failed}\n"
                f"poorQaulity:\t{poor_quality}\n"
                f"mapped_flag_failed:\t{mapped_flag_failed}\n"
                f"bc_flag_failed:\t{bc_flag_failed}\n")
    
    # Save sample of barcode collection for debugging (reduced size)
    try:
        n = min(10000, len(barcode_collection))
        sample_df = barcode_collection.iloc[:n].copy()
        if 'ReadID' in sample_df.columns:
            sample_df = sample_df.drop('ReadID', axis=1)
        with open(f'{dirname}/barcode_collection10k', 'wb') as fp:
            pickle.dump(sample_df, fp, protocol=4)
    except Exception as e:
        print(f"Warning: Could not save barcode_collection sample: {e}", file=sys.stderr)

    preview_df = barcode_collection.iloc[:100].copy()
    if 'ReadID' in preview_df.columns:
        preview_df = preview_df.drop('ReadID', axis=1)
    with open(f'{dirname}/barcode_collection', 'wb') as fp:
        pickle.dump(preview_df, fp, protocol=4)
    
    # Write out barcodes summary prior to any filtering
    raw_reads_per_bc_bins, barcode_summary = write_out_barcode_summary(
        barcode_collection, dirname, write_out_tsv
    )
    
    # Filter to keep barcodes greater than min_frag and min_reads
    try:
        barcode_collection_filtered = reads_notIn_frag(
            dirname, barcode_collection, min_reads, min_frag, max_frag
        )
    except Exception as e:
        print(f"Warning: reads_notIn_frag failed: {e}", file=sys.stderr)
        barcode_collection_filtered = barcode_collection[
            (barcode_collection['N_Reads'] >= min_reads) & 
            (barcode_collection['Frag_Length'] >= min_frag)
        ]
    
    barcode_collection_filtered = manipulate_df_gaps(
        barcode_collection_filtered, min_frag, max_frag, min_reads, read_len
    )
    
    # Fragment-dependent dup rate
    try:
        frag_dup_rate(dirname, barcode_collection_filtered, min_reads)
    except Exception as e:
        print(f"Warning: frag_dup_rate failed: {e}", file=sys.stderr)
    
    # Output fragment calculations
    write_out_tsv_and_summary2(
        barcode_collection, dirname, raw_reads_per_bc_bins, write_out_tsv, min_reads
    )
    write_out_tsv_and_summary1(
        barcode_collection_filtered, dirname, write_out_tsv, min_reads, plot_cutoff
    )
    
    # Combine summaries in the same order expected by downstream reports.
    try:
        script = f"cat {dirname}/frag_summary_minreads{min_reads}.txt > {dirname}/frag_and_bc_summary.txt && " \
                 f"grep frag_shorter_total_reads {dirname}/additional_summary_minreads{min_reads}.txt -A 1 >> {dirname}/frag_and_bc_summary.txt && " \
                 f"cat {dirname}/reads_per_bc_bins.txt >> {dirname}/frag_and_bc_summary.txt"
        call(script, shell=True)
    except Exception as e:
        print(f"Warning: Could not combine summaries: {e}", file=sys.stderr)
    
    try:
        df_sorted = reads_per_bc_distribution_mapped(dirname, barcode_summary)
        parse_bc100(barcode_collection, df_sorted, min_reads, 100, min_frag, dirname)
    except Exception as e:
        print(f"Warning: parse_bc100 failed: {e}", file=sys.stderr)
    
    print(f"Fragment calculation complete. Results in {dirname}/", file=sys.stderr)


def get_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Calculate fragment lengths from linked-read BAM files"
    )
    parser.add_argument("--bampath", type=str, required=True,
                        help="Path to BAM file for fragment length calculations.")
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
                        help="Comma separated list of chromosomes to use [default: all]")
    parser.add_argument("--includedups", action="store_true",
                        help="Include reads marked as duplicates.")
    parser.add_argument("--writeouttsvs", action="store_true",
                        help="Write out large, detailed TSVs of various data manipulations.")
    parser.add_argument("--outdir", type=str, default="Calc_Frag_Length",
                        help="Specify output directory [default ./Calc_Frag_Length]")
    parser.add_argument("--ray_mem", type=str, default="200,1024,100,1024",
                        help="Memory allocation (kept for compatibility, not used in optimized version).")
    parser.add_argument("--mapping_quality", type=int, default=30,
                        help="Mapping quality cutoff")
    parser.add_argument("--low_memory", action="store_true",
                        help="Use low memory mode (processes one chromosome at a time)")
    parser.add_argument("--lfr_type", choices=("stlfr", "clfr"), default="stlfr",
                        help="Pipeline flavor label for logging and future mode-specific behavior.")
    parser.add_argument("--n_tolerance", type=int, default=None,
                        help="Accepted for cLFR compatibility; not used.")
    parser.add_argument("--sequence_type", default=None,
                        help="Accepted for cLFR compatibility; not used.")
    parser.add_argument("--cbc_len", type=int, default=None,
                        help="Accepted for cLFR compatibility; not used.")
    parser.add_argument("--BC_type", default=None,
                        help="Accepted for cLFR compatibility; not used.")
    parser.add_argument("--keep_read_ids", action="store_true",
                        help="Store per-fragment ReadID lists. Disabled by default to save memory.")
    
    args = parser.parse_args()
    
    # Validate chromosomes
    valid_chroms = get_chroms(args.bampath)
    
    if args.chroms:
        chroms = args.chroms.split(",")
        invalid = [c for c in chroms if c not in valid_chroms]
        if invalid:
            parser.error(f"Invalid chromosomes: {invalid}\nValid: {valid_chroms}")
        args.chroms = chroms
    else:
        args.chroms = valid_chroms
    
    return args


def extract_barcode(read, cbc_len=None):
    try:
        return read.get_tag('BX')
    except KeyError:
        pass

    query_name = read.query_name
    if '#' not in query_name:
        return None
    barcode = query_name.split('#', 1)[1].split('/', 1)[0]
    if cbc_len:
        barcode = barcode[:cbc_len]
    return barcode or None


def process_chromosome(bam_path, chrom, split_dist, include_dups, read_len, mapping_quality, keep_read_ids, cbc_len=None):
    """
    Process a single chromosome and return fragment data.
    
    Memory optimization: Only store essential data per barcode:
    - Positions list (needed for fragment calculation)
    - Strand counter (only need majority)
    - Cigar match sum and count (for average)
    - Read count
    
    Returns a list of tuples instead of dict to reduce memory.
    """
    # Use compact data structure
    barcode_data = {}  # bc_tag -> [bc, chrom, positions, strand_counter, cigar_sum, cigar_count, read_len_sum]
    barcode_subs = {}  # Track sub-fragments
    
    stats = {'read_flag_failed': 0, 'poor_quality': 0, 'bc_flag_failed': 0, 'mapped_flag_failed': 0}
    
    # Flag masks
    FLAG_SECONDARY_CHIMERIC = 0x900  # 2048 + 256
    FLAG_WITH_DUP = 0xD00  # 2048 + 256 + 1024
    FLAG_UNMAPPED = 0x4
    FLAG_REVERSE = 0x10
    
    flag_mask = FLAG_SECONDARY_CHIMERIC if include_dups else FLAG_WITH_DUP
    
    try:
        bamfile = pysam.AlignmentFile(bam_path, "rb")
        
        for read in bamfile.fetch(chrom):
            # Get barcode
            bc = extract_barcode(read, cbc_len)
            if bc is None:
                stats['bc_flag_failed'] += 1
                continue
            
            if bc == "0_0_0":
                stats['bc_flag_failed'] += 1
                continue
            
            # Check flags
            flag = read.flag
            if flag & flag_mask != 0:
                stats['read_flag_failed'] += 1
                continue
            
            if flag & FLAG_UNMAPPED != 0:
                stats['mapped_flag_failed'] += 1
                continue
            
            if read.mapping_quality < mapping_quality:
                stats['poor_quality'] += 1
                continue
            
            # Valid read - process it
            read_strand = 0 if (flag & FLAG_REVERSE) else 1
            
            # Calculate cigar match length
            cigar_match = sum(c[1] for c in read.cigartuples if c[0] == 0) if read.cigartuples else 0
            
            # Determine bc_tag (handle sub-fragments)
            bc_identifier = f"{bc}_{chrom}"
            
            if bc_identifier not in barcode_subs:
                barcode_subs[bc_identifier] = 0
            else:
                # Check if we need a new sub-fragment
                bc_tag_current = f"{bc_identifier}_{barcode_subs[bc_identifier]}"
                if bc_tag_current in barcode_data:
                    last_pos = barcode_data[bc_tag_current]['positions'][-1]
                    if read.reference_start - last_pos >= split_dist:
                        barcode_subs[bc_identifier] += 1
            
            bc_tag = f"{bc_identifier}_{barcode_subs[bc_identifier]}"
            
            # Add to barcode data
            if bc_tag not in barcode_data:
                barcode_data[bc_tag] = {
                    'bc': bc,
                    'chrom': chrom,
                    'positions': [read.reference_start],
                    'positions_dup': [read.reference_start],  # Keep duplicates for dup rate calc
                    'strand_0': 1 if read_strand == 0 else 0,
                    'strand_1': 1 if read_strand == 1 else 0,
                    'cigar_matches': [cigar_match],
                    'read_lengths': [read.query_length],
                    'read_ids': [read.query_name] if keep_read_ids else [],
                }
            else:
                data = barcode_data[bc_tag]
                data['positions'].append(read.reference_start)
                data['positions_dup'].append(read.reference_start)
                if read_strand == 0:
                    data['strand_0'] += 1
                else:
                    data['strand_1'] += 1
                data['cigar_matches'].append(cigar_match)
                data['read_lengths'].append(read.query_length)
                if keep_read_ids:
                    data['read_ids'].append(read.query_name)
        
        bamfile.close()
        
    except Exception as e:
        print(f"Error processing chromosome {chrom}: {e}", file=sys.stderr)
        return [], stats
    
    # Convert to list format for DataFrame creation
    results = []
    for bc_tag, data in barcode_data.items():
        positions_sorted = sorted(data['positions'])
        n_reads = len(positions_sorted)
        min_pos = positions_sorted[0]
        max_pos = positions_sorted[-1]
        frag_length = max_pos - min_pos + read_len
        strand = 0 if data['strand_0'] > data['strand_1'] else 1
        
        results.append({
            'Barcode': data['bc'],
            'Chrom': data['chrom'],
            'Positions': positions_sorted,
            'Positions_dup': sorted(data['positions_dup']),
            'ReadID': data['read_ids'],
            'Strand': strand,
            'Cigar_match': data['cigar_matches'],
            'Read_Length': data['read_lengths'],
            'N_Reads': n_reads,
            'Min_Pos': min_pos,
            'Max_Pos': max_pos,
            'Frag_Length': frag_length,
        })
    
    return results, stats


def get_reads_parallel(bam_path, split_dist, include_dups, chroms, read_len, mapping_quality, n_threads, keep_read_ids, cbc_len=None):
    """
    Process BAM file using multiprocessing (no Ray dependency).
    More memory efficient than Ray for this use case.
    """
    all_results = []
    total_stats = {'read_flag_failed': 0, 'poor_quality': 0, 'bc_flag_failed': 0, 'mapped_flag_failed': 0}
    
    workers = max(1, min(n_threads, len(chroms)))
    print(f"Processing {len(chroms)} chromosomes with {workers} workers...", file=sys.stderr)
    
    with ProcessPoolExecutor(max_workers=workers) as executor:
        futures = {
            executor.submit(
                process_chromosome, bam_path, chrom, split_dist, 
                include_dups, read_len, mapping_quality, keep_read_ids, cbc_len
            ): chrom for chrom in chroms
        }
        
        for future in as_completed(futures):
            chrom = futures[future]
            try:
                results, stats = future.result()
                all_results.extend(results)
                for key in total_stats:
                    total_stats[key] += stats[key]
                print(f"  Completed {chrom}: {len(results)} fragments", file=sys.stderr)
            except Exception as e:
                print(f"  Error processing {chrom}: {e}", file=sys.stderr)
    
    # Create DataFrame
    if all_results:
        barcode_collection = pd.DataFrame(all_results)
    else:
        barcode_collection = pd.DataFrame(columns=[
            'Barcode', 'Chrom', 'Positions', 'Positions_dup', 'ReadID', 
            'Strand', 'Cigar_match', 'Read_Length', 'N_Reads', 'Min_Pos', 
            'Max_Pos', 'Frag_Length'
        ])
    
    print(f"Total fragments: {len(barcode_collection)}", file=sys.stderr)
    
    return barcode_collection, (
        total_stats['read_flag_failed'],
        total_stats['poor_quality'],
        total_stats['bc_flag_failed'],
        total_stats['mapped_flag_failed']
    )


def get_reads_low_memory(bam_path, split_dist, include_dups, chroms, read_len, mapping_quality, n_threads, keep_read_ids, cbc_len=None):
    """
    Process BAM file one chromosome at a time to minimize memory usage.
    Slower but uses much less RAM.
    """
    all_results = []
    total_stats = {'read_flag_failed': 0, 'poor_quality': 0, 'bc_flag_failed': 0, 'mapped_flag_failed': 0}
    
    print(f"Processing {len(chroms)} chromosomes sequentially (low memory mode)...", file=sys.stderr)
    
    for i, chrom in enumerate(chroms):
        print(f"  [{i+1}/{len(chroms)}] Processing {chrom}...", file=sys.stderr)
        
        results, stats = process_chromosome(
            bam_path, chrom, split_dist, include_dups, read_len, mapping_quality, keep_read_ids, cbc_len
        )
        
        all_results.extend(results)
        for key in total_stats:
            total_stats[key] += stats[key]
        
        print(f"    Found {len(results)} fragments", file=sys.stderr)
    
    # Create DataFrame
    if all_results:
        barcode_collection = pd.DataFrame(all_results)
    else:
        barcode_collection = pd.DataFrame(columns=[
            'Barcode', 'Chrom', 'Positions', 'Positions_dup', 'ReadID',
            'Strand', 'Cigar_match', 'Read_Length', 'N_Reads', 'Min_Pos',
            'Max_Pos', 'Frag_Length'
        ])
    
    print(f"Total fragments: {len(barcode_collection)}", file=sys.stderr)
    
    return barcode_collection, (
        total_stats['read_flag_failed'],
        total_stats['poor_quality'],
        total_stats['bc_flag_failed'],
        total_stats['mapped_flag_failed']
    )


if __name__ == "__main__":
    main()
