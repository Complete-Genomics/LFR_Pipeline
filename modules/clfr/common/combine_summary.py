'''
combine frag_and_bc_summary{}.txt
combine reads_per_bc_bins.txt
'''
import pandas as pd
from subprocess import call 
import argparse, os, sys

parser = argparse.ArgumentParser()
parser.add_argument("--chroms", type=str,
                        help="Comma seperated list of chromosomes to use [default: all]")
parser.add_argument("--outdir", type=str, default="Calc_Frag_Length",
                        help="Specify output directory [default ./Calc_Frag_Length]")

chroms = args.chroms
dirname = args.outdir


def combine_summary(chroms, dirname):
    if len(chroms)==1:
        script = f"cat {dirname}/frag_and_bc_summary{chrom[0]}.txt > {dirname}/frag_and_bc_summary.txt && cat {dirname}/reads_per_bc_bins{chroms[0]}.txt >> {dirname}/frag_and_bc_summary.txt"
        call(script, shell=True)
    
    else:
        infile=f'{dirname}/reads_per_bc_bins{chroms[0]}.txt'
        reads_per_bc_bins = pd.read_csv(infile, index_col=False, header=None, sep='\t')
        for chrom in chroms[1:]:
            infile=f'{dirname}/reads_per_bc_bins{chroms[i]}.txt'
            tmp =pd.read_csv(infile, index_col=False, header=None, sep='\t')[1]
            reads_per_bc_bins[1] = reads_per_bc_bins[1]+tmp
        out_stats = dirname + "/reads_per_bc_bins.txt"
        df.to_csv(out_stats, index=False, header=False)

        infile=f'{dirname}/frag_and_bc_summary{chroms[0]}.txt'
        reads_per_bc_bins = pd.read_csv(infile, index_col=False, header=None, sep='\t')
        for chrom in chroms[1:]:
            infile=f'{dirname}/frag_and_bc_summary{chroms[i]}.txt'
            tmp =pd.read_csv(infile, index_col=False, header=None, sep='\t')[1]
            reads_per_bc_bins[1] = reads_per_bc_bins[1]+tmp
            
        out_stats = dirname + "/frag_and_bc_summary.txt"
        df.to_csv(out_stats, index=False, header=False)

if __name__ == "__main__":
    # script = f"cat {dirname}/frag_summary_minreads{min_reads}.txt > {dirname}/frag_and_bc_summary.txt && cat {dirname}/reads_per_bc_bins.txt >> {dirname}/frag_and_bc_summary.txt"
    # call(script, shell=True)