import pysam
import io
import subprocess
import os
import sys
import re
from Bio import SeqIO
import argparse
import gzip
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

parser = argparse.ArgumentParser()
parser.add_argument("--chr_name", type=str, required=False)
parser.add_argument("--module", type=str, required=False)
parser.add_argument("--seq_type", type=str, required=False)
parser.add_argument("--umi_len", type=int, required=True)
parser.add_argument("--input_bam", type=str, required=False)
parser.add_argument("--output_bam", type=str, required=False)
args = parser.parse_args()
module = args.module
chr_name = args.chr_name
seq_type = args.seq_type
BC_LEN = args.umi_len

def reformat_readid(chr_name, seq_type, input_bam, output_bam):
    if input_bam is None:
        input_bam = f"Make_Vcf/step3_hapcut/step1_modify_bam/data_sort.markdup_{chr_name}.bam"
    if output_bam is None:
        output_bam = f"Align/tmp/data_{chr_name}.name.bam"

    samfile = pysam.AlignmentFile(input_bam, "rb")
    outfile = pysam.AlignmentFile(output_bam, "wb", template=samfile)
    
    for read in samfile:
        readid = read.query_name
        umi = readid.split('#')[-1][:BC_LEN]
        pos = str(read.pos)
        if seq_type=='pe':
            read.query_name = '#'.join([umi,pos,read.query_name[:-2]])
        elif seq_type=='se':
            read.query_name = '#'.join([umi,pos,read.query_name])
        else:
            print('needs SEQUENCE_TYPE')
        outfile.write(read)
    outfile.close()
    return 

def fix_fasta(input_fasta, output_fasta):

    with open(input_fasta, "rt") as handle, open(output_fasta, 'w') as corrected:
        for record in SeqIO.parse(handle, "fasta"):
            new_seq = str(record.seq).replace("N", "").replace("n", "")
            if len(new_seq)>0:
                corrected.write('>'+record.id+'\n')
                corrected.write(new_seq+'\n')

def compareN(input_fasta, output_pdf):
    def parse_fasta(file_path):
        """
        Parses a FASTA file and yields (header, sequence) tuples.
        Handles both plain and gzipped FASTA files.
        """
        if file_path.endswith('.gz'):
            opener = gzip.open
        else:
            opener = open

        with opener(file_path, 'rt') as f:
            header = None
            sequence_lines = []
            for line in f:
                line = line.strip()
                if not line: # Skip empty lines
                    continue
                if line.startswith('>'):
                    if header is not None:
                        yield header, "".join(sequence_lines)
                    header = line
                    sequence_lines = []
                else:
                    sequence_lines.append(line)
            if header is not None: # Yield the last record
                yield header, "".join(sequence_lines)

    def count_n_in_sequences(fasta_file):
        """
        Counts 'N' occurrences in each sequence of a FASTA file.

        Args:
            fasta_file (str): Path to the input FASTA file.

        Returns:
            list: A list of integers, where each integer is the count of 'N's
                in a sequence.
        """
        n_counts = []
        print(f"Processing FASTA file: {fasta_file}")
        for i, (header, sequence) in enumerate(parse_fasta(fasta_file)):
            n_count = sequence.upper().count('N') # Convert to upper to count 'n' as well
            n_counts.append(n_count)
            if (i + 1) % 10000 == 0:
                print(f"  Processed {i + 1} sequences...")
        print(f"Finished processing {len(n_counts)} sequences.")
        return n_counts

    def plot_n_count_distribution(n_counts, output_pdf_path):
        """
        Plots the distribution of 'N' counts and saves it as a PDF.

        Args:
            n_counts (list): A list of 'N' counts.
            output_pdf_path (str): Path to save the output PDF file.
        """
        if not n_counts:
            print("No 'N' counts to plot. The input FASTA might be empty or contain no sequences.")
            return

        plt.figure(figsize=(10, 6))
        sns.histplot(n_counts, kde=True, stat="density", bins="auto", color='skyblue', edgecolor='black')
        plt.xlim(0, 100)
        plt.title('Distribution of N_Counts Per Sequence')
        plt.xlabel('Number_of_N in Sequence')
        plt.ylabel('Density')
        plt.grid(axis='y', alpha=0.75)

        # Add descriptive statistics to the plot (optional)
        mean_n = np.mean(n_counts)
        median_n = np.median(n_counts)
        max_n = np.max(n_counts)
        min_n = np.min(n_counts)
        total_sequences = len(n_counts)

        stats_text = (f'Total Sequences: {total_sequences}\n'
                    f'Mean N count: {mean_n:.2f}\n'
                    f'Median N count: {median_n:.0f}\n'
                    f'Min N count: {min_n}\n'
                    f'Max N count: {max_n}')
        plt.text(0.95, 0.95, stats_text, transform=plt.gca().transAxes,
                fontsize=10, verticalalignment='top', horizontalalignment='right',
                bbox=dict(boxstyle='round,pad=0.5', fc='wheat', alpha=0.5))

        plt.tight_layout()
        plt.savefig(output_pdf_path)
        print(f"Distribution plot saved to {output_pdf_path}")

    n_counts = count_n_in_sequences(input_fasta)
    plot_n_count_distribution(n_counts, output_pdf)
    return


if __name__ == "__main__":
    if module == 'reformat_readid':
        reformat_readid(chr_name, seq_type, args.input_bam, args.output_bam)
    elif module == 'fix_fasta':
        fix_fasta()
    elif module == 'compareN':
        input_fasta = 'down10/Align/consensus/consensus.fixRC.fasta'
        output_pdf = 'down10/Align/consensus/consensus.fixRC.fasta.N.pdf'
        compareN(input_fasta, output_pdf)

        
