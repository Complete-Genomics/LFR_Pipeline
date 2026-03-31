'''

usage:
chrom=chrX
bam=bx_sorted.bam
dir=
python $src --chrom ${chrom} --dir ${dir} --Nreads 10 --min_frag 300
'''
import pandas as pd
import pysam
import os
import argparse
import numpy as np
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime
from Bio import SeqIO
import pickle
import subprocess
import glob

def bed_sum(bed_file_merged):
    total_len_merged = 0
    with open(bed_file_merged, 'r') as f:
        for line in f:
            info = line.strip().split('\t')
            chrom = info[0]
            start = int(info[1])
            end = int(info[2])
            total_len_merged +=end-start+1
    # print(total_len_merged)
    return total_len_merged

def bc_bed_metrics(df, read_len, min_reads):
    # df = pd.read_csv(
    #     bed_file,
    #     sep='\t',
    #     header=None,
    #     names=['chrom', 'start', 'end', 'bc', 'Nreads']
    # )
    dirname = 'Align'

    df['frag_length'] = read_len+ df['end'] - df['start']

    unique_bc_count = df['bc'].nunique()
    fragment_count = len(df)
    total_reads_in_frag = df['N_Reads'].sum()
    avg_frag_per_bc = fragment_count / unique_bc_count
    avg_bc_read_count = total_reads_in_frag / unique_bc_count
    avg_frag_read_count = total_reads_in_frag / fragment_count
    avg_frag_length = df['frag_length'].mean()
    median_frag_length = df['frag_length'].median()


    bed_file = "Align/all_bc.bed"
    # df_bed.to_csv(bed_file, index=False, header=False, sep='\t')
    bed_file_merged = dirname + "/frag_minreads"+str(min_reads)+"_merged.bed"
    bed_sorted = dirname + "/frag_minreads"+str(min_reads)+"_sorted.bed"
    cmd = f'bedtools sort -i {bed_file} > {bed_sorted}'
    os.system(cmd) 
    bashCommand = f"bedtools merge -i {bed_sorted} > {bed_file_merged}"
    os.system(bashCommand)       
    ## frag coverage
    total_len_merged = bed_sum(bed_file_merged)

    GCA_000001405_15_GRCh38_len= 2934876545
    fragments_length_total_pct_genome = total_len_merged / GCA_000001405_15_GRCh38_len

    metrics = {
        "Unique_BC_Count": unique_bc_count,
        "Avg_Frag_Per_BC": avg_frag_per_bc,
        "Avg_BC_Read_Count": avg_bc_read_count,
        "Fragment_Count": fragment_count,
        "Avg_Frag_Length": avg_frag_length,
        "Median_Frag_Length": median_frag_length,
        "Avg_Frag_Read_Count": avg_frag_read_count,
        "Total_Reads_in_Frag": total_reads_in_frag,
        "Fragments_Length_Total_Pct_Genome": fragments_length_total_pct_genome,
        "Fragments_Length_Total": total_len_merged
    }
    
    output_file = dirname + "/frag_summary_minreads"+str(min_reads)+".txt"
    with open(output_file, "w") as f:
        for k,v in metrics.items():
            f.write(f"{k}: {v}\n")

    # print(f"Metrics saved to {output_file}")
    return 


def bc_bed():
    ## output all_bc.bed, N_reads as reads in frag
    ## end is end of mapped region, no need to add read_len
    df = pd.read_csv(
        'Align/all_N100_intChr.txt',
        sep='\t',
        header=None,
        names=['chrom', 'start', 'end', 'bc', 'strand', '']
    )
    

    grouped_df = df.groupby(['chrom', 'bc'], as_index=False).agg({
        'start': 'min',
        'end': 'max'
    })

    nreads_df = df.groupby(['chrom', 'bc']).size().reset_index(name='N_Reads')
    grouped_df = grouped_df.merge(nreads_df, on=['chrom', 'bc'], how='left')

    bed_file = 'Align/all_bc.bed'
    output_df = grouped_df[['chrom', 'start', 'end', 'bc', 'N_Reads']]
    if not grouped_df.empty:
        output_df.to_csv(bed_file, sep='\t', index=False, header=False)

    return output_df

def raw_bed_chr2bx(chrom):
    # chr='chr21'
    inputfile = f'Align/{chrom}.bed'
    outputfile = f'Align/{chrom}.chr2bx.bed'
    outfile = open(outputfile, 'w')
    with open(inputfile, 'r') as f:
        for line in f:
            info = line.strip().split('\t')
            bx = info[3].split('#')[1][:-2]
            new_line = '\t'.join([bx]+ info[1:3])
            outfile.write(new_line+'\n')
    outfile.close()
    return

def parallel_filterBED(chrom, read_len):
    split_distance, Nreads, min_frag = 30000, 100, 5000

    bx_dict = defaultdict(list)
    with open(infile, 'r') as f:
        lst =f.read().splitlines()

    return



# def merge_exon():
#     # if input_fasta=='' and output_fasta=='':
#     input_fasta=f'Align/all_tmp.fasta'
#     output_fasta=f'Align/all_N100.fasta'

#     bc_fa_dict = {}
#     bc_id_dict = {}
#     last_bc = 'NNNNN'
#     for record in SeqIO.parse(input_fasta, "fasta"):
#         info = record.id.split('::')
#         bc = info[0]
#         fa = record.seq
#         if bc != last_bc:
#             last_bc = bc
#             id = info[1]
#             bc_fa = fa
#             bc_fa_dict[bc]= bc_fa
#             bc_id = id
#             bc_id_dict[bc]=id
#         else:
#             bc_fa = bc_fa_dict[bc]+fa
#             bc_fa_dict[bc] =bc_fa
#             id = bc_id_dict[bc]+'+'+info[1]
#             bc_id_dict[bc]=id
    
#     with open(output_fasta, 'w') as f:
#         for bc in bc_fa_dict.keys():
#             id = '>'+bc+' '+bc_id_dict[bc]
#             f.write(f'{id}\n')
#             f.write(f'{bc_fa_dict[bc]}\n')
#     return

from collections import defaultdict
from Bio import SeqIO

def merge_exon(bc_reads_dict, input_fasta='Align/all_tmp.fasta', output_fasta='Align/all_N100.fasta'):
    ## output >bc exon_pos (reads count)
    bc_fa_dict = defaultdict(str)  
    bc_id_dict = defaultdict(str)  

    last_bc = None  

    for record in SeqIO.parse(input_fasta, "fasta"):
        info = record.id.split('::')
        bc = info[0]  
        fa = str(record.seq)  

        if bc != last_bc:
            bc_fa_dict[bc] = fa
            bc_id_dict[bc] = info[1]
        else:
            bc_fa_dict[bc] += fa
            bc_id_dict[bc] += '+' + info[1]

        last_bc = bc  

    with open(output_fasta, 'w') as f:
        for bc, seq in bc_fa_dict.items():
            id = f'>{bc} {bc_id_dict[bc]} {bc_reads_dict.get(bc, 0)}'
            f.write(f'{id}\n')
            f.write(f'{seq}\n')

    return

def filterNreads(cutoff):
    ## filter all_N100.fasta with Nread/bc NOT reads in frag !!
    cmd = f"awk '/^>/{'{'}split($0,a,\" \"); if(a[length(a)] > {cutoff}) p=1; else p=0{'}'} p' Align/all_N100.fasta > Align/all_N{cutoff}.fasta"
    subprocess.call(cmd, shell=True)
    return

def parse_bc(bc_path):
        bc_lst = []
        with open(bc_path) as f:
            for line in f:
                info = line.strip().split('\t')
                bc = info[0]
                bc_lst.append(bc)
        return bc_lst

def convert_stLFR(bc_lst):
    bc_dict = {}
    replacements = ['A', 'T', 'C', 'G']
    for bc in bc_lst:
        for i in range(len(bc)):
            for replacement in replacements:                    
                new_bc = bc[:i] + replacement + bc[i+1:]
                bc_dict[new_bc]= bc
    return bc_dict


def parallel_filterN100_dict(dir, chrom, split_distance, Nreads, min_frag, read_len):
    
    def process_read(read):
        bx = read.get_tag('BX')
        pos = read.pos
        

        nonlocal n, bx_seen, idx_lst, bx_idx_lst, read_lst, read_idx_lst, pos_last
        n += 1

        if bx != bx_seen:
            if idx_lst:
                if len(idx_lst[-1]) < Nreads or (read_len+idx_lst[-1][-1]-idx_lst[-1][0]<min_frag):
                    idx_lst = idx_lst[:-1]
                    read_lst = read_lst[:-1]
                if len(idx_lst) > 0:
                    # tmp_idx_lst = np.array([item for sublist in idx_lst for item in sublist])
                    # bx_idx_lst = np.append(bx_idx_lst, tmp_idx_lst)
                    tmp_read_lst = np.array([item for sublist in read_lst for item in sublist])
                    read_idx_lst = np.append(read_idx_lst, tmp_read_lst)
            bx_seen = bx
            idx_lst, read_lst = [], []
            idx_lst.append([pos])
            read_lst.append([read])
        else:
            if idx_lst:
                if pos - pos_last < split_distance:
                    idx_lst[-1].append(pos)
                    read_lst[-1].append(read)

                else:
                    if len(idx_lst[-1]) < Nreads:
                        idx_lst = idx_lst[:-1]
                        read_lst = read_lst[:-1]
                    else:
                        idx_lst.append([pos])
                        read_lst.append([read])
            else:
                idx_lst.append([pos])
                read_lst.append([read])
        pos_last = pos
        return

    threads = 10
    n = 0
    # min_frag = 5000
    read_idx, idx_lst, read_lst, pos_last = [], [[0]], [['x']], 0
    bx_idx_lst, read_idx_lst = np.array([]), np.array([])
    bx_seen = 'x'

    samfile = pysam.AlignmentFile(f"{dir}/data_sort.markdup.bx_sorted.{chrom}.bam", "rb")
    with ThreadPoolExecutor(max_workers=threads) as executor:
        for read in samfile:
            executor.submit(process_read, read)

    # print(f'len bx_idx_lst {len(bx_idx_lst)}')
    print(f'len read_idx_lst {len(read_idx_lst)}')
    # print(bx_idx_lst)

    samfile = pysam.AlignmentFile(f"{dir}/data_sort.markdup.bx_sorted.{chrom}.bam", "rb")
    outfile = pysam.AlignmentFile(f"{dir}/data_{chrom}.N100.bam", "wb", template=samfile)

    for read in read_idx_lst:
        outfile.write(read)

    outfile.close()
    return

def single_filterN100_dict(dir, chrom, split_distance, Nreads, min_frag, read_len):
    threads = 50
    n = 0
    # min_frag = 5000
    read_idx, idx_lst, read_lst, pos_last = [], [[0]], [['x']], 0
    bx_idx_lst, read_idx_lst = np.array([]), np.array([])
    bx_seen = 'x'

    samfile = pysam.AlignmentFile(f"{dir}/data_sort.markdup.bx_sorted.{chrom}.bam", "rb")
    print(f'start calculation {datetime.now()}')
    for read in samfile:
        n+=1
        bx = read.get_tag('BX')
        # bc = read.query_name
        # bx = bc.split('#')[1]
        pos = read.pos

        if bx != bx_seen:
            if idx_lst:
                if len(idx_lst[-1]) < Nreads or (read_len+idx_lst[-1][-1]-idx_lst[-1][0]<min_frag):
                    idx_lst = idx_lst[:-1]
                    read_lst = read_lst[:-1]
                if len(idx_lst) > 0:
                    # tmp_idx_lst = np.array([item for sublist in idx_lst for item in sublist])
                    # bx_idx_lst = np.append(bx_idx_lst, tmp_idx_lst)
                    tmp_read_lst = np.array([item for sublist in read_lst for item in sublist])
                    read_idx_lst = np.append(read_idx_lst, tmp_read_lst)
            bx_seen = bx
            idx_lst, read_lst = [], []
            idx_lst.append([pos])
            read_lst.append([read])
        else:
            if idx_lst:
                if pos - pos_last < split_distance:
                    idx_lst[-1].append(pos)
                    read_lst[-1].append(read)

                else:
                    if len(idx_lst[-1]) < Nreads:
                        idx_lst = idx_lst[:-1]
                        read_lst = read_lst[:-1]
                    else:
                        idx_lst.append([pos])
                        read_lst.append([read])
            else:
                idx_lst.append([pos])
                read_lst.append([read])
        pos_last = pos
    
    print(f'len read_idx_lst {len(read_idx_lst)}')
    print(f'n {n}')
    print(f'finish calculation {datetime.now()}')

    
    outfile = pysam.AlignmentFile(f"{dir}/data_{chrom}.N100.bam", "wb", template=samfile)
    
    for read in read_idx_lst:
        outfile.write(read)

    outfile.close()
    samfile.close()
    print(f'finish writing .bam {datetime.now()}')
    return

def replace_bx_chrom(chrom):
    # chr='chr21'
    inputfile = f'Align/{chrom}.merged.bed'
    outputfile = f'Align/{chrom}.intChr.bed'
    outfile = open(outputfile, 'w')
    with open(inputfile, 'r') as f:
        for line in f:
            info = line.strip().split('\t')
            bx = info[0]
            new_line = '\t'.join([chrom]+ info[1:]+[bx])
            outfile.write(new_line+'\n')
    outfile.close()
    return
    
def write_bam_pos(chrom):
    infile=f'Align/{chrom}.bam'
    bx_lst =[]
    bf = pysam.AlignmentFile(infile, 'rb')

    for read in bf:
        bc = read.query_name
        bx = bc.split('#')[1]
        bx_lst.append(bx)
        
    with open(f'Align/{chrom}.sam', 'w') as f:
        f.write("@HD"+'\t'+"VN:1.4" + "  "+ "GO:none\tSO:unsorted"+'\n')
        for bx in set(bx_lst):
            f.write("@SQ"+ "\t" +"SN:" + bx + "\t" + "LN:1"+'\n')
            
    tmp = f'Align/{chrom}.sam'
    chr_header = pysam.AlignmentFile(tmp).header
    samfile = pysam.AlignmentFile(f'Align/{chrom}.bam', "rb")
    with pysam.AlignmentFile(f'Align/{chrom}.bx.bam', "wb", header=chr_header) as outfile:
        for read in samfile:
            name=read.query_name
            a = pysam.AlignedSegment(outfile.header)
            bc = read.query_name
            bx = bc.split('#')[1]
            a.query_name = read.query_name
            a.pos = read.pos
            a.query_sequence = read.query_sequence
            a.reference_name = bx
            a.flag = read.flag
            a.reference_start = read.reference_start
            a.mapping_quality = read.mapping_quality
            a.cigar = read.cigar
            a.template_length = read.template_length
            a.query_qualities = read.query_qualities
            a.tags = read.tags
            
            outfile.write(a)

    return

def filter_2bc_count(cutoff, write2dict):
    _dict = defaultdict(int)
    total_reads = 0

    def count(log):
        with open(log, 'r') as f:
            [next(f) for x in range(4)]
            for line in f:
                info = line.strip().split('\t')
                bc = info[2]
                cnt = int(info[1])
                _dict[bc]+=cnt
    ## all logs, link other logs to current path for UMI data
    logs = glob.glob("*split_stat_read1.log")
    print(logs)
    for log in logs:
        # log = f'{batch}/split_stat_read1.log' 
        count(log)

    # out = open('./bc_count.txt', 'w')
    bc_reads_dict = {}
    for k,v in _dict.items():
        if v>=cutoff:
            # total_reads+=v
            bc_reads_dict[k]= v
    # out.close()
    # print(f'total reads >={cutoff}={total_reads}')
    if write2dict:
        with open('Align/bc_reads_dict', 'wb') as fp:
            pickle.dump(bc_reads_dict, fp)
    return bc_reads_dict

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--module", type=str, help='filterN100, merge_exon')
    parser.add_argument("--chrom", type=str, required=False)
    parser.add_argument("--dir", type=str, required=False)
    parser.add_argument("--Nreads", type=int, required=False)
    parser.add_argument("--threads", type=int, required=False)
    parser.add_argument("--min_frag", type=int, required=False)
    parser.add_argument("--read_len", type=int, required=False)
    parser.add_argument("--minreads_fasta", type=int, required=False)
    parser.add_argument("--min_reads", type=int, required=False)
    parser.add_argument("--name", type=str, required=False)
    args = parser.parse_args()
    
    module = args.module
    chrom=args.chrom
    threads=args.threads
    min_reads= args.min_reads

    if module=='filterN100':
        
        dir=args.dir
        Nreads=args.Nreads
        
        min_frag = args.min_frag
        read_len = args.read_len
        # ins = BX_bam(chr, bam_file)
        # ins.convert_back_BCsortedBAM()
        # ins.convert_BX()

        split_distance = 50000
        # Nreads = 100
        # min_frag = 5000
        single_filterN100_dict(dir, chrom, split_distance, Nreads, min_frag, read_len)
    
    elif module == 'filter_2bc_count':
        minreads_fasta = args.minreads_fasta
        filter_2bc_count(minreads_fasta, True)
    
    elif module =='merge_exon':
        name=args.name
        with open ('Align/bc_reads_dict', 'rb') as fp:
            bc_reads_dict = pickle.load(fp)
        # if name:
        merge_exon(bc_reads_dict)
        # else:
        #     merge_exon(f'{name}_tmp.fasta', f'{name}.fasta')

    elif module =='bc_bed':
        read_len = args.read_len
        output_df = bc_bed()
        bc_bed_metrics(output_df, read_len, min_reads)

    elif module =='filterNreads':
        cutoff = args.minreads_fasta
        filterNreads(cutoff)

    elif module=='write_pos':
        write_bam_pos(chrom)

    elif module=='replace_bx_chrom':
        replace_bx_chrom(chrom)

    elif module =='raw_bed_chr2bx':
        raw_bed_chr2bx(chrom)

    else:
        print('no module')