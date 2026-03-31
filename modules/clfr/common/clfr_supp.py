
import gzip

from collections import defaultdict
import pickle
import pysam
import argparse
import gzip
from Bio import SeqIO
import pandas as pd
from concurrent.futures import ThreadPoolExecutor, as_completed
import datetime


class count_mismatch(object):

    def __init__(self):
        self.dup = defaultdict(int)
        # self.mismatch_by1_perBC = defaultdict(list)
        return

    def set_attributes(self, all_bc_dict, gap_cutoff):
        self.all_bc_dict = all_bc_dict
        self.gap_cutoff = gap_cutoff
        return

    def count_dup(self, bc):
        all_bc_dict_key = list(self.all_bc_dict.keys())
        # mismatch_by1_perBC = defaultdict(list)
        for i in range(len(bc)):
            for base in ["A", "T", "C", "G"]:
                mismatch = bc[0:i]+base+bc[i+1:]
                if mismatch in all_bc_dict_key:
                    self.dup[mismatch]+=1
                    ## remove the bc, avoid double count
                    # self.all_bc_dict.pop(mismatch, None)
        return 

def unmapped_fq(_dir):
    # samtools view -b -f 4 $bam > unmapped.bam
    ## for R2

    infile=f'{_dir}/analysis/Align/unmapped.bam'
    id_lst = []

    bf = pysam.AlignmentFile(infile, 'rb')
    for r in bf:
        readid=str(r.query_name)
        id_lst.append(readid)
    # print(id_lst[:4])
    with open('unmapped_readid', 'wb') as fp:
        pickle.dump(id_lst, fp)

    # with open ('unmapped_readid', 'rb') as fp:
    #     id_lst = pickle.load(fp)
    id_lst = set(id_lst)

    cnt1, cnt2 = 0, 0
    cutoff = 200

    with gzip.open(f"{_dir}/analysis/data/split_read.2.fq.gz", "rt") as handle, open('unmapped.fq', 'w') as f1, open('mapped.fq', 'w') as f2:
        for record in SeqIO.parse(handle, "fastq"):
            readid=str(record.id)[:-2]

            # print(readid)
            if readid in id_lst: 
                cnt1+=1
                SeqIO.write(record, f1, 'fastq')
            else:
                cnt2+=1
                SeqIO.write(record, f2, 'fastq')
                
            if cnt1 >=cutoff and cnt2 >=cutoff:
                break
    return 


def count_singleton(chrom):
    
    with open('singleton_'+chrom, 'rb') as fp:
        singleton = pickle.load(fp)

    print(len(singleton))


## single-thread, cost 4.75hr to process bc_cnt 
def count_bc_single(infile, outfile):
    start=datetime.datetime.now()
    bc_cnt = defaultdict(int)
    
    with gzip.open(infile, "rt") as handle:
        for record in SeqIO.parse(handle, "fastq"):
            bc = str(record.id).split('#')[1][:-2]
            bc_cnt[bc]+=1
    # df = pd.DataFrame.from_dict(bc_cnt, orient='index').reset_index()
    # df.to_csv(outfile, index=False, header=False)
    print(f'time cost={datetime.datetime.now()-start}')
    return

## multi-threading, cost 10hr 
# def count_bc_chunk(records):
#     bc_cnt = defaultdict(int)
#     for record in records:
#         bc = str(record.id).split('#')[1][:-2]
#         bc_cnt[bc] += 1
#     return bc_cnt

# def merge_counts(counts_list):
#     merged_counts = defaultdict(int)
#     for counts in counts_list:
#         for bc, count in counts.items():
#             merged_counts[bc] += count
#     return merged_counts

# def process_file_in_chunks(infile, chunk_size=10000, max_workers=4):
#     start=datetime.datetime.now()
#     def chunk_generator(handle, chunk_size):
#         records = []
#         for record in SeqIO.parse(handle, "fastq"):
#             records.append(record)
#             if len(records) == chunk_size:
#                 yield records
#                 records = []
#         if records:
#             yield records

#     counts_list = []
#     with gzip.open(infile, "rt") as handle:
#         with ThreadPoolExecutor(max_workers=max_workers) as executor:
#             futures = [executor.submit(count_bc_chunk, chunk) for chunk in chunk_generator(handle, chunk_size)]
#             for future in as_completed(futures):
#                 counts_list.append(future.result())
    
#     total_counts = merge_counts(counts_list)
#     # print(total_counts)
#     # df = pd.DataFrame.from_dict(total_counts, orient='index').reset_index()
#     # df.to_csv(outfile, index=False, header=False)
#     print(f'time cost={datetime.datetime.now()-start}')
#     return

    
    # with open('all_bc_dict', 'rb') as fp:
    #     all_bc_dict = pickle.load(fp)
    # with open('singleton_'+chrom, 'rb') as fp:
    #     singleton = pickle.load(fp)

    # s = count_mismatch() 
    # s.set_attributes(all_bc_dict, gap_cutoff)
    # j = 0
    # for bc in singleton:
    #     print(j)
    #     j+=1
    #     s.count_dup(bc)

    # print(len(s.dup))
    # print(sum(list(s.dup.values())))
    # print(sum([i for i in list(s.dup.values()) if i>1]))
    # print(sum([i for i in list(s.dup.values()) if i>2]))

import itertools
import math
def combinations_pick(N, n):
    '''
    N = frag_per_bc, n = n_sample
    odds of N frag in different among n sample
    E150023620_E150023623
    30M BC, 160M frag (20*8samples)= 5.3 frag/BC
    Odds of 5 frag belong to 5 different sample = 0.21 = (pick5from8 * 5!)/8^5
    '''
    items = [i for i in range(n)] 
    combinations = list(itertools.combinations(items, N))
    a = len(combinations)

    def factorial_recursive(n):
        if n < 0:
            raise ValueError("Factorial is not defined for negative numbers.")
        if n == 0:
            return 1
        else:
            return n * factorial_recursive(n - 1)
            
    b = factorial_recursive(N)
    c = n**N
    res = a*b/c
    print(f'{a}, {b}, {c}')
    print(f'odds of {N} frag all in different samples among {n} samples: {res}')
    # return combinations


def main():
    gap_cutoff = 50000
    chrom = 'chr19'
    # count_singleton()

    with open ('bc_pos_dic_'+chrom, 'rb') as fp:
        bc_pos_dic = pickle.load(fp)
    
    total_reads= 0
    for k,v in bc_pos_dic.items():
        cnt = len(v)
        total_reads+=cnt
    print(total_reads)

if __name__ == "__main__":
    # main()

    module = 'demulx_bc_frag'

    if module =='demulx_bc_frag':
        ## demulx sample, how many sample/frag for 30M BC
        # frag as million
        frag_per_sample = 10
        n_sample = 10
        frag = n_sample * frag_per_sample
        ## bc as 30M
        frag_per_bc = int(frag/30)
        combinations_pick(frag_per_bc, n_sample)
    elif module == 'count_bc_single':
        ## count bc
        infile = '/prod/hustor-04/zebra/ycai/results/V350216292/L01/analysis/data/split_read.2.fq.gz'
        outfile = 'test2.csv'
        count_bc_single(infile, outfile)
    elif module =='process_file_in_chunks':
        process_file_in_chunks(infile, outfile)
    else:
        print('need a module name')