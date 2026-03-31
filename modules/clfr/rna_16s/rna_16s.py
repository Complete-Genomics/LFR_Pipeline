"""
all functions used for 16s analysis, including frag denovo
output: 
summary stats
denovo_frag_length_distribution.pdf
"""
from pathlib import Path
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
# matplotlib.use('Agg')
from collections import Counter, defaultdict
import argparse
pd.set_option('display.precision', 4)
import json
import gzip
import linecache
import subprocess
from Bio import SeqIO
from statistics import mean, median
import plotly.graph_objs as go
import plotly
import pysam
import pickle

class rna_16s(object):
    def __init__(self, fasta):
        self.assembly_len = []
        self.fasta = fasta
        self.bc_frag_dict = defaultdict(int)
        return 

    def summary_denovo(self, bc):
        try:
            infile=f"quast/{bc}/summary/TXT/Total_length.txt"
            with open(infile) as f:
                lst =f.read().splitlines() 
            self.assembly_len.append(int(lst[1].split(' ')[2]))
        except Exception as e: print(e)
        return 

    def assembly_len_distribution(self, CBC_LEN):
        try:
            infile = self.fasta
            records = SeqIO.parse(infile, 'fasta')
            for record in records:
                bc = record.id[1:(CBC_LEN+1)]
                self.bc_frag_dict[bc]+=1
                self.assembly_len.append(len(record.seq))

        except Exception as e: print(e)
        return

'''
json2txt
parse_bc_index_txt
run_write2fq_index_txt
'''

def replace_chr_with_bx(name):
    inputfile = f'{name}'
    outfile = open(f'{name}.BXasChr', 'w')
    with open(inputfile, 'r') as f:
        for line in f:
            info = line.strip().split('\t')
            bx= info[3].split('#')[1].split('/')[0]
            new_line = '\t'.join([bx]+ info[1:])
            outfile.write(new_line+'\n')
    outfile.close()
    return

def replace_bx_with_chr(name, chrom):
    inputfile = f'{name}'
    outfile = open(f'{name}.intChr', 'w')
    with open(inputfile, 'r') as f:
        for line in f:
            info = line.strip().split('\t')
            bx = info[0]
            new_line = '\t'.join([chrom]+ info[1:]+[bx])
            outfile.write(new_line+'\n')
    outfile.close()
    return


def run_per_base_density(chrom):
    n_frag = 3
    num_1k, num_2k, num_3k, num_6k = 0,0,0,0
    try:
        bed = f'Align/{chrom}.chr2bx'
        bc_dict = defaultdict(list)
        
        with open(bed, 'r') as f:
            last_bc =next(f).split('\t')[0]
            for line in f:
                bc = line.strip().split('\t')[0]
                if bc !=last_bc and 'N' not in bc and "AAAAAAAAAAAAAA" not in bc:
                    total_exon_len = per_base_density(bc_dict[last_bc], num_1k, num_2k, num_3k, num_6k, n_frag, chrom)
                    last_bc =bc
                    ## for ERCC
                    if chrom.startswith('chr'):
                        if num_1k>=n_frag and num_2k>=n_frag and num_3k>=n_frag and num_6k>=2:
                            break
                        else:
                            if 1000<=total_exon_len<2000:
                                num_1k+=1
                            elif 2000<=total_exon_len<3000:
                                num_2k+=1
                            elif 3000<=total_exon_len<4000:
                                num_3k+=1
                            elif 6000<=total_exon_len<7000:
                                num_6k+=1
                    elif chrom.startswith('ERCC'):
                        if num_1k>=n_frag:
                            break
                        else:
                            if 1000<=total_exon_len<2000:
                                num_1k+=1

                else:
                    bc_dict[bc].append(line.split('\t'))
    except Exception as e: print(e)
    return

def per_base_density(bc_lst, num_1k, num_2k, num_3k, num_6k, n_frag, chrom):
    # bc_bed_file={bc}.chr2bx
    ## frag coverage density plot
    # chr22.chr2bx
    # AAGACTGCATGGAGG	17036195	17036332	E200008241_L01_C004R010_4869428#AAGACTGCATGGAGG/2	255	+

    # df = pd.read_csv(bc_bed_file, sep='\t', header=None, names=['chrom', 'start', 'end', 'id', 'score', 'strand'])
    df= pd.DataFrame(bc_lst, columns=['chrom', 'start', 'end', 'id', 'score', 'strand'])
    df['start']= df['start'].astype(int)
    df['end']= df['end'].astype(int)
    df = df.sort_values('start')
    bc = df.chrom[0]
    frag_start = df.start[0]
    # _mean = int(0.8*(list(df.end)[-1]-list(df.start)[0]))
    df.start = df.start+1-frag_start
    df.end = df.end+1-frag_start
    # df['len']= df.end-df.start
    total_exon_len = 0
    # df.to_csv(f'{bc}.tsv', index=False)
    lst = list(df.start)
    block = [0]
    # print(df.head())

    def plot_exon_frag_density(lst, df, block, total_exon_len):
        for i in range(len(lst)-1):
            if lst[i+1]-lst[i]>=1000:
                block.append(i+1)
        block.append(len(lst))
        
        m = len(block)
        # TODO: fix bin length, not bin_num=100
        if m-1>1:
            fig, axs = plt.subplots(1, m-1, figsize=(15, 5))
            for i in range(m-1):
                tmp= df.iloc[block[i]:block[i+1],]
                strand = df.strand[0]
                len_k = (tmp['end'].to_list()[-1]-tmp['start'].to_list()[0])/1000
                base_positions = pd.concat([pd.Series(row['chrom'],range(row['start'], row['end'])) for _, row in tmp.iterrows()])
                # print(base_positions.head())
                df2 = base_positions.index.tolist()
                bin_num = int(len_k*100)
                if bin_num<1:
                    bin_num=1
                sns.histplot(df2, ax=axs[i], bins=bin_num, stat='count')
                axs[i].set_title(f'pseudo-exon{i}-strand{strand}')

            plt.tight_layout()
            plt.savefig( f"Align/{chrom}_{bc}_base_density_len{total_exon_len}.pdf")
            plt.close()
        else:
            base_positions = pd.concat([pd.Series(row['chrom'],range(row['start'], row['end'])) for _, row in df.iterrows()])
            df2 = base_positions.index.tolist()
            strand = df.strand[0]
            len_k =  (df['end'].to_list()[-1]-df['start'].to_list()[0])/1000
            # total_exon_len = list(df.end)[-1]-list(df.start)[0]
            bin_num = int(len_k*100)
            plt.figure(figsize=(15, 6))
            sns.histplot(df2,bins=bin_num,stat='count')
            plt.xlabel('Base Position')
            plt.ylabel('Density')
            plt.title(f'Per-Base Density Plot strand{strand}')
            plt.savefig( f"Align/{chrom}_{bc}_base_density_len{total_exon_len}.pdf")
            plt.close()
        return 

    def merge_sum_bed(df):
        current_start = df.iloc[0]['start']
        current_end = df.iloc[0]['end']
        region_sum = 0
        for i in range(1, len(df)):
            row = df.iloc[i]
            if row['start'] <= current_end:
                current_end = max(current_end, row['end'])
                if i==len(df)-1:
                    region_sum += current_end-current_start
            else:
                region_sum += current_end-current_start
                current_start = row['start']
                current_end = row['end']
        return region_sum

    total_exon_len = merge_sum_bed(df)

    if 1000<=total_exon_len<2000 and num_1k<=n_frag:
        print(f'bc {bc}, num_1k {num_1k+1}, num_2k {num_2k}, num_3k {num_3k}, num_6k {num_6k}')
        plot_exon_frag_density(lst, df, block, total_exon_len)
    elif 2000<=total_exon_len<3000 and num_2k<=n_frag:
        print(f'bc {bc}, num_1k {num_1k}, num_2k {num_2k+1}, num_3k {num_3k}, num_6k {num_6k}')
        plot_exon_frag_density(lst, df, block, total_exon_len)
    elif 3000<=total_exon_len<4000 and num_3k<=n_frag:
        print(f'bc {bc}, num_1k {num_1k}, num_2k {num_2k}, num_3k {num_3k+1}, num_6k {num_6k}')
        plot_exon_frag_density(lst, df, block, total_exon_len)
    elif 6000<=total_exon_len<7000 and num_6k<=2:
        print(f'bc {bc}, num_1k {num_1k}, num_2k {num_2k}, num_3k {num_3k}, num_6k {num_6k+1}')
        plot_exon_frag_density(lst, df, block, total_exon_len)

    return total_exon_len


def parse_chr2bx(chrom, n_frag):
    ## frag coverage density plot
    bed = f'Align/{chrom}.chr2bx'
    bc_dict = defaultdict(list)
    
    with open(bed, 'r') as f:
        for line in f:
            bc = line.strip().split('\t')[0]
            if len(bc_dict)>n_frag:
                break
            else:
                bc_dict[bc].append(line)
    # remove last
    del bc_dict[bc]
    
    for bc in list(bc_dict.keys())[:n_frag]:
        # print(len(v))
        with open(f'Align/{bc}.chr2bx', 'w') as f:
            for line in bc_dict[bc]:
                f.write(f"{line}")
    # print(bc_dict.keys())
    return bc_dict.keys()

def pos_dict_SE(inputfile):
    # to get frag.fasta as ref for frag.denovo eval  
    # readID:pos dict for SE, as bxtools bx sort lost pos
    pos_dict = {}
    samfile = pysam.AlignmentFile(inputfile, "rb")
    for read in samfile:
        readid = read.query_name
        pos = read.pos
        pos_dict[readid] = pos

    with open('pos_dict_SE', 'wb') as fp:
        pickle.dump(pos_dict, fp)
    return pos_dict

def fix_chr_pos_PE(bx, chrom):
    ## switch pos for bxtools index .bed using mate pos, r1 mate pos is for r2
    out = open(bx+f'.fix.sam','w')
    with open(bx+'.tmp.sam', 'r') as f:
        lst =f.read().splitlines()
    for i in range(0,len(lst),2):
        r1 = lst[i].split('\t')
        r2 = lst[i+1].split('\t')
        r1[2]=chrom
        r2[2]=chrom
        r1[3] = r2[7]
        r2[3] = r1[7]
        out.write('\t'.join(r1)+'\n')
        out.write('\t'.join(r2)+'\n')
    out.close()
    return

def merge_exon_ref(bx): 
    fa_seq, name = [], []
    for read in SeqIO.parse(bx+'.fasta', "fasta"):
        name.append(read.id)
        fa_seq.append(str(read.seq))
    if len(fa_seq)>1:
        new_seq = ''.join(fa_seq)
        new_name = '-'.join(name)
        with open( bx+'.1line.fasta','w') as f:
            f.write('>'+new_name+'\n')
            f.write(new_seq+'\n')
    else:
        cmd=f'mv {bx}.fasta {bx}.1line.fasta '
        subprocess.call(cmd, shell=True)
    return



def sort_N50df(infile, min_reads, min_frag):
    # infile='Calc_Frag_Length_50000/chr20.frag_and_bc_dataframe.tsv'
    bx_lst = set()
    df = pd.read_csv(infile, index_col=False, sep='\t')
    df = df[df.Frag_Length>=5000 and df.N_Reads>=50]
    ## todo: set min_reads, min_frag
    for bx in df.iloc[:, 0]:
        bx_lst.add(bx)
    bx_lst_sorted = sorted(list(bx_lst))

    start =(len(bx_lst_sorted)-500)//2
    bx_lst = bx_lst_sorted[start:start+500]
    # remove dup bc
    bx_set = set(bx_lst)

    with open('bx_lst_sorted.txt', 'w') as f:
        for line in bx_set:
            f.write(f"{line}\n")

    return bx_lst

def add_bx(inputfile, outputfile):
    samfile = pysam.AlignmentFile(inputfile, "rb")
    outfile = pysam.AlignmentFile(outputfile, "wb", template=samfile)
    for read in samfile:
        bc = read.query_name
        bx = bc.split('#')[1]
        read.set_tag('BX',bx,replace=True)
        outfile.write(read)
    outfile.close()

def run_write2fq_index_txt(r1r2, bc_index, start, chunk):
    infile=f'data/split_read.{r1r2}.fq.gz'
    with gzip.open(infile, 'rb') as file:
        content = file.read()
        decoded_content = content.decode('utf-8')
        fq = decoded_content.split('\n')
        
    for bc in list(bc_index.keys())[start:start+chunk]:
        idx = bc_index[bc]
        write2fq_index(bc, idx, fq, r1r2)
    return

def parse_bc_index_txt(bc_index_txt, if_string, Nreads_cutoff, num_frag):
    '''
    get read index of bc with >100 reads
    '''
    bc_index ={}
    stop_reading = False  # Flag to stop reading
    
    with open(bc_index_txt,'r') as f:
        for line in f:
            if stop_reading:
                break
            
            if if_string==False:
                info = line.strip().split('\t')
                bc, idx = info[0], info[1].split(',')
                idx = [int(i) for i in idx]  
            else:
                info = line.strip().split('\t')
                bc, idx = info[0], info[1][1:-1].split(',') 
                idx = [int(i) for i in idx]   
            if len(idx)>=Nreads_cutoff:
                bc_index[bc]=idx

            if len(bc_index) > num_frag:
                stop_reading = True
                
    return bc_index

def json2txt(bc_index_json):
    '''
    bc_index.json to bc_index.txt
    '''
    f = open('bc_index.txt', 'w')
    with open(bc_index_json, 'r') as j:
        bc_index = json.loads(j.read())
    for k,v in bc_index.items():
        f.write(f"{k}\t{v}\n")
    f.close()
    return


def parse_bc4denovo(threads, infile):
    # infile = 'bc4denovo_n.txt'
    with open(infile) as f:
        lst =f.read().splitlines() 
    block_len=len(lst)//threads+1

    for i in range(threads):
        outfile=f'./fq/bc_{i+1}.txt'
        with open(outfile, 'w') as f:
            for bc in lst[i*block_len:(i+1)*block_len]:
                f.write(f"{bc}\n")
    # i=i+1
    # outfile=f'./fq/bc_{i+1}.txt'
    # with open(outfile, 'w') as f:
    #     for bc in lst[i*block_len:]:
    #         f.write(f"{bc}\n")
    return

def reads_per_bc_distribution_regardless_mapping(dirname):
    split_log = dirname +'/split_stat_read1.log'
    summary = pd.read_csv(split_log, header=None, index_col=None, sep='\t', skiprows=4)        
    summary.columns = ['index', 'Reads', 'BC']

    df = pd.DataFrame.from_dict(Counter(summary.Reads), orient='index').reset_index()
    df.columns = ['reads_per_bc','bc_count']
    df_sorted = df.sort_values('reads_per_bc')
    df_sorted['reads_count'] = df_sorted.bc_count*df_sorted.reads_per_bc
    total_reads = sum(df_sorted.reads_count)
    df_sorted['reads_count_pct'] =df_sorted['reads_count']/total_reads
    with open(dirname+f'additional_summary.txt', 'a') as f:
        print(f"total_reads_regardless_mapping:\t{total_reads}",file = f)
    df_sorted.to_csv(dirname + f"reads_per_bc_distribution_regardless_mapping.tsv", sep='\t', index=False, float_format='%.7f')
    return

def write2fq_index(bc, index_list, fq, r1r2):
    with gzip.open(f'./fq/{bc}.{r1r2}.fq.gz', 'wt') as f1:
        for idx in index_list:
            ## idx is seq row, line1=1
            f1.write('\n'.join(fq[idx-2:idx+2])+'\n')
    # with gzip.open(f'{bc}.2.fq.gz', 'wt') as f2:
    #     for idx in index_list:
    #         f2.write('\n'.join(fq1[idx-2:idx+2]))   
    return

def fix_fq(bc, r1r2):
    infile=f'./fq/{bc}.{r1r2}.fq.gz'
    with gzip.open(infile, 'rb') as file:
        content = file.read()
        decoded_content = content.decode('utf-8')
        fq = decoded_content.split('\n')

    sep='@V350170623'
    with gzip.open(f'./fq/{bc}.{r1r2}.fix.fq.gz', 'wt') as f1:
        f1.write('\n'.join(fq[:3])+'\n')
        
        for i in range(3,len(fq)-1):
            if i%3==0:
                lst = fq[i].split(sep)
                new_list = [lst[0], sep+lst[1]]
                f1.write('\n'.join(new_list)+'\n')
            else:
                f1.write(fq[i]+'\n')
        f1.write(fq[-1])
    
def run_write2fq_index(r1r2, bc_list):
    ### output bc.2.fq.gz
    infile='bc_index.json'
    with open(infile, 'r') as f:
        bc_dict = json.load(f)
    # infile='bc4denovo_n.txt'
    with open(bc_list, 'r') as f:
        bc_enough_reads = f.read().splitlines() 

    infile=f'data/split_read.{r1r2}.fq.gz'
    with gzip.open(infile, 'rb') as file:
        content = file.read()
        decoded_content = content.decode('utf-8')
        fq = decoded_content.split('\n')
    # fq = np.array(fq)numpy.core._exceptions.MemoryError: Unable to allocate 798. GiB for an array
    # infile='data/split_read.1.fq.gz'
    # with gzip.open(infile, 'rb') as file:
    #     content = file.read()
    #     decoded_content = content.decode('utf-8')
    #     fq1 = decoded_content.split('\n')

    # for bc in bc_enough_reads:
    #     index_list = bc_dict.get(bc,'')
    #     gzip_file='data/split_read.2.fq.gz'
    #     _line = []
    #     # with gzip.open(gzip_file, 'rt') as file:
    #     for idx in index_list:
    #         for i in range(idx-2, idx+2):
    #             _line.append(linecache.getline(gzip_file, i)+'\n')
    #     with gzip.open(f'{bc}.2.fq.gz', 'wt') as f2:
    #         for line in _line:
    #             f2.write(line+'\n') 
    # bc_enough_reads = ["GGCCAAGCTCGAGCT", "GGTACGGGATAGTTG"]

    for bc in bc_enough_reads:
        try:
            index_list = bc_dict.get(bc,'')
            write2fq_index(bc, index_list, fq, r1r2)
        except Exception as e: print(e)
    return


        
def pct_denovoFrag_ref(dir_name):
    rna_16s_len_dict ={
    "Bacillus_subtilis_16S_1": 1558,
    "Bacillus_subtilis_16S_2": 1558,
    "Bacillus_subtilis_16S_4": 1558,
    "Bacillus_subtilis_16S_5": 1558,
    "Bacillus_subtilis_16S_6": 1558,
    "Bacillus_subtilis_16S_7": 1558,
    "Enterococcus_faecalis_16S_1": 1562,
    "Escherichia_coli_16S_1": 1542,
    "Escherichia_coli_16S_2": 1542,
    "Escherichia_coli_16S_4": 1542,
    "Lactobacillus_fermentum_16S_1": 1568,
    "Lactobacillus_fermentum_16S_3": 1578,
    "Lactobacillus_fermentum_16S_5": 1568,
    "Listeria_monocytogenes_16S_1": 1552,
    "Listeria_monocytogenes_16S_2": 1552,
    "Listeria_monocytogenes_16S_3": 1552,
    "Pseudomonas_aeruginosa_16S_1": 1526,
    "Salmonella_enterica_16S_1": 1534,
    "Salmonella_enterica_16S_2": 1534,
    "Staphylococcus_aureus_16S_1": 1556,
    "Staphylococcus_aureus_16S_2": 1556,
    "Staphylococcus_aureus_16S_3": 1556, 
    }
    
    denovo_frag_len_lst = []
    denovo_frag_seen = set()
    infile=f'{dir_name}/quast/contigs_reports/minimap_output/all-contigs_max.coords'
    with open(infile, 'r') as f:
        for line in f:
            info = line.strip().split(' ')
            frag_name = info[12]
            ref_name = info[11]
            ref_len = rna_16s_len_dict[ref_name]
            frag_len = int(info[7])
            pct = frag_len/ref_len
            if frag_name not in denovo_frag_seen:
                denovo_frag_seen.add(frag_name)
                denovo_frag_len_lst.append(pct)
    
    print(f'mean={round(mean(denovo_frag_len_lst),4)}, median={median(denovo_frag_len_lst)}')
    print(f'min={min(denovo_frag_len_lst)}, max={max(denovo_frag_len_lst)}')
    return denovo_frag_len_lst

def coverage_bias(infile, name):
    ## visualize frag cover region on 16s ref 
    
    start_list, end_list =[], []
    denovo_frag_seen = set()
    with open(infile, 'r') as f:
        for line in f:
            info = line.strip().split('\t')
            strand = info[4]
            ref_name=info[5]
            frag_name = info[0]
            if frag_name not in denovo_frag_seen:
                denovo_frag_seen.add(frag_name)
                ref_len = int(info[6])
                # if strand=='+':
                frag_start, frag_end = int(info[7]),int(info[8])-int(info[7])
                start_pct, end_pct = frag_start/ref_len, frag_end/ref_len

                start_list.append(start_pct)
                end_list.append(end_pct)
                # elif strand=='-':
                #     frag_start, frag_end = (ref_len-int(info[8]))/ref_len,(ref_len-int(info[7]))/ref_len
                #     start_list.append(frag_start)
                #     end_list.append(frag_end)

    ## ref https://community.plotly.com/t/create-bar-between-start-and-end-time-on-chart/16750
    fig = go.FigureWidget()
    fig.add_bar(x=list(range(len(start_list))),
            y=end_list,
            base=start_list)
    fig.update_traces(marker=dict(color="black"))
    fig
    plotly.offline.plot(fig, filename=f'./{name}.coverage_bias.html')
    # fig.write_image( "coverage_bias.png")
    return

# def cdna_coverage_bias(infile, name):
#     # AAGACTGCATGGAGG	17036195	17036332	E200008241_L01_C004R010_4869428#AAGACTGCATGGAGG/2	255	+
#     start_list, end_list =[], []
#     denovo_frag_seen = set()
#     with open(infile, 'r') as f:
#         for line in f:
#             info = line.strip().split('\t')
#             strand = info[5]
#             # ref_name=info[5]
#             frag_name = info[0]
#             start = info[1]
#             end = info[2]
#             if frag_name not in denovo_frag_seen:
#                 denovo_frag_seen.add(frag_name)
#                 ref_len = int(info[6])
#                 # if strand=='+':
#                 # frag_start, frag_end = int(info[7]),int(info[8])-int(info[7])
#                 # start_pct, end_pct = frag_start/ref_len, frag_end/ref_len

#                 start_list.append(start)
#                 end_list.append(end)


#     ## ref https://community.plotly.com/t/create-bar-between-start-and-end-time-on-chart/16750
#     fig = go.FigureWidget()
#     fig.add_bar(x=list(range(len(start_list))),
#             y=end_list,
#             base=start_list)
#     fig.update_traces(marker=dict(color="black"))
#     fig
#     plotly.offline.plot(fig, filename=f'./{name}.coverage_bias.html')
#     # fig.write_image( "coverage_bias.png")
#     return

### all.fq
def metrics_basic(fasta, outdir, cutoff):
    
    s = rna_16s(fasta)
    s.assembly_len_distribution(CBC_LEN)   
    len_cutoff = 7000
    assembly_len = s.assembly_len
    assembly_len_filtered = [i for i in s.assembly_len if i<=len_cutoff]
    # ## plot denovo_frag_length_distribution
    sns.set(rc={'figure.figsize':(15,8)})
    sns.set_style("whitegrid")
    _fig = sns.distplot(assembly_len_filtered,bins=50,kde=True)
    plt.title(f"frag_length_distribution")
    plt.savefig( f"{outdir}/frag_length_distribution_N{cutoff}.pdf")
    plt.clf()
    # print(f'count {len(assembly_len)} frag')
    # print(f'mean={round(mean(assembly_len),4)}, median={median(assembly_len)}')
    # print(f'min={min(assembly_len)}, max={max(assembly_len)}')
    try:
        with open(f"{outdir}/frag_length_distribution.txt", 'a') as f:
            f.write(f'{fasta}\n')
            f.write(f'frag count {len(assembly_len)}, bc count {len(s.bc_frag_dict)}\n')
            f.write(f'mean={round(mean(assembly_len),4)}, median={median(assembly_len)}\n')
            f.write(f'min={min(assembly_len)}, max={max(assembly_len)}\n')
    except Exception as e: print(e)
    return

def run_per_base_density_bc(chrom, bc):
    num_1k, num_2k, num_3k,num_6k, n_frag = 1, 1, 1, 1, 1
    try:
        bed = f'Align/{chrom}.{bc}.chr2bx'
        bc_dict = defaultdict(list)
        
        with open(bed, 'r') as f:
            for line in f:
                bc_dict[bc].append(line.split('\t'))
        # print(bc_dict[bc])
        total_exon_len = per_base_density(bc_dict[bc], num_1k, num_2k, num_3k,num_6k,n_frag, chrom)
        # print(total_exon_len)    
    except Exception as e: print(e)
    return

# ### frag position on 16s ref
def metrics_ref_zymo():
    pct_denovoFrag_ref('./')
    
    # infile = f'{dir_name}/quast/contigs_reports/minimap_output/all-contigs_max.coords_tmp'
    infile='./quast/all_contigs/contigs_reports/minimap_output/all-contigs.coords_tmp'
    name='all_contigs'
    coverage_bias(infile, name)
    return

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--module", type=str, required=False)
    parser.add_argument("--fasta", type=str, required=False)
    parser.add_argument("--chrom", type=str, required=False)
    parser.add_argument("--bx", type=str, required=False)
    parser.add_argument("--outdir", type=str, required=False)
    parser.add_argument("--minreads_fasta", type=int)
    args = parser.parse_args() 
    module = args.module
    fasta = args.fasta
    bx = args.bx
    chrom = args.chrom
    outdir = args.outdir
    CBC_LEN = 15
    # r1r2 = args.r1r2
    # print(module)
    
    n=3
    ### get read/bc distribution
    # reads_per_bc_distribution_regardless_mapping(dirname)

    ### fix .fq no new line bug
    # infile = 'bc4denovo_n.txt'
    # with open(infile) as f:
    #     lst =f.read().splitlines() 
    # # for bc in lst:
    # #     fix_fq(bc, 1)
    # #     fix_fq(bc, 2)

    # ### quast per bc
    # for bc in lst[-52:-2]:
    #     run_quast(bc)

    ### summary plot
    # s = rna_16s()
    dir_name='./'
    ## per bc
    # infile = 'bc4denovo_n.txt'
    # with open(infile) as f:
    #     lst =f.read().splitlines() 
    # for bc in lst:
    #     s.summary_denovo(bc)
    
    if module=='metrics_basic':
        # fasta = 'all.contigs_max.fasta'
        if fasta:
            cutoff = args.minreads_fasta
            metrics_basic(fasta, outdir, cutoff)
        else:
            print('needs fasta')
    
    elif module=='metrics_ref_zymo':
        metrics_ref_zymo()

    elif module=='frag_coverage':
        # chrom='chr21'
        # n_frag = 50
        # bc_lst = parse_chr2bx(chrom, n_frag)
        # # bc_lst = ['AAATAAGGGGAATGG']
        # for bc in bc_lst:
        #     print(bc)
        #     bc_bed_file = f'Align/{bc}.chr2bx'
        #     per_base_density(bc_bed_file)
        #     replace_bx_with_chr(f'Align/{bc}.chr2bx', chrom)
        for _chrom in [chrom, 'ERCC-00130', 'ERCC-00096', 'ERCC-00116', 'ERCC-00025']: 
            run_per_base_density(_chrom)
        cmd = f'touch Align/frag_coverage_done'
        subprocess.call(cmd, shell=True)

    elif module=='fix_chr_pos_PE':
        fix_chr_pos_PE(bx, chrom)

    elif module=='add_bx':
        inputfile = f'Calc_Frag_Length_50000/data_chr22.N100.bam'
        outputfile = f'Calc_Frag_Length_50000/data_chr22.N100.bx.bam'
        add_bx(inputfile, outputfile)

    elif module == "replace_bx_with_chr":
        outfile = open(f'Align/{bx}.merged.intChr', 'w')

        with open(f'Align/{bx}.merged', 'r') as f:
            for line in f:
                # print(line)
                info = line.strip().split('\t')
                bx = info[0]
                new_line = '\t'.join([chrom]+ info[1:]+[bx])
                outfile.write(new_line+'\n')
                # print(new_line)       
        outfile.close()
    elif module =='frag_coverage_1bc':
        ### plot for chrom.bc.chr2bx
        chrom = 'ERCC-00116'
        bc='CCGCCGACCCAATAG'
        run_per_base_density_bc(chrom, bc)
    else:
        print('no module')


