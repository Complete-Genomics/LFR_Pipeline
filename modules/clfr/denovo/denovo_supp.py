'''
1/ using mapped.fasta to fix denovo.fasta RC
python $src --module fix_fa_rc
2/ get RC of a seq
python $src --module get_rc_seq --fasta $fa
'''
import subprocess
import pickle
from Bio import SeqIO
from collections import defaultdict
import itertools
from statistics import mean, median
import seaborn as sns
import glob
import pysam
import numpy as np
import gzip
from sklearn.linear_model import LinearRegression
import warnings
warnings.filterwarnings("ignore")
import matplotlib.pyplot as plt
# matplotlib.use('Agg')
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--module", type=str, required=True)
parser.add_argument("--fasta", type=str, required=False)     
parser.add_argument("--adapter_seq", type=str, required=False)     
parser.add_argument("--flank_end", type=int, required=False)     
parser.add_argument("--outdir", type=str, required=False) 
parser.add_argument("--dir", type=str, required=False)  
parser.add_argument("--minreads_fasta", type=int, required=False)     
parser.add_argument("--batch_list", type=str, required=False)  

args = parser.parse_args()

###### combine 2 lanes 
def filter_2bc(batches, cutoff):

    #combine 2 lanes to increase coverage, for umi
    #output:  filtered_barcode_freq.txt

    _dict = defaultdict(int)
    total_reads = 0
    def count(batch):
        with open(batch, 'r') as f:
            [next(f) for x in range(4)]
            for line in f:
                info = line.strip().split('\t')
                bc = info[2]
                cnt = int(info[1])
                _dict[bc]+=cnt

    for i in batches:
        batch = DIR+ i+'/split_stat_read1.log' 
        # print(batch)
        count(batch)

    # cmd = 'mkdir -p denovo'
    # subprocess.call(cmd, shell=True)
    out = open('./filtered_barcode_freq.txt', 'w')
    for k,v in _dict.items():
        if v>=cutoff:
            total_reads+=v
            out.write('BX:Z:{}\n'.format(k))
    out.close()
    print(f'total reads >={cutoff}={total_reads}')

    return

def sum_bc(batches, cutoff):
    _dict = defaultdict(int)
    total_reads = 0
    def count(batch):
        with open(batch, 'r') as f:
            [next(f) for x in range(4)]
            for line in f:
                info = line.strip().split('\t')
                bc = info[2]
                cnt = int(info[1])
                _dict[bc]+=cnt

    for i in batches:
        batch = DIR+ i+'/split_stat_read1.log' 
        count(batch)
    for k,v in _dict.items():
        if v>=cutoff:
            total_reads+=v
    print(f'total reads >={cutoff}={total_reads}')

    return    



###### get direction, correct RC
def get_rc_seq(seq):
    complement = str.maketrans('ATCG', 'TAGC')
    complement_sequence = seq.translate(complement)
    rc_seq = complement_sequence[::-1]
    return rc_seq

def correct_direction_mapping(mapped_fasta, denovo_fasta, output, infile_paf):
    '''
    get frag direction from mapped.fasta (using r1 direction) to correct denovo.fasta direction (megahit output RC when palindrome found in assembly)
    '''
    ### get direction from map.frag.fasta
    neg, pos=0, 0
    mapped_bc_strand_dict = {}
    with open(mapped_fasta, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            strand = '-'
            description = record.description
            _id = record.id

            if '(+)' in description:
                strand = '+'
                pos+=1
            else:
                neg+=1
            # bc_len=15
            bc= _id[:CBC_LEN]
            mapped_bc_strand_dict[bc]=strand
    print(f'neg{neg}, pos{pos}')
    print(f'mapped_bc_strand_dict_len {len(mapped_bc_strand_dict)}')

    ## {bc: denovo_seq}
    bc_fasta_dict = {}
    dup_bc = []
    with open(denovo_fasta, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            _id = record.id
            bc_len =CBC_LEN+6
            bc = _id[:(bc_len)]
            bc_fasta_dict[bc]=(_id, str(record.seq))

    print(f'denovo_fasta_dict_len {len(bc_fasta_dict)}')
    # try:
    #     with open('denovo/dup_bc.txt', 'w') as f:
    #         for line in dup_bc:
    #             f.write(f"{line}\n")
    # except Exception as e:
    #     pass

    # {bc: strand}
    denovo_bc_strand_dict = defaultdict(int)
    with open(infile_paf, 'r') as f:
        [next(f) for x in range(7)]
        for line in f:
            try:
                info = line.strip().split('\t')
                bc = info[0][:(CBC_LEN+6)]
                strand = info[4] # error: index out of boundary
                denovo_bc_strand_dict[bc]=strand
            except Exception as e: 
                # print(e)
                print(line)

    # with open(infile_paf, 'r') as f:
    #     # skip header
    #     [next(f) for x in range(7)]
    #     for line in f:
    #         info = line.strip().split('\t')
    #         bc = info[0][:15]
    #         strand = info[4]
    #         denovo_bc_strand_dict[bc+strand]+=1

    ## {bc: strand}
    # Align/all_N100.fasta
    # >AAAAAAAAAGTATCC chr1:23557938-23558415(-)+chr1:23558934-23559019(-)+chr1:23559126-23559436(-)
    # mapped_denovo_bc_strand_dict = {}
    # with open(mapped_fasta, "r") as handle:
    #     for record in SeqIO.parse(handle, "fasta"):
    #         _id = record.id
    #         bc_len=15
    #         bc= _id[:bc_len]
    #         if '(-)' in _id:
    #             mapped_denovo_bc_strand_dict[bc]= '-'
    #         else:
    #             mapped_denovo_bc_strand_dict[bc]= '+'

    ## correct direction
    # output = 'denovo.fixRC.fasta'
    # with open(output,'w') as f:
    #     for bc in bc_fasta_dict:
    #         if denovo_bc_strand_dict[bc+'-'] or denovo_bc_strand_dict[bc+'+']:
    #             # if r1 denovo same strand
    #             if not denovo_bc_strand_dict[bc+'-']:
    #                 f.write('>'+(bc_fasta_dict[bc][0])+'\n')
    #                 f.write(bc_fasta_dict[bc][1]+'\n')
    #             # if r1 denovo same strand
    #             if not denovo_bc_strand_dict[bc+'+']:
    #                 f.write('>'+(bc_fasta_dict[bc][0])+'\n')
    #                 f.write(bc_fasta_dict[bc][1]+'\n')
    #             # if r1 denovo diff strand 
    #             if denovo_bc_strand_dict[bc+'+']>denovo_bc_strand_dict[bc+'-']:
    #                 seq = get_rc_seq(bc_fasta_dict[bc][1])
    #                 f.write('>'+(bc_fasta_dict[bc][0])+'\n')
    #                 f.write(seq+'\n')

    cnt_rc, frag_not_in_mapped =0, 0
    with open(output,'w') as f, open('./denovo/frag_not_in_mapped.fa', 'w') as f2:
        for bc in bc_fasta_dict:   
            try:
                # same direction
                if mapped_bc_strand_dict[bc[:CBC_LEN]]==denovo_bc_strand_dict[bc]:
                    _id,seq = bc_fasta_dict[bc]
                    f.write('>'+_id+'\n')
                    f.write(seq+'\n')
                    # pass
                else:
                    cnt_rc+=1
                    _id,seq = bc_fasta_dict[bc]
                    new_seq = get_rc_seq(seq)
                    f.write('>'+_id+'\n')
                    f.write(new_seq+'\n')
            except Exception as e: 
                frag_not_in_mapped+=1
                _id,seq = bc_fasta_dict[bc]
                f2.write('>'+_id+'\n')
                f2.write(seq+'\n')
                # f.write('>'+_id+'\n')
                # f.write(seq+'\n')
                
                # TODO: how to handle frag not in mapped.fasta
                # print(e)
                # print(bc)
    print(f'flip RC {cnt_rc}')
    print(f'denovo frag not in mapped.frag.fasta {frag_not_in_mapped}')
    return 

def correct_direction_adapter(out_tmp, output):
    cnt, cnt_rc, cnt_both = 0, 0, 0
    line_count = 0 
    max_lines = 1000 
    
    with open(out_tmp) as infile, open(output, 'w') as f:
        records = SeqIO.parse(infile, 'fasta')
        for record in records:
            _seq = str(record.seq)
            # needs correct rc
            if ADAPTER_RC in _seq[-FLANK_END:]:
                cnt_rc += 1
                new_seq = get_rc_seq(_seq)
                f.write('>'+record.id+'\n')
                f.write(new_seq+'\n')
                if ADAPTER in _seq[:FLANK_END]:
                    cnt_both += 1
            else:
                if ADAPTER in _seq[:FLANK_END]:
                    cnt += 1
                # no ADAPTER_RC, just print
                f.write('>'+record.id+'\n')
                f.write(_seq+'\n')
                      
            line_count += 1
    #         if line_count >= max_lines:
    #             break 
                
    print(f'adapter: {cnt}, adapter_RC: {cnt_rc}, both: {cnt_both}, total seq: {line_count}')
    return

def frag_10x_coverage(fasta):
    ## total frag
    frag_cnt = sum(1 for _ in SeqIO.parse(fasta, "fasta"))

    ## get cov10.bed
    cpus = 10
    bam = "Align/data.sort.removedup_rm000.bam"  
    output_bed = "Align/cov10.bed"  
    samtools="samtools"
    cmd = f"{samtools} depth -@ {cpus} {bam} | awk -v cov=10 '$3 < cov {{print $1\"\\t\"$2\"\\t\"$2}}' | bedtools merge > {output_bed}"
    subprocess.call(cmd, shell=True)
    print(cmd)

    ## get frag.10x.bed
    cmd = f"cat Align/*N100_intChr > Align/tmp"
    subprocess.call(cmd, shell=True)
    print(cmd)

    grep_cmd = f"grep -v ERCC Align/tmp > Align/frag.bed"
    subprocess.call(grep_cmd, shell=True)

    cmd = f"bedtools sort -i Align/frag.bed > Align/frag.sorted.bed"
    subprocess.call(cmd, shell=True)
    with open('Align/frag.sorted.bed', 'r') as infile, open('Align/frag.sorted.tag.bed', 'w') as outfile:
        for line in infile:
            fields = line.strip().split('\t')
            outfile.write('\t'.join(fields + [fields[1], fields[2]]) + '\n')

    bedtools_intersect_cmd = f"bedtools intersect -wao -a Align/frag.sorted.tag.bed -b Align/cov10.bed > Align/frag.cov10.bed"
    subprocess.call(bedtools_intersect_cmd, shell=True)

    ## count <10x region in frag
    gap_lst = defaultdict(int)
    cnt_5, cnt_3 = 0,0
    infile='Align/frag.cov10.bed'
    ## format chr1	11566	12506	GGTTATATTTCGG	+	+	chr1	11308	11602	36 
    # with open(infile, 'r') as f:
    #     for line in f:
    #         info = line.strip().split('\t')
    #         if int(info[9])>0:
    #             chrom = info[0]
    #             bc = info[3]
    #             gap_lst[bc+chrom] +=1
    ## format chr1	11566	11602(overlap s/e)	GGTTATATTTCGG	+	+	11566	12506(frag start/end)
    with open(infile, 'r') as f:
        for line in f:
            info = line.strip().split('\t')   
            chrom = info[0]
            bc = info[3]
            strand = info[4]
            gap_lst[bc+chrom] +=1
            ## bias to 5 or 3 prime
            o_start, o_end = int(info[1]),int(info[2])
            frag_start, frag_end = int(info[6]), int(info[7])
            frag_midpoint = (frag_end-frag_start)//2
            overlap_midpoint = (o_end-o_start)//2
            if frag_midpoint >=overlap_midpoint:
                if strand=='+':
                    cnt_5+=1
                else:
                    cnt_3+=1
            else:
                if strand=='+':
                    cnt_3+=1
                else:
                    cnt_5+=1

    per_frag_10x_cov = [v for v in gap_lst.values()]
    frag_with_gap_cnt = len(per_frag_10x_cov)
    with open(f"Align/mapping_frag_10x_cov.txt", 'a') as f:
        f.write(f'mean={round(mean(per_frag_10x_cov),4)}, median={median(per_frag_10x_cov)}\n')
        f.write(f'min={min(per_frag_10x_cov)}, max={max(per_frag_10x_cov)}\n')
        f.write(f'5 prime gap bias={cnt_5}, 3 prime gap bias={cnt_3}\n')
        f.write(f'frag with gap={frag_with_gap_cnt}, total frag={frag_cnt}\n')
        f.write(f'frag with gap ratio= {round(frag_with_gap_cnt/frag_cnt, 4)}\n')
    return

class frag_denovo(object):
    def __init__(self, fasta):
        self.assembly_len = []
        self.fasta = fasta
        self.bc_frag_dict = defaultdict(int)
        return 

    def assembly_len_distribution(self, CBC_LEN):
        try:
            infile = self.fasta
            records = SeqIO.parse(infile, 'fasta')
            for record in records:
                bc = record.id[0:(CBC_LEN)]
                self.bc_frag_dict[bc]+=1
                self.assembly_len.append(len(record.seq))

        except Exception as e: print(e)
        return

def metrics_basic(fasta, outdir):
    len_cutoff = 7000
    
    s = frag_denovo(fasta)
    s.assembly_len_distribution(CBC_LEN)  
    assembly_len = s.assembly_len
    assembly_len_filtered = [i for i in s.assembly_len if i<=len_cutoff]
    frag_cnt = len(assembly_len)

    ## method1, filter 1bc1frag, 1bcNfrag
    # too slow 
    # ONEfrag_id = [k for k,v in s.bc_frag_dict.items() if v==1]  
    # try:
    #     with open(f'{outdir}/denovo.fixRC.fasta') as original_file, open(f'{outdir}/denovo.fixRC.1bc1frag.fasta', 'w') as frag1, open(f'{outdir}/denovo.fixRC.1bcNfrag.fasta', 'w') as fragN:
    #         records = SeqIO.parse(original_file, 'fasta')
    #         for record in records:
    #             bc = record.id[0:(CBC_LEN)]
    #             if bc in ONEfrag_id:
    #                 SeqIO.write(record, frag1, 'fasta')
    #             else:
    #                 SeqIO.write(record, fragN, 'fasta')
    # except Exception as e: print(e)

    ## method2, faster, but frag dup because of chunk 200000 remains
    bc_frag_dict = defaultdict(list)
    with open(f'{outdir}/denovo.fixRC.fasta') as original_file, open(f'{outdir}/denovo.fixRC.1bc1frag.fasta', 'w') as frag1, open(f'{outdir}/denovo.fixRC.1bcNfrag.fasta', 'w') as fragN:

        current_id = 'space_holder'
        first_seq = 'AAA'
        bc_frag_dict[current_id].append((current_id, first_seq))
        # print(f'first id {bc_frag_dict}')
        records = SeqIO.parse(original_file, 'fasta')
        for record in records:
            bc = record.id[0:(CBC_LEN)]
            seq = str(record.seq)
            if bc !=current_id:
                # print(len(bc_frag_dict[current_id]))
                if len(bc_frag_dict[current_id])>1:
                    for v in bc_frag_dict[current_id]:
                        fragN.write(f'>{v[0]}\n{v[1]}\n')
                else:
                    # print(f'bc !=current_id{bc_frag_dict[current_id][0]}')
                    frag1.write(f'>{bc_frag_dict[current_id][0][0]}\n{bc_frag_dict[current_id][0][1]}\n')
                current_id = bc
                bc_frag_dict = defaultdict(list)
                bc_frag_dict[bc].append((str(record.id), seq))
            else:
                bc_frag_dict[bc].append((str(record.id), seq))
                print(len(bc_frag_dict))
                # print(f'bc ==current_id{bc_frag_dict}')

        # Handle the last barcode
        if len(bc_frag_dict[current_id]) > 1:
            for v in bc_frag_dict[current_id]:
                fragN.write(f'>{v[0]}\n{v[1]}\n')
        else:
            frag1.write(f'>{bc_frag_dict[current_id][0][0]}\n{bc_frag_dict[current_id][0][1]}\n')

    ## frag count per bc
    Nfrag_cnt_per_bc = [i for i in list(s.bc_frag_dict.values()) if i>1]  
    # sns.set(rc={'figure.figsize':(15,8)})
    # sns.set_style("whitegrid")
    # _fig = sns.distplot(Nfrag_cnt_per_bc,bins=50,kde=True)
    # plt.title(f"Nfrag_cnt_per_bc")
    # plt.savefig( f"{outdir}/Nfrag_cnt_per_bc.pdf")
    # plt.clf()
    try:
        with open(f"{outdir}/Nfrag_cnt_per_bc.txt", 'a') as f:
            f.write(f'{fasta}\n')
            f.write(f'count {len(Nfrag_cnt_per_bc)} frag with gap that <10x coverage\n')
            f.write(f'BC with multiple frag ratio={round(len(Nfrag_cnt_per_bc)/frag_cnt,4)}\n')            
            f.write(f'mean={round(mean(Nfrag_cnt_per_bc),4)}, median={median(Nfrag_cnt_per_bc)}\n')
            f.write(f'min={min(Nfrag_cnt_per_bc)}, max={max(Nfrag_cnt_per_bc)}\n')
    except Exception as e: print(e)

    ## plot denovo_frag_length_distribution
    sns.set(rc={'figure.figsize':(15,8)})
    sns.set_style("whitegrid")
    _fig = sns.distplot(assembly_len_filtered,bins=50,kde=True)
    plt.title(f"frag_length_distribution")
    plt.savefig( f"{outdir}/frag_length_distribution.pdf")
    plt.clf()
    plt.close()
    # print(f'count {len(assembly_len)} frag')
    # print(f'mean={round(mean(assembly_len),4)}, median={median(assembly_len)}')
    # print(f'min={min(assembly_len)}, max={max(assembly_len)}')
    try:
        with open(f"{outdir}/frag_length_distribution.txt", 'a') as f:
            f.write(f'{fasta}\n')
            # f.write(f'count {len(assembly_len)} frag\n')
            f.write(f'frag count {len(assembly_len)}, bc count {len(s.bc_frag_dict)}\n')
            f.write(f'mean={round(mean(assembly_len),4)}, median={median(assembly_len)}\n')
            f.write(f'min={min(assembly_len)}, max={max(assembly_len)}\n')
    except Exception as e: print(e)

    # 1bc1frag metrics
    s = frag_denovo(f'{outdir}/denovo.fixRC.1bc1frag.fasta')
    s.assembly_len_distribution(CBC_LEN)  
    assembly_len = s.assembly_len
    assembly_len_filtered = [i for i in s.assembly_len if i<=len_cutoff]
    frag_cnt = len(assembly_len)

    ## plot denovo_frag_length_distribution
    sns.set(rc={'figure.figsize':(15,8)})
    sns.set_style("whitegrid")
    _fig = sns.distplot(assembly_len_filtered,bins=50,kde=True)
    plt.title(f"frag_length_distribution")
    plt.savefig( f"{outdir}/frag_length_distribution.1bc1frag.pdf")
    plt.clf()
    plt.close()

    # print(f'count {len(assembly_len)} frag')
    # print(f'mean={round(mean(assembly_len),4)}, median={median(assembly_len)}')
    # print(f'min={min(assembly_len)}, max={max(assembly_len)}')
    try:
        with open(f"{outdir}/frag_length_distribution.txt", 'a') as f:
            f.write(f'{outdir}/denovo.fixRC.1bc1frag.fasta\n')
            # f.write(f'count {len(assembly_len)} frag\n')
            f.write(f'frag count {len(assembly_len)}, bc count {len(s.bc_frag_dict)}\n')
            f.write(f'mean={round(mean(assembly_len),4)}, median={median(assembly_len)}\n')
            f.write(f'min={min(assembly_len)}, max={max(assembly_len)}\n')
    except Exception as e: print(e)

    ## cumulative plot
    try:
        plt.figure(figsize=(8, 5))
        plt.hist(assembly_len_filtered, bins=100, histtype="stepfilled", alpha=0.6, color="C0", cumulative=True, density=True)
        plt.gca().yaxis.set_major_formatter(plt.FuncFormatter(lambda y, _: f"{y * 100:.0f}%"))
        plt.xlabel("Value Range")
        plt.ylabel("Cumulative Percentage (%)")
        plt.title("Cumulative Percentage Plot")
        plt.grid(True, linestyle="--", alpha=0.6)
        plt.savefig( f"{outdir}/cumulative_length_distribution.1bc1frag.pdf")
        plt.clf()
        plt.close()
    except Exception as e: print(e)
    
    return

def filter_plot(cutoff, fasta, batches):
    ## from Align/all_N100.fasta filter N200.fasta
    bc_reads_dict = filter_2bc_count(batches, 200, True)
    bc_list = [v for v in bc_reads_dict.values() if v>=cutoff ]

    ## filter N read/bc, plot distribution
    frag_len = []
    with open(fasta, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            _frag_len = len(record.seq)
            if record.id in bc_list:
                frag_len.append(_frag_len)

    ## plot denovo_frag_length_distribution
    sns.set(rc={'figure.figsize':(15,8)})
    sns.set_style("whitegrid")
    _fig = sns.distplot(assembly_len_filtered,bins=50,kde=True)
    plt.title(f"frag_length_distribution")
    plt.savefig( f"Align/frag_length_distribution.N{cutoff}.pdf")
    plt.clf()
    plt.close()
    return

def refilter_fasta(input_fasta, n_bp):
    output_fasta = f"{input_fasta[:-6]}.filter{n_bp}.fasta"
    with open(output_fasta, "w") as output_handle:
        for record in SeqIO.parse(input_fasta, "fasta"):
            if len(record.seq) > n_bp:
                SeqIO.write(record, output_handle, "fasta")
    return

# def associate_vc_reads(bam, vcf):

#     bamfile = pysam.AlignmentFile(bam, "rb")
#     vcf_in = pysam.VariantFile(vcf)

#     with open("vc_reads.bed", "w") as bed_out:
#         for rec in vcf_in.fetch():
#             chrom = rec.chrom
#             pos = rec.pos - 1  # convert 1-based VCF position to 0-based BED position

#             for pileupcolumn in bamfile.pileup(chrom, pos, pos + 1):
#                 if pileupcolumn.reference_pos == pos:
#                     for pileupread in pileupcolumn.pileups:
#                         if not pileupread.is_del and not pileupread.is_refskip:
#                             read_base = pileupread.alignment.query_sequence[pileupread.query_position]
#                             if read_base in rec.alts:  # Check if the read base matches any of the alternate alleles
#                                 bed_out.write(f"{chrom}\t{pos}\t{pos + 1}\t{pileupread.alignment.query_name}\t{read_base}\t{rec.ref}\n")

#     return


def associate_vc_reads2(bam, vcf):
    bamfile = pysam.AlignmentFile(bam, "rb")
    vcf_in = pysam.VariantFile(vcf)

    # Open BED file for writing
    with open("reads_supporting_indels.bed", "w") as bed_out:
        # Process each variant in the VCF file
        for rec in vcf_in.fetch():
            chrom = rec.chrom
            pos = rec.pos - 1  # VCF is 1-based, BED is 0-based

            # indel ?
            for pileupcolumn in bamfile.pileup(chrom, pos, pos + 1):
                if pileupcolumn.reference_pos == pos:
                    for pileupread in pileupcolumn.pileups:
                        if not pileupread.is_refskip:
                            read_base = pileupread.alignment.query_sequence[pileupread.query_position]
                            if read_base in rec.alts:  # Match alternate allele for INDEL
                                # Write to BED file: chrom, start, end, read_id, variant type (INDEL)
                                bed_out.write(f"{chrom}\t{pos}\t{pos + 1}\t{pileupread.alignment.query_name}\t{read_base}\t{rec.ref}\n")

    return

def remove_bz(fa, out_file):
    with open(out_file, 'w') as out: 
        for record in SeqIO.parse(fa, "fasta"):
            tmp = str(record.id).split(' ')[0]
            readid = tmp[:28]+tmp[-2:]
            out.write(readid+'\n')
            out.write(str(record.seq)[:128]+'\n')
    return

def fragLen_readsCount_corr():
    # outdir = './'
    # denovo_fasta = f'{outdir}/denovo.fixRC.1bc1frag.1k.fasta'
    # frag_len_dict = {}
    # records = SeqIO.parse(denovo_fasta, 'fasta')
    # for record in records:
    #     bc = str(record.id)[:CBC_LEN]
    #     frag_len = len(record.seq)
    #     frag_len_dict[bc] = frag_len

    # # batches = ['']
    # # bc_reads_dict = filter_2bc_count(batches, cutoff)
    # with open('frag_len_dict', 'wb') as fp:
    #     pickle.dump(frag_len_dict, fp)

    with open ('frag_len_dict', 'rb') as fp:
        frag_len_dict = pickle.load(fp)

    with open ('bc_reads_dict1k', 'rb') as fp:
        bc_reads_dict = pickle.load(fp)
    # bc_reads_dict = {}
    # for k,v in tmp.items():
    #     if v<1000:
    #         bc_reads_dict[k] = v
    # with open('bc_reads_dict1k', 'wb') as fp:
    #     pickle.dump(bc_reads_dict, fp)

    keys = set(bc_reads_dict.keys()).intersection(frag_len_dict.keys())
    x_values = np.array([frag_len_dict[key] for key in keys])
    y_values = np.array([bc_reads_dict[key] for key in keys])

    bins = [(1000, 2000), (2000, 3000), (3000, 4000), (4000, 5000), (5000, 6000)]

    ## Plot for each bin
    for i, (low, high) in enumerate(bins, start=1):
        bin_x = [x for x in x_values if low <= x < high]
        bin_y = [y for x, y in zip(x_values, y_values) if low <= x < high]
        _min, _mean, _max = min(bin_y), mean(bin_y), max(bin_y)
        plt.figure(figsize=(6, 4))
        plt.scatter(bin_x, bin_y, color='blue', alpha =0.5, s=30)

        plt.title(f"denovo_fragLen_readsCount_min{_min}_mean{_mean}_max{_max}")
        plt.xlabel(f"denovo_fragLen")
        plt.ylabel("readsCount")
        plt.legend()
        plt.grid(True)
        plt.savefig( f"{outdir}/{low}-{high}_denovo_fragLen_readsCount.pdf")
        plt.legend()
        plt.grid(True)

    # k1 = list(bc_reads_dict.keys())
    # k2 = list(frag_len_dict.keys())
    # print(k1[1:5])
    # print(k2[1:5])

    # m, b = np.polyfit(x_values, y_values, 1)
    # plt.plot(x_values, m * x_values + b, color='red')

    x_values_reshaped = x_values.reshape(-1, 1)
    # Predict y values based on the regression model
    reg_model = LinearRegression()
    reg_model.fit(x_values_reshaped, y_values)
    # Predict y values based on the regression model
    y_pred = reg_model.predict(x_values_reshaped)

    plt.figure(figsize=(8, 6))
    plt.scatter(x_values, y_values, color='blue', alpha =0.5, s=30)
    plt.plot(x_values, y_pred, color='red', label=f'Regression Line (y = {reg_model.coef_[0]:.2f}x + {reg_model.intercept_:.2f})')

    # plt.plot([min(x_values), max(x_values)], [min(x_values), max(x_values)], color='red', linestyle='--', label='y=x Line')
    plt.title("denovo_fragLen_readsCount")
    plt.xlabel("denovo_fragLen")
    plt.ylabel("readsCount")
    plt.legend()
    plt.grid(True)
    plt.savefig( f"{outdir}/denovo_fragLen_readsCount.pdf")
    return

def filter_2bc_count(batches, cutoff, write2dict):

    _dict = defaultdict(int)
    total_reads = 0
    def count(batch):
        with open(batch, 'r') as f:
            [next(f) for x in range(4)]
            for line in f:
                info = line.strip().split('\t')
                bc = info[2]
                cnt = int(info[1])
                _dict[bc]+=cnt

    for i in batches:
        batch = DIR+ i+'/split_stat_read1.log' 
        count(batch)
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
    module = args.module
    fasta = args.fasta
    outdir = args.outdir
    minreads_fasta = args.minreads_fasta
    batch_list = args.batch_list
    DIR = args.dir

    CBC_LEN = 15

    if module == 'metrics_basic':
        if fasta:
            metrics_basic(fasta, outdir)
        else:
            print('needs fasta')
    elif module =='fragLen_readsCount_corr':
        fragLen_readsCount_corr()

    elif module == 'sum_bc':
        batches = ['V350155392/L04/analysis', 'V350151742/L02/analysis']
        cutoff= 200
        sum_bc(batches, cutoff)

    elif module =='remove_bz':
        remove_bz('read.fa', 'read.nobz.fa')
    
    elif module == 'associate_vc_reads':
        bam, vcf = "sorted.bam", "variants.vcf.gz"
        associate_vc_reads2(bam, vcf)

    elif module =='filter_2bc':
        # batches = [f'{_dir}/L0{lane}/analysis' for lane in [1,2,3,4]]
        with open(batch_list) as f:
            batches =f.read().splitlines()
        filter_2bc(batches, minreads_fasta)

    ## see filterN100.py
    # elif module =='filter_2bc_count':
    #     # batches = ['/V350155392_V350151742/L02/analysis', '/V350155392_V350151742/L04/analysis']
    #     with open(batch_list) as f:
    #         batches =f.read().splitlines()
    #     minreads_fasta = 50
    #     filter_2bc_count(batches, minreads_fasta, False)

    elif module == 'filter_plot':
        fasta = "Align/all_N100.fasta"
        batches = ['/V350155392_V350151742/L02/analysis', '/V350155392_V350151742/L04/analysis']
        filter_plot(200, fasta, batches)
    elif module=='refilter_fasta':
        n_bp=400
        refilter_fasta(fasta, n_bp)

    elif module == 'get_rc_seq':
        info = fasta.split('.')
        out_file = info[0]+'_fixed.'+info[1]

        with open(fasta, "r") as handle, open(out_file, 'w') as out:
            for record in SeqIO.parse(handle, "fasta"):
                _seq = str(record.seq)
                seq_fixed = get_rc_seq(_seq)
                out.write('>'+record.id+'\n')
                out.write(seq_fixed+'\n')   

    elif module =='fix_fa_rc':
        ADAPTER = args.adapter_seq
        ADAPTER_RC = get_rc_seq(ADAPTER)
        FLANK_END = args.flank_end
        # denovo_paf = 'denovo.100.paf'
        # dir='/prod/hustor-04/zebra/ycai/results/V350136898/L01/analysis'
        _dir = f'./'
        infile_paf=_dir+'denovo/denovo.paf'
        denovo_fasta= _dir+'denovo/final_contigs_0.fa'
        out_tmp = _dir+'denovo/denovo.fixRCmapping.fasta'
        output0 = _dir+'denovo/denovo.fixAdapter.fasta'
        mapped_fasta = _dir+'/Align/all_N100.fasta'
    

    #     

        # ins = denovo_supp()
        # # chroms = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"]
        # chroms = ["chr21"] 
        # for chrom in chroms:
        #     try:
        #         ins.get_mapped_denovo_bc_strand_dict( mapped_fasta)
        #     except Exception as e: print(e)
        
        # mapped_denovo_bc_strand_dict = ins.mapped_denovo_bc_strand_dict
        # print(len(mapped_denovo_bc_strand_dict))
        # with open('mapped_denovo_bc_strand_dict', 'wb') as fp:
        #     pickle.dump(mapped_denovo_bc_strand_dict, fp)

        correct_direction_mapping(mapped_fasta, denovo_fasta, out_tmp, infile_paf)
        correct_direction_adapter(out_tmp, output0)
        input_notmapped ='denovo/frag_not_in_mapped.fa'
        output_notmapped = 'denovo/frag_not_in_mapped.fixAdapter.fa'
        correct_direction_adapter(input_notmapped, output_notmapped)

        output = _dir+'denovo/denovo.fixRC.fasta'
        cmd = f'cp {output0} {output} && cat {output_notmapped} >> {output}'
        subprocess.call(cmd, shell=True)

    elif module=='frag_10x_coverage':
        frag_10x_coverage(fasta)

    else:
        print('input a correct module')


