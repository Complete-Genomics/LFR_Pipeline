'''
# supporting functions for LFR data
# usage
# python $src --module combine_split_log
'''
import gzip
import numpy as np
import math
from collections import defaultdict
import pickle
from statistics import mean, median, stdev
import pandas as pd
from Bio import SeqIO
import sys
import re
import argparse
import os
from collections import Counter
import matplotlib.pyplot as plt
import seaborn as sns

class combine_demulx(object):
    def __init__(self):
        self.bc_dict = defaultdict(int)
        self.Reads_pair_num = 0
        self.Reads_pair_num_split = 0
        return

    def combine_split_log(self, log_list, minreads_fasta, module):
        '''
        combine split_stat_read1.log from different lane (#47, #73 etc.), count BC for 47_73_74_75
        '''
        for log in log_list:
            line = []
            with open(log, 'r') as f:
                for i in range(4):
                    line.append(f.readline().strip())
                self.Reads_pair_num += int(line[2].split(' ')[3])
                self.Reads_pair_num_split += int(line[3].split(' ')[3])
            
        for log in log_list:
            with open(log, 'r') as f:
                [next(f) for x in range(4)]
                for line in f:
                    info = line.strip().split('\t')
                    bc = info[2]
                    cnt = int(info[1])
                    if cnt > minreads_fasta:
                        self.bc_dict[bc] +=cnt
        
        bc_types = len(self.bc_dict)
        with open('split_stat_read1.log', 'a') as f:
            f.write('Barcode_types = 1536 * 1536 * 1536 \n')
            f.write(f'Real_Barcode_types_readsMoreThan{minreads_fasta} = {bc_types} ({bc_types/3623878656} %)\n')
            f.write(f'Reads_pair_num  = {self.Reads_pair_num}\n')
            f.write(f'Reads_pair_num(after split) = {self.Reads_pair_num_split} ({self.Reads_pair_num_split/self.Reads_pair_num} %)\n')
        
        ## merge n lanes, see filterN100.py filter_2bc_count
        if module == 'combine_split_log_exon2fasta':
            with open('Align/bc_reads_dict', 'wb') as fp:
                pickle.dump(self.bc_dict, fp)

        return

def across_sample(lanes):
    ## use frags_reads_per_bc.tsv to get frags_reads_per_bc across multiple lane (no split by sample BC)
    frags_reads_per_bc_dict = defaultdict(int)
    global frag_count
    global bc_count
    frag_count , bc_count= 0, 0
    

    def count_frag(lane):
        bc_set = set()
        global frag_count 
        global bc_count
        cnt = 0 
        infile = f'{lane}/analysis/Calc_Frag_Length_50000/frag_and_bc_dataframe.tsv'
        with open(infile, 'r') as f:
            next(f)
            for line in f:
                cnt+=1
                info = line.strip().split('\t')
                bc = info[0]
                # frags_reads_per_bc_dict[bc]+=cnt
                bc_set.add(bc)
                # frag_count +=cnt

        unique_bc = len(bc_set)
        print(f'id:{lane}')
        print(f'frag_cnt:{cnt}')
        print(f'unique_bc_cnt:{unique_bc}')
        frag_count+=cnt
        bc_count+=unique_bc

    for lane in lanes:
        count_frag(lane)
    
    # with open('frag_and_bc_dataframe', 'wb') as fp:
    #     pickle.dump(frags_reads_per_bc_dict, fp)

    # frag_cnt_val = frags_reads_per_bc_dict.values()
    # ave_frag_bc = mean(frag_cnt_val)
    # median_frag_bc = median(frag_cnt_val)
    # print(f'ave_frag_bc:{ave_frag_bc}')
    # print(f'frag_count:{frag_count}')
    # print(f'median_frag_bc:{median_frag_bc}')

    ## bc count across samples
    # bc_cnt_set = set()
    # def count_bc(lane):
    #     bc_set = set()
    #     infile = f'{lane}/analysis/Calc_Frag_Length_50000/frags_reads_per_bc.tsv'
    #     with open(infile, 'r') as f:
    #         next(f)
    #         for line in f:
    #             info = line.strip().split('\t')
    #             bc = info[0]
    #             bc_cnt_set.add(bc)
    #             bc_set.add(bc)
    #     print(f'unique_bc:{len(bc_set)}')
    
    # for lane in lanes:
    #     count_bc(lane)

    # bc_count = len(bc_cnt_set)
    print(f'bc_count:{bc_count}')
    print(f'frag_count:{frag_count}')
    print(f'ave_frag_bc:{frag_count/bc_count}')

    return

def calculate_expected_beads(num_beads, per_sample_frag, num_sample, ):
    ## for demulx, apply 30M beads to 8 samples (20M frag each, total 160M), how many beads are expected per sample
    total_frag = per_sample_frag * num_sample
    # 每个珠子结合的 fragment 数量（假设服从均匀分布 U(1, 10)）
    k_values = np.random.randint(1, 11, size=num_beads)  # len(k_values)= num_beads

    # 计算珠子落入某个 sample 的概率
    # 对于某个珠子，其 k 个 fragment 都不属于某个 sample 的概率为 (1 - 1/S)^k
    # 因此，至少一个 fragment 属于某个 sample 的概率为 1 - (1 - 1/S)^k
    P = 1 - (1 - 1/num_sample) ** k_values

    # 每个 sample 能结合的珠子种类数（期望值）
    beads_per_sample = np.sum(P)

    print(f"each sample expect to see {int(beads_per_sample):,} beads")
    return 

def frag_bed2summary(minreads):
    """
    V350152037, get frag_and_bc_summary.txt from frag_minreads4.bed
    chr1	9999	15145	1226_211_1080	9
    """
    infile = f'Calc_Frag_Length_50000/frag_minreads{minreads}.bed'

    test_barcodes = pd.read_csv(
        infile,
        sep='\t',
        index_col=False,
        names=['Chrom', 'start', 'end', 'Barcode', 'N_Reads'],
        dtype={'start': int, 'end': int, 'N_Reads': int}
    )

    test_barcodes['Frag_Length'] = test_barcodes['end']-test_barcodes['start']
    unique_barcode_count  = test_barcodes['Barcode'].nunique()
    fragment_count = test_barcodes.shape[0]
    frags_per_bc = fragment_count/unique_barcode_count
    avg_frag_len = test_barcodes['Frag_Length'].sum()/fragment_count
    med_frag_len = test_barcodes['Frag_Length'].median()
    avg_frag_read_count = test_barcodes['N_Reads'].sum()/fragment_count
    avg_bc_read_count = test_barcodes['N_Reads'].sum()/unique_barcode_count
    total_frag_read_count = test_barcodes['N_Reads'].sum()

    # min_reads = 4
    bed_file_merged =  f"Calc_Frag_Length_50000/frag_minreads"+str(minreads)+"_merged.bed"

    total_len_merged = bed_sum(bed_file_merged)
    GCA_000001405_15_GRCh38_len= 2934876545
    total_len_merged_pct = total_len_merged/GCA_000001405_15_GRCh38_len

    outfile = 'Calc_Frag_Length_50000/frag_and_bc_summary1.txt'

    # out_stats = dirname + "/frag_summary_minreads"+str(min_reads)+".txt"
    with open(outfile, "w") as frag_stats:
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
            f"Ave_read_per_kbp:\t{1000*avg_frag_read_count/med_frag_len}",
            file = frag_stats)

    # Generate distribution plot for fragments per barcode
    plt.style.use('default')
    plt.rcParams['figure.figsize'] = (12, 8)
    plt.rcParams['font.size'] = 12

    # Count fragments per barcode
    frag_counts_per_bc = test_barcodes['Barcode'].value_counts().sort_values(ascending=False)

    # Calculate statistics
    frag_per_bc_stats = frag_counts_per_bc.describe()

    # Create figure with subplots
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(16, 12))

    # 1. Histogram of fragments per barcode
    ax1.hist(frag_counts_per_bc.values, bins=50, alpha=0.7, color='skyblue', edgecolor='black')
    ax1.set_xlabel('Number of Fragments per Barcode')
    ax1.set_ylabel('Number of Barcodes')
    ax1.set_title('Distribution of Fragments per Barcode')
    ax1.grid(True, alpha=0.3)

    # 2. Box plot
    ax2.boxplot(frag_counts_per_bc.values, vert=False, patch_artist=True,
                boxprops=dict(facecolor='lightgreen', color='darkgreen'),
                medianprops=dict(color='red', linewidth=2))
    ax2.set_xlabel('Number of Fragments per Barcode')
    ax2.set_title('Box Plot: Fragments per Barcode')
    ax2.grid(True, alpha=0.3)

    # Add jitter points
    y_jitter = np.random.normal(1, 0.05, size=len(frag_counts_per_bc.values))
    ax2.scatter(frag_counts_per_bc.values, y_jitter, alpha=0.6, color='darkblue', s=20, edgecolors='none')

    # 3. Cumulative distribution
    sorted_counts = np.sort(frag_counts_per_bc.values)
    cumulative = np.arange(1, len(sorted_counts) + 1) / len(sorted_counts) * 100

    ax3.plot(sorted_counts, cumulative, 'b-', linewidth=2, alpha=0.8)
    ax3.set_xlabel('Number of Fragments per Barcode')
    ax3.set_ylabel('Cumulative Percentage (%)')
    ax3.set_title('Cumulative Distribution')
    ax3.grid(True, alpha=0.3)
    ax3.axhline(y=50, color='r', linestyle='--', alpha=0.7, label='Median')
    ax3.axhline(y=90, color='orange', linestyle='--', alpha=0.7, label='90th percentile')
    ax3.legend()

    # 4. Top 20 barcodes by fragment count
    top_20 = frag_counts_per_bc.head(20)
    bars = ax4.bar(range(len(top_20)), top_20.values, color='coral', alpha=0.7)
    ax4.set_xlabel('Barcode Rank')
    ax4.set_ylabel('Number of Fragments')
    ax4.set_title('Top 20 Barcodes by Fragment Count')
    ax4.set_xticks(range(len(top_20)))
    ax4.set_xticklabels(top_20.index.tolist(), rotation=45, ha='right')

    # Add value labels on bars
    for bar, count in zip(bars, top_20.values):
        height = bar.get_height()
        ax4.text(bar.get_x() + bar.get_width()/2., height + max(top_20.values)*0.01,
                f'{count}', ha='center', va='bottom', fontsize=9)

    plt.tight_layout()

    # Save plot
    plot_file = 'Calc_Frag_Length_50000/fragments_per_barcode_distribution.png'
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    plt.close()

    # Save detailed statistics
    stats_file = 'Calc_Frag_Length_50000/fragments_per_barcode_stats.txt'
    with open(stats_file, 'w') as f:
        f.write("Fragments per Barcode Distribution Statistics\n")
        f.write("=" * 50 + "\n\n")
        f.write(f"Total unique barcodes: {len(frag_counts_per_bc)}\n")
        f.write(f"Total fragments: {frag_counts_per_bc.sum()}\n")
        f.write(f"Average fragments per barcode: {frag_counts_per_bc.mean():.2f}\n")
        f.write(f"Median fragments per barcode: {frag_counts_per_bc.median():.0f}\n")
        f.write(f"Standard deviation: {frag_counts_per_bc.std():.2f}\n")
        f.write(f"Min fragments per barcode: {frag_counts_per_bc.min()}\n")
        f.write(f"Max fragments per barcode: {frag_counts_per_bc.max()}\n\n")

        f.write("Percentiles:\n")
        for p in [10, 25, 50, 75, 90, 95, 99]:
            f.write(f"{p}th percentile: {frag_counts_per_bc.quantile(p/100):.0f}\n")

        f.write(f"\nNumber of barcodes with different fragment counts:\n")
        count_dist = frag_counts_per_bc.value_counts().sort_index()
        for frag_count, barcode_count in count_dist.items():
            f.write(f"{frag_count} fragments: {barcode_count} barcodes\n")

    print(f"Distribution plot saved to: {plot_file}")
    print(f"Detailed statistics saved to: {stats_file}")
    print(f"Summary statistics:")
    print(f"  Total barcodes: {len(frag_counts_per_bc)}")
    print(f"  Average fragments per barcode: {frag_counts_per_bc.mean():.2f}")
    print(f"  Median fragments per barcode: {frag_counts_per_bc.median():.0f}")
    print(f"  Barcodes with 1 fragment: {sum(frag_counts_per_bc == 1)}")
    print(f"  Barcodes with ≥10 fragments: {sum(frag_counts_per_bc >= 10)}")

def bed_sum(bed_file_merged):
    total_len_merged = 0
    with open(bed_file_merged, 'r') as f:
        for line in f:
            info = line.strip().split('\t')
            start = int(info[1])
            end = int(info[2])
            total_len_merged += end - start + 1
    return total_len_merged

def adaptor_frequency(adaptor_seq, start, end, fq_file):
    start_idx, end_idx = start-1, end
    cnt, adaptor_cnt = 0, 0

    with gzip.open(fq_file, 'rt') as handle:
        for record in SeqIO.parse(handle, "fastq"):
            cnt +=1
            subseq = str(record.seq)[start_idx:end_idx]
            if subseq == adaptor_seq:
                adaptor_cnt+=1
    print(f'{fq_file} adaptor_cnt/total_cnt = {adaptor_cnt}/{cnt} = {round(adaptor_cnt/cnt, 2)}')
    return 

def summary_hapblock():
    # infile hapblock
    infile = 'data_hapblock'
    phased_match, span_match, span_phased_lst = [], [], []
    with open(infile, 'r') as f:
        for line in f:
            if line.startswith('BLOCK:'):
                input_string = line.strip()
                info = input_string.split(': ')
                phased_value = int(info[4].split(' ')[0])
                span_value = int(info[5].split(' ')[0])
                span_phased = round(span_value/phased_value, 2)

                phased_match.append(phased_value)
                span_match.append(span_value)
                span_phased_lst.append(span_phased)
    
    print(f'phased_count, mean={mean(phased_match)}, min={min(phased_match)}, max={max(phased_match)} std={stdev(phased_match)}')
    print(f'span_match, mean={mean(span_match)}, min={min(span_match)}, max={max(span_match)} std={stdev(span_match)}')
    print(f'span_phased, mean={mean(span_phased_lst)}, min={min(span_phased_lst)}, max={max(span_phased_lst)} std={stdev(span_phased_lst)}')

def swap_gc_in_fastq(input_fastq_path, output_fastq_path, BASE1, BASE2, BASE1_LOWER, BASE2_LOWER):
    """
    Swaps all 'G's with 'C's and 'C's with 'G's in the sequences of a FASTQ file.
    Preserves case (G -> C, g -> c).
    Handles both plain and gzipped FASTQ files.
    """
    if input_fastq_path.endswith('.gz'):
        in_opener = gzip.open
        if not output_fastq_path.endswith('.gz'):
            output_fastq_path += '.gz' # Ensure gzipped output if input is gzipped
    else:
        in_opener = open

    if output_fastq_path.endswith('.gz'):
        out_opener = gzip.open
    else:
        out_opener = open

    print(f"Processing '{input_fastq_path}'...")
    print(f"Output will be written to '{output_fastq_path}'...")

    try:
        with in_opener(input_fastq_path, 'rt') as infile, \
             out_opener(output_fastq_path, 'wt') as outfile:
            line_count = 0
            for line in infile:
                # line_count += 1
                # FASTQ format:
                # Line 1: Header (starts with @)
                # Line 2: Sequence
                # Line 3: Separator (starts with +)
                # Line 4: Quality scores

                # if line_count % 4 == 2:  # This is the sequence line
                    # Perform the swap: G -> temp, C -> G, temp -> C
                    # Using replace is tricky due to order, so we can do it char by char or use translate
                    # Char by char replacement:
                new_sequence_chars = []
                for char in line.strip():
                    if char == BASE1:
                        new_sequence_chars.append(BASE2)
                    elif char == BASE2:
                        new_sequence_chars.append(BASE1)
                    elif char == BASE2_LOWER:
                        new_sequence_chars.append(BASE1_LOWER)
                    elif char == BASE1_LOWER:
                        new_sequence_chars.append(BASE2_LOWER)
                    else:
                        new_sequence_chars.append(char)
                outfile.write("".join(new_sequence_chars) + "\n")
                    # outfile.write(line)
        print("Processing complete.")

    except FileNotFoundError:
        print(f"Error: Input file '{input_fastq_path}' not found.")
        return
    except Exception as e:
        print(f"An error occurred: {e}")
        # Clean up partially created output file if an error occurs
        if os.path.exists(output_fastq_path):
            os.remove(output_fastq_path)
        return

def base_swap(toolsdir):
    BASE1, BASE2 = 'A', 'C'
    BASE1_LOWER, BASE2_LOWER = BASE1.lower(),BASE2.lower()
    input_path = f'{toolsdir}/barcode.list'
    output_path = f'{toolsdir}/barcode.{BASE1}{BASE2}swap.list'

    if not output_path:
        base_name, ext = os.path.splitext(input_path)
        if ext.lower() == '.gz':
            base_name, gz_ext = os.path.splitext(base_name)
            output_path = f"{base_name}_GCswapped{gz_ext}.gz"
        else:
            output_path = f"{base_name}_GCswapped{ext}"

    if input_path == output_path:
        print("Error: Input and output file paths cannot be the same. Please specify a different output file.")
        return

    swap_gc_in_fastq(input_path, output_path, BASE1, BASE2, BASE1_LOWER, BASE2_LOWER)

def idxstats():
    input_file_path = "Align/idxstats_removedup_rm000.bam.txt"
    total_reads = 0
    # df = pd.read_csv(infile, index_col=False, sep='\t', header=None)
    chromosome_pattern = re.compile(r"^(chr([1-9]|1[0-9]|2[0-2])|chrX|chrY)$")

    try:
        with open(input_file_path, 'r') as f:
            for line in f:
                line = line.strip() # Remove leading/trailing whitespace
                if not line: # Skip empty lines
                    continue

                parts = line.split('\t') # Split by tab delimiter

                if len(parts) >= 3: # Ensure there are at least 3 columns
                    chromosome_name = parts[0]
                    try:
                        value_to_sum = int(parts[2]) # Convert the 3rd column to an integer
                    except ValueError:
                        # Skip lines where the 3rd column is not a valid integer
                        print(f"Warning: Skipping line due to non-integer value in 3rd column: '{line}'", file=sys.stderr)
                        continue

                    # Check if the chromosome name matches the criteria
                    if chromosome_pattern.match(chromosome_name):
                        total_reads += value_to_sum
                else:
                    print(f"Warning: Skipping line with fewer than 3 columns: '{line}'", file=sys.stderr)

    except FileNotFoundError:
        print(f"Error: Input file not found at '{input_file_path}'", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"An unexpected error occurred: {e}", file=sys.stderr)
        sys.exit(1)

    hg38_len=2934876545
    df = pd.read_csv('Align/picard_align_metrics_BCreads_only.txt', index_col=False, sep='\t')
    tmp = df[df['Sample'] == 'Uniq bases(bp)']
    uniq_dedup_withBC_bases = int(tmp['data'].item())
    with open("Align/mapped_uniq_bc_bases_count.txt", 'w') as f:
        f.write(f'mapped_primary_withBC_reads {total_reads}\n')
        f.write(f'avg_mapped_primary_withBC_bases_per_read {round(uniq_dedup_withBC_bases/(total_reads), 2)}\n')
    return  

def recreate_split_log(fastq_file):
    # input data/fq2
    umi_counts = Counter()
    Reads_pair_num = 0
    # 正则表达式用于从 BX:Z: 标签中提取 UMI
    # 匹配 BX:Z: 后面的非空格字符序列，直到遇到空格或行尾
    umi_pattern = re.compile(r'BX:Z:(\S+)')
    
    # if fastq_file.endswith('.gz'):
    #     # 使用 gzip.open 处理 .gz 文件
    #     # 'rt' 模式表示以文本模式读取 GZIP 文件
    #     _open = gzip.open
    #     open_mode = 'rt'
    # else:
    #     # 对于普通文本文件，使用内置的 open
    #     _open = open
    #     open_mode = 'r'

    # print(f"开始处理 FASTQ 文件: {fastq_file}...")
    try:
        with gzip.open(fastq_file, 'rt') as f:
            line_num = 0
            for line in f:
                line_num += 1
                if line_num % 4 == 1:  # FASTQ 文件的第一行是头信息
                    # 查找 BX:Z: 标签
                    match = umi_pattern.search(line)
                    if match:
                        umi = match.group(1)
                        umi_counts[umi] += 1
                    else:
                        # 如果没有 BX:Z: 标签，尝试从 # 后面提取 UMI
                        # 假设 # 后面的第一个下划线分隔的字段是 UMI
                        # 例如: E250068551_L01_C001R001_0#0_0_0/2
                        parts = line.split('#')
                        if len(parts) > 1:
                            sub_parts = parts[1].split(' ')[0].split('/')
                            if len(sub_parts) > 0:
                                # 假设 UMI 是 # 后面的第一个部分，直到 /
                                umi_from_id = sub_parts[0]
                                if umi_from_id: # 确保不是空字符串
                                    umi_counts[umi_from_id] += 1
                                else:
                                    print(f"警告: 第 {line_num} 行未能提取 UMI (ID部分): {line.strip()}")
                            else:
                                print(f"警告: 第 {line_num} 行未能提取 UMI (ID部分): {line.strip()}")
                        else:
                            print(f"警告: 第 {line_num} 行未能提取 UMI (无BX:Z:也无#): {line.strip()}")

                # if line_num % 100000 == 0:
                #     print(f"  已处理 {line_num} 行...")

        bc_types = len(umi_counts)
        num_bc000 = umi_counts['0_0_0']
        Reads_pair_num = line_num
        Reads_pair_num_split = Reads_pair_num-num_bc000
        # print(f"FASTQ 文件处理完成。共找到 {len(umi_counts)} 种不同的 UMI。")

        # 将结果写入文件
        with open('split.log', 'w') as out_f:
            out_f.write('Barcode_types = 1536 * 1536 * 1536 = 3623878656 \n')
            out_f.write(f'Real_Barcode_types = {bc_types} ({bc_types/3623878656} %)\n')
            out_f.write(f'Reads_pair_num  = {Reads_pair_num}\n')
            out_f.write(f'Reads_pair_num(after split) = {Reads_pair_num_split} ({Reads_pair_num_split/Reads_pair_num} %)\n')
        
            idx = 1 # 手动维护索引
            for umi_type, count in umi_counts.items():
                out_f.write(f"{idx}\t{count}\t{umi_type}\n")
                idx += 1


        # print(f"UMI 计数结果已保存到: {output_file}")

    except FileNotFoundError:
        print(f"错误：输入文件 '{fastq_file}' 未找到。")
    except Exception as e:
        print(f"处理文件时发生错误：{e}")


if __name__ == "__main__":

    # module ='adaptor_frequency'
    parser = argparse.ArgumentParser()
    parser.add_argument("--module", type=str)
    parser.add_argument("--fq_file", type=str, required=False)
    parser.add_argument("--minreads", type=int, required=False)
    parser.add_argument("--minreads_fasta", type=int, required=False)
    parser.add_argument("--batch", type=str, required=False)
    parser.add_argument("--dir", type=str, required=False, help="Base results directory")
    parser.add_argument("--toolsdir", type=str, required=False, help="Path to Data_and_Tools directory")

    args = parser.parse_args()
    module = args.module
    fq_file = args.fq_file
    minreads = args.minreads
    minreads_fasta = args.minreads_fasta
    batch = args.batch
    # batch='V350134870'

    bc_list = ['L01', 'L02', 'L03', 'L04']
    _dir = f'{args.dir}/{batch}' if args.dir else f'./{batch}'

    if module =='combine_split_log_summary':
        # count_per_bc = 4
        ins = combine_demulx()
        log_list = [f'{_dir}/{bc}/analysis/split_stat_read1.log' for bc in bc_list]
        ins.combine_split_log(log_list, minreads_fasta, module)
    
    elif module == 'recreate_split_log':
        fastq_file = "data/split_read.2.fq.gz"
        recreate_split_log(fastq_file)

    elif module =='idxstats':
        idxstats()

    elif module == 'base_swap':
        base_swap(args.toolsdir)

    elif module == 'combine_split_log_exon2fasta':
        ins = combine_demulx()
        log_list = [f'{_dir}/{bc}/analysis/split_stat_read1.log' for bc in bc_list]
        ins.combine_split_log(log_list, minreads_fasta, module )

    elif module =='across_sample':
        lanes = [43, 57, 58, 59, 60, 61, 63, 64]
        # lanes = [43, 57]
        across_sample(lanes)
    
    elif module == 'calculate_expected_beads':
        num_sample = 8  
        per_sample_frag = 20_000_000  
        num_beads = 30_000_000  
        calculate_expected_beads(num_beads, per_sample_frag, num_sample, )

    elif module =='frag_bed2summary':
        frag_bed2summary(minreads)

    elif module =='adaptor_frequency':
        start, end = 6, 27
        adaptor_seq = 'GAGACGTTCTCGACTCAGCAGA'
        # fq_file ='test.fq.gz'
        adaptor_frequency(adaptor_seq, start, end, fq_file)

    elif module == 'summary_hapblock':
        summary_hapblock()

    else:
        print('need a module name')