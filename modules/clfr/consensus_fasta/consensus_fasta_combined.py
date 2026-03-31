'''
usage:
python $src --bam $bam --ref_fasta $ref --output_fasta output.fa --start_index 2000 --end_index 4000
'''

import subprocess
import os
import tempfile
import pysam
import sys
import shutil
import glob
import itertools
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--bam", type=str, required=False)
parser.add_argument("--ref_fasta", type=str, required=False)
parser.add_argument("--output_fasta", type=str, required=False)
parser.add_argument("--chrom", type=str, required=False)
# parser.add_argument("--start_index", type=int, required=False)
parser.add_argument("--dict_file", type=str, required=False)
parser.add_argument("--split_index", type=int, required=False)
parser.add_argument("--min_reads", type=int, required=False)
parser.add_argument("--samtools", type=str, required=False, default="samtools")
parser.add_argument("--stringtie", type=str, required=False, default="stringtie")

args = parser.parse_args()
input_bam_file = args.bam
reference_fasta_file = args.ref_fasta
output_fasta_file = args.output_fasta

chrom = args.chrom
split_index = args.split_index
dict_file = args.dict_file
MIN_READS = args.min_reads

def get_idx(dict_file, split_index, chrom):
    chr_dict = {}
    with open(dict_file, 'r') as f:
        for line in f:
            fields = line.split('\t')
            chr_name = fields[0]
            if chr_name == chrom:
                mapped_reads = int(fields[2])
                chunk = mapped_reads//5
                start_index = split_index*chunk+1
                end_index = (split_index+1)*chunk

    return start_index, end_index

start_index, end_index  = get_idx(dict_file, split_index, chrom)


# --- Configuration ---
SAMTOOLS_PATH = args.samtools
STRINGTIE_PATH = args.stringtie
MIN_FRAG_LEN = 400
# 使用一个唯一的临时目录，确保不会与其他进程冲突
TEMP_BASE_DIR = "/dev/shm"
TEMP_DIR_NAME = f"consensus_single_thread_tmp_{os.getpid()}"
FALLBACK_TEMP_DIR = "/tmp"

REF = reference_fasta_file
subprocess.call(f'mkdir -p /dev/shm/consensus', shell=True)
TEMP_DIR = f"/dev/shm/consensus/{chrom}_{split_index}"
subprocess.call(f'mkdir -p {TEMP_DIR}', shell=True)
current_temp_dir = TEMP_DIR

# --- Global variable for FASTA reference ---
_fasta_ref_file = None

# def get_actual_temp_dir():
#     """Determine the actual temporary directory to use, prefer /dev/shm/, fallback to /tmp."""
#     candidate_temp_dir = os.path.join(TEMP_BASE_DIR, TEMP_DIR_NAME)
    
#     # 尝试在 /dev/shm 创建目录
#     try:
#         os.makedirs(candidate_temp_dir, exist_ok=True)
#         if os.path.exists(candidate_temp_dir) and os.access(candidate_temp_dir, os.W_OK):
#             stat = shutil.disk_usage(candidate_temp_dir)
#             # 至少保留 500MB 空间，根据实际 UMI 组数量和大小调整
#             if stat.free > 1024 * 1024 * 500:
#                 sys.stderr.write(f"Using {candidate_temp_dir} for temporary files.\n")
#                 return candidate_temp_dir
#             else:
#                 sys.stderr.write(f"WARNING: {candidate_temp_dir} space is low ({stat.free / (1024*1024):.2f} MB available). Falling back to {FALLBACK_TEMP_DIR}.\n")
#         else:
#             sys.stderr.write(f"WARNING: {candidate_temp_dir} not accessible or writable. Falling back to {FALLBACK_TEMP_DIR}.\n")
#     except OSError as e:
#         sys.stderr.write(f"WARNING: Could not create directory {candidate_temp_dir} ({e}). Falling back to {FALLBACK_TEMP_DIR}.\n")
    
#     # 如果 /dev/shm 不可用或空间不足，使用 /tmp
#     fallback_dir = os.path.join(FALLBACK_TEMP_DIR, TEMP_DIR_NAME)
#     os.makedirs(fallback_dir, exist_ok=True)
#     sys.stderr.write(f"Using {fallback_dir} for temporary files (fallback).\n")
#     return fallback_dir

def initialize_fasta_ref(ref_path):
    """Initializes the pysam FastaFile object once."""
    global _fasta_ref_file
    if _fasta_ref_file is None:
        if not os.path.exists(f"{ref_path}.fai"):
            sys.stderr.write(f"WARNING: Reference FASTA file {ref_path} does not have a .fai index. Creating one now...\n")
            try:
                subprocess.run([SAMTOOLS_PATH, "faidx", ref_path], check=True, capture_output=True)
            except subprocess.CalledProcessError as e:
                sys.stderr.write(f"ERROR: Failed to create FASTA index for {ref_path}: {e.stderr.decode()}\n")
                raise
        _fasta_ref_file = pysam.FastaFile(ref_path)

def get_umi_from_read_id(read_id):
    """
    Extract UMI from read ID (format: Vxxx#UMI_id)
    """
    parts = read_id.split('#')
    if len(parts) > 1:
        return parts[-1]
    return None

def generate_bed_data(reads_list_obj, header_dict, stringtie_path, current_temp_dir):
    """Generate BED data in memory from BAM objects using stringtie."""
    
    # Create temp BAM for StringTie
    with tempfile.NamedTemporaryFile(suffix=".bam", dir=current_temp_dir, delete=False) as temp_bam:
        header = pysam.AlignmentHeader.from_dict(header_dict)
        with pysam.AlignmentFile(temp_bam.name, "wb", header=header) as bam_out:
            for read_obj in reads_list_obj:
                try:
                    bam_out.write(read_obj)
                except Exception as e:
                    sys.stderr.write(f"Error writing read to temp BAM for StringTie: {e}\n")
                    continue
        temp_bam_path = temp_bam.name

    if os.path.getsize(temp_bam_path) == 0:
        sys.stderr.write(f"ERROR: Temporary BAM file {temp_bam_path} is empty. Skipping StringTie.\n")
        os.unlink(temp_bam_path)
        return []

    with tempfile.NamedTemporaryFile(suffix=".gtf", dir=current_temp_dir, delete=False) as temp_gtf:
        temp_gtf_path = temp_gtf.name
    # Run stringtie (no indexing, as sorting is prohibited)
    # env = os.environ.copy()
    # env["TMPDIR"] = current_temp_dir # Set TMPDIR for StringTie if needed
    stringtie_cmd = [
        stringtie_path,
        temp_bam_path,
        "-o", temp_gtf_path, # Output to stdout
        "-p", "1"  # Use 1 thread for stringtie
    ]
    try:
        result = subprocess.run(stringtie_cmd, capture_output=True, text=True, check=True)
        # gtf_lines = result.stdout.strip().split("\n")
        with open(temp_gtf_path, 'r') as gtf_file:
            gtf_lines = gtf_file.read().strip().split("\n")
        
        # this concatenate all isoforms 
        # exon_lines = [line for line in gtf_lines if "exon" in line]
        # only use 1st isoforms, 
        transcript_id = None
        for line in gtf_lines:
            if "transcript" in line and "transcript_id" in line:
                # Extract the first transcript ID
                transcript_id = line.split('transcript_id "')[-1].split('"')[0]
                break

        # Filter for exons that belong to the first transcript
        exon_lines = []
        if transcript_id:
            for line in gtf_lines:
                if "exon" in line and f'transcript_id "{transcript_id}"' in line:
                    exon_lines.append(line)
        else:
            # Fallback to original behavior if no transcript is found
            exon_lines = [line for line in gtf_lines if "exon" in line]

        bed_lines = []
        for line in exon_lines:
            # if "exon" in exon_lines: # Only interested in exons
            fields = line.strip().split("\t")
            if len(fields) < 8:
                sys.stderr.write(f"WARNING: GTF line format invalid: {line}\n")
                continue
            try:
                chrom = fields[0]
                # GTF is 1-based start, 1-based end. BED is 0-based start, 1-based end.
                start = int(fields[3]) # Convert to 0-based
                end = int(fields[4])       # Already 1-based end
                length = end - start
                bed_lines.append(f"{chrom}\t{start}\t{end}\t{length}")
            except (ValueError, IndexError) as e:
                sys.stderr.write(f"Error parsing GTF line: {line}, Error: {e}\n")
                continue
    except subprocess.CalledProcessError as e:
        sys.stderr.write(f"Error running stringtie: {e.stderr}\n")
        bed_lines = []
    finally:
        if os.path.exists(temp_bam_path):
            os.unlink(temp_bam_path)

    return bed_lines

def get_2bp_sequence_pysam(chrom, start, end):
    """Extract 2bp sequences from start and end of a region using pysam.FastaFile."""
    if _fasta_ref_file is None:
        raise RuntimeError("Reference FASTA file not initialized for pysam.")
    
    # Ensure chromosome naming matches the FASTA reference.
    # Assuming the FASTA reference has 'chr' prefix (e.g., 'chr22').
    # If not, you might need to convert chrom here or modify the FASTA index.
    if not chrom.startswith('chr'):
        chrom = 'chr' + chrom

    try:
        # Pysam.fetch(contig, start, end) -> 0-based start, 0-based end-exclusive
        # Extracting the 2bp sequences based on the original bedtools logic
        seq1 = _fasta_ref_file.fetch(chrom, start - 3, start - 1)
        seq2 = _fasta_ref_file.fetch(chrom, end, end + 2)
        # print(f'{start},{end}')
        return seq1, seq2
    except Exception as e:
        sys.stderr.write(f"WARNING: Error extracting 2bp sequence from {chrom}:{start}-{end} using pysam: {e}\n")
        return "", ""

def process_umi_group_single_thread(umi_id, reads_list_obj, header_dict_data, ref_fasta_path, current_temp_dir):
    """Process a UMI group's reads to generate a consensus sequence."""
    
    if not reads_list_obj:
        sys.stderr.write(f"WARNING: UMI ID '{umi_id}' has no reads data, skipping.\n")
        return None

    if len(reads_list_obj) < MIN_READS:
        # sys.stderr.write(f"WARNING: UMI ID '{umi_id}' read count {len(reads_list_obj)} < {MIN_READS}, skipping.\n")
        return None

    bed_lines = generate_bed_data(reads_list_obj, header_dict_data, STRINGTIE_PATH, current_temp_dir)
    if not bed_lines:
        # sys.stderr.write(f"WARNING: UMI ID '{umi_id}' has no valid BED regions, skipping.\n")
        return None

    # Create temp BAM for Samtools consensus
    temp_bam_path = None # Initialize to None for finally block
    try:
        with tempfile.NamedTemporaryFile(suffix=".bam", dir=current_temp_dir, delete=False) as temp_bam:
            header = pysam.AlignmentHeader.from_dict(header_dict_data)
            with pysam.AlignmentFile(temp_bam.name, "wb", header=header) as bam_out:
                for read_obj in reads_list_obj:
                    try:
                        # read = pysam.AlignedSegment.fromstring(read_obj.decode('utf-8'), header=header)
                        bam_out.write(read_obj)
                    except Exception as e:
                        sys.stderr.write(f"Error writing read to temp BAM for UMI {umi_id}: {e}\n")
                        continue
            temp_bam_path = temp_bam.name

        if os.path.getsize(temp_bam_path) == 0:
            sys.stderr.write(f"ERROR: Temporary BAM file for UMI {umi_id} ({temp_bam_path}) is empty. Skipping consensus.\n")
            return None
        else:
            # Index the temp BAM file
            try:
                subprocess.run([SAMTOOLS_PATH, "index", temp_bam_path], check=True, capture_output=True)
            except subprocess.CalledProcessError as e:
                sys.stderr.write(f"ERROR: Failed to index temporary BAM {temp_bam_path} for UMI {umi_id}: {e.stderr.decode()}\n")
                return None


        combined_seq = ''
        for line in bed_lines:
            if line.strip() and not line.startswith("#"):
                try:
                    fields = line.strip().split()
                    if len(fields) < 4:
                        sys.stderr.write(f"WARNING: BED line format invalid (UMI {umi_id}): {line}\n")
                        continue
                    chrom, start_0based, end_1based, length = fields[:4]
                    
                    # Samtools consensus expects 1-based start and end for region
                    # And usually expects 'chr' prefix if your FASTA has it
                    consensus_start = int(start_0based) + 2 # +1 for 1-based, +2 for clipping
                    consensus_end = int(end_1based) - 2 # -2 for clipping
                    
                    # Determine correct chromosome name for samtools consensus
                    samtools_chrom = chrom
                    if not samtools_chrom.startswith('chr'):
                        samtools_chrom = 'chr' + samtools_chrom

                except (ValueError, IndexError) as e:
                    sys.stderr.write(f"Error parsing BED line (UMI {umi_id}): {line}, Error: {e}\n")
                    continue
                
                # Skip if region becomes invalid after clipping
                if consensus_end <= consensus_start:
                    sys.stderr.write(f"WARNING: Clipped region for UMI {umi_id} is invalid: {samtools_chrom}:{consensus_start}-{consensus_end}. Skipping.\n")
                    continue
                # -aa, Output absolutely all positions, including references with no data aligned against them
                region = f"{samtools_chrom}:{consensus_start}-{consensus_end}"
                cmd = [
                    SAMTOOLS_PATH,
                    "consensus",
                    "-aa",
                    "-r", region,
                    "--show-ins", "yes",
                    temp_bam_path
                ]
                try:
                    result = subprocess.run(
                        cmd,
                        capture_output=True,
                        text=True,
                        check=True
                    )
                    consensus_fasta = result.stdout.strip()
                    consensus_lines = consensus_fasta.split("\n")
                    if len(consensus_lines) < 2:
                        sys.stderr.write(f"WARNING: Samtools consensus output for UMI {umi_id} is empty or invalid (region {region}).\n")
                        continue
                    consensus_seq_core = "".join(consensus_lines[1:])
                    
                    # Extract 2bp sequences from reference FASTA using pysam
                    start_2bp, end_2bp = get_2bp_sequence_pysam(chrom, consensus_start, consensus_end)
                    
                    combined_seq += start_2bp + consensus_seq_core + end_2bp
                except subprocess.CalledProcessError as e:
                    sys.stderr.write(f"Error running samtools consensus (UMI {umi_id}, region {region}): {e.stderr}\n")
                    continue
                except Exception as e:
                    sys.stderr.write(f"Unexpected error processing UMI {umi_id}, region {region}: {e}\n")
                    continue
    finally:
        # Cleanup temporary BAM and its index
        if temp_bam_path and os.path.exists(temp_bam_path):
            os.unlink(temp_bam_path)
        if temp_bam_path and os.path.exists(temp_bam_path + ".bai"): # Samtools index creates a .bai file
            os.unlink(temp_bam_path + ".bai")

    if len(combined_seq) < MIN_FRAG_LEN:
        # sys.stderr.write(f"WARNING: UMI ID '{umi_id}' combined sequence {len(combined_seq)}bp < {MIN_FRAG_LEN}bp, skipping.\n")
        return None

    return f">{umi_id}_{chrom}\n{combined_seq}\n"

def generate_consensus_sequential_to_single_file(input_bam, reference_fasta, output_fasta_file, start_index, end_index):
    """Generate consensus sequences for each UMI from a read ID-sorted BAM file in a single thread."""
    if not os.path.exists(input_bam):
        sys.stderr.write(f"ERROR: Input BAM file '{input_bam}' does not exist.\n")
        sys.exit(1)
    if not os.path.exists(reference_fasta):
        sys.stderr.write(f"ERROR: Reference FASTA file '{reference_fasta}' does not exist.\n")
        sys.exit(1)

    # Initialize global FASTA reference (once)
    initialize_fasta_ref(reference_fasta)
    
    # Get the unique temporary directory for this run
    current_temp_dir = TEMP_DIR

    try:
        infile = pysam.AlignmentFile(input_bam, "rb")
    except pysam.SamtoolsError as e:
        sys.stderr.write(f"ERROR: Could not open BAM file '{input_bam}': {e}\n")
        sys.exit(1)

    output_dir = os.path.dirname(output_fasta_file)
    os.makedirs(output_dir or '.', exist_ok=True)

    header_dict = infile.header.to_dict()
    current_umi = None
    current_reads_buffer = [] # Store pysam.AlignedSegment objects directly
    read_count = 0
    umi_count = 0

    sys.stderr.write(f"Reading BAM file '{input_bam}' and processing UMI groups (single-threaded)...\n")

    try:
        with open(output_fasta_file, 'wt') as output_file_handle:
            # for read in infile:
            for i, read in enumerate(itertools.islice(infile, start_index, end_index)):
                read_count += 1
                umi_id = get_umi_from_read_id(read.query_name)
                if umi_id is None:
                    sys.stderr.write(f"WARNING: Skipping read ID '{read.query_name}' due to invalid UMI format.\n")
                    continue

                if current_umi is None:
                    current_umi = umi_id
                
                if umi_id != current_umi:
                    if current_reads_buffer:
                        fasta_record = process_umi_group_single_thread(
                            current_umi, current_reads_buffer, header_dict, reference_fasta, current_temp_dir
                        )
                        if fasta_record:
                            output_file_handle.write(fasta_record)
                            umi_count += 1
                            # if umi_count % 1000 == 0: # Print progress more frequently for single thread
                                # sys.stderr.write(f"Processed {umi_count} UMI groups so far...\n")
                    current_umi = umi_id
                    current_reads_buffer = [read]
                else:
                    current_reads_buffer.append(read)

            # Process the last UMI group
            if current_umi is not None and current_reads_buffer:
                fasta_record = process_umi_group_single_thread(
                    current_umi, current_reads_buffer, header_dict, reference_fasta, current_temp_dir
                )
                if fasta_record:
                    output_file_handle.write(fasta_record)
                    umi_count += 1
                    sys.stderr.write(f"Processed {umi_count} UMI groups (final count).\n")

    except IOError as e:
        sys.stderr.write(f"ERROR: Could not write to FASTA file '{output_fasta_file}': {e}\n")
        sys.exit(1)
    except Exception as e:
        sys.stderr.write(f"ERROR: An unexpected error occurred: {e}\n")
        sys.exit(1)
    finally:
        try:
            infile.close()
        except Exception as e:
            sys.stderr.write(f"Error closing input BAM file: {e}\n")

    sys.stderr.write(f"Processing complete. All consensus sequences written to '{output_fasta_file}'.\n")

    # Clean up the specific temporary directory created by this run
    sys.stderr.write(f"Cleaning up temporary directory: {current_temp_dir}...\n")
    if os.path.exists(current_temp_dir):
        shutil.rmtree(current_temp_dir, ignore_errors=True)
    sys.stderr.write("Temporary directory cleaned.\n")


if __name__ == "__main__":
    # if len(sys.argv) != 4:
    #     sys.stderr.write("用法: python script.py <输入 readID 排序的 BAM> <参考 FASTA> <输出 FASTA>\n")
    #     sys.exit(1)

    # input_bam_file = sys.argv[1]
    # reference_fasta_file = sys.argv[2]
    # output_fasta_file = sys.argv[3]



    generate_consensus_sequential_to_single_file(input_bam_file, reference_fasta_file, output_fasta_file, start_index, end_index)