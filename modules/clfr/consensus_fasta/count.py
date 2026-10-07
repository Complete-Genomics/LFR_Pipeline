#!/usr/bin/env python3
"""Count successful consensus UMI groups and their mapped reads.

The consensus workflow writes headers as ``>{umi}_{read_count}_{chrom}``. This
script uses those headers as the successful consensus groups, then counts BAM reads whose
UMI and reference chromosome match one of those groups. The total mapped-read
count is read from ``samtools idxstats`` output. It does not assemble sequences
or create a FASTA file. Chromosome names with and without a ``chr`` prefix are
treated as equivalent.

usage:
<python> modules/clfr/consensus_fasta/count.py \
  --consensus-fasta consensus/consensus.fasta \
  --bam Align/<sample>.sort.removedup_rm000.bam \
  --idxstats Align/samtools_idx.txt
"""

import argparse
import sys

import pysam


def normalized_chromosome(chrom):
    return chrom[3:] if chrom.startswith("chr") else chrom


def successful_groups(consensus_fasta):
    groups = set()
    with open(consensus_fasta) as handle:
        for line in handle:
            if not line.startswith(">"):
                continue
            record_id = line[1:].strip().split()[0]
            parts = record_id.rsplit("_", 2)
            if len(parts) == 3 and parts[1].isdigit():
                umi, _, chrom = parts
            elif len(parts) == 2:
                umi, chrom = parts
            else:
                raise ValueError(
                    "Consensus FASTA header must use the '<umi>_<read_count>_<chrom>' format: "
                    f"{record_id}"
                )
            groups.add((umi, normalized_chromosome(chrom)))
    return groups


def umi_from_read_name(read_name):
    if "#" not in read_name:
        return None
    umi = read_name.rsplit("#", 1)[1]
    if umi.endswith(("/1", "/2")):
        umi = umi[:-2]
    return umi or None


def mapped_reads_from_idxstats(idxstats_path):
    mapped_reads = 0
    with open(idxstats_path) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 4:
                raise ValueError(f"Invalid samtools idxstats line: {line.rstrip()}")
            if fields[0] != "*":
                mapped_reads += int(fields[2])
    return mapped_reads


def count_assembled_reads(bam_path, groups):
    assembled_reads = 0

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam.fetch(until_eof=True):
            if read.is_unmapped:
                continue
            umi = umi_from_read_name(read.query_name)
            if umi is not None and (umi, normalized_chromosome(read.reference_name)) in groups:
                assembled_reads += 1

    return assembled_reads


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--consensus-fasta", required=True)
    parser.add_argument("--bam", required=True,
                        help="BAM used to generate the consensus FASTA.")
    parser.add_argument("--idxstats", required=True,
                        help="samtools idxstats output for the consensus source BAM.")
    args = parser.parse_args()

    groups = successful_groups(args.consensus_fasta)
    mapped_reads = mapped_reads_from_idxstats(args.idxstats)
    assembled_reads = count_assembled_reads(args.bam, groups)

    print(f"mapped_reads\t{mapped_reads}")
    print(f"consensus_assembled_reads\t{assembled_reads}")
    print(f"consensus_assembled_umis\t{len(groups)}")


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, pysam.SamtoolsError) as error:
        sys.exit(f"ERROR: {error}")
