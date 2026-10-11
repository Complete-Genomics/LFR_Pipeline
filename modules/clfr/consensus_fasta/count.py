#!/usr/bin/env python3
"""Count successful consensus UMI groups and their reads.

The consensus workflow writes headers as ``>{umi}_{read_count}_{chrom}``. This
script sums ``read_count`` from those headers and reads the total mapped-read
count from ``samtools idxstats`` output. It does not assemble sequences or
create a FASTA file.

usage:
<python> modules/clfr/consensus_fasta/count.py \
  --consensus-fasta consensus/consensus.fasta \
  --idxstats Align/samtools_idx.txt
"""

import argparse
import sys

def consensus_counts(consensus_fasta):
    assembled_reads = 0
    assembled_umi_groups = 0
    with open(consensus_fasta) as handle:
        for line in handle:
            if not line.startswith(">"):
                continue
            record_id = line[1:].strip().split()[0]
            parts = record_id.rsplit("_", 2)
            if len(parts) == 3 and parts[1].isdigit():
                assembled_reads += int(parts[1])
                assembled_umi_groups += 1
                continue
            raise ValueError(
                "Consensus FASTA header must use the '<umi>_<read_count>_<chrom>' format: "
                f"{record_id}"
            )
    return assembled_reads, assembled_umi_groups


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


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--consensus-fasta", required=True)
    parser.add_argument("--idxstats", required=True,
                        help="samtools idxstats output for the mapped-read total.")
    args = parser.parse_args()

    mapped_reads = mapped_reads_from_idxstats(args.idxstats)
    assembled_reads, assembled_umi_groups = consensus_counts(args.consensus_fasta)

    print(f"mapped_reads\t{mapped_reads}")
    print(f"consensus_assembled_reads\t{assembled_reads}")
    print(f"consensus_assembled_umis\t{assembled_umi_groups}")


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError) as error:
        sys.exit(f"ERROR: {error}")
