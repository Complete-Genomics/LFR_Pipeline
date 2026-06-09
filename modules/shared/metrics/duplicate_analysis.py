#!/usr/bin/env python3
"""Barcode-aware duplicate analysis.

Python replacement for the legacy barcode duplicate helper.  It reads a BAM
directly and writes the same three output files expected by the workflow.
"""

import argparse
import os
from collections import defaultdict

import pysam


def parse_read_identity(query_name):
    parts = query_name.split("#", 1)
    if len(parts) == 1:
        return query_name, "NA"
    return parts[0], parts[1]


def analyze_duplicates(bam_path, outdir):
    os.makedirs(outdir, exist_ok=True)
    locate_hash = defaultdict(lambda: defaultdict(dict))
    dup_pe = 0

    pe_dup_reads_path = os.path.join(outdir, "PE_dup_reads")
    with pysam.AlignmentFile(bam_path, "rb") as bam, open(pe_dup_reads_path, "w") as pe_out:
        for read in bam.fetch(until_eof=True):
            if not read.is_duplicate:
                continue

            read_id, barcode = parse_read_identity(read.query_name)
            chrom = read.reference_name or "*"
            pos = read.reference_start + 1 if read.reference_start >= 0 else 0
            chr_pos = f"{chrom}:{pos}"
            seq = read.query_sequence or ""

            if read_id not in locate_hash[chr_pos][barcode]:
                locate_hash[chr_pos][barcode][read_id] = seq
            else:
                dup_pe += 1
                pe_out.write(
                    f"{chr_pos}\t{barcode}\t{read_id}\t"
                    f"1:{locate_hash[chr_pos][barcode][read_id]}\t2:{seq}\n"
                )

    dup_info_path = os.path.join(outdir, "dup_info")
    duplicate_rate_path = os.path.join(outdir, "duplicate_rate")

    bar_total = 0
    reads_total = 0
    dup_reads = 0
    un_snp_total = 0
    bar_reads_total = 0

    with open(dup_info_path, "w") as info_out:
        info_out.write("#duplicate location\tbar_in_reads_rate\tdupbar_in_reads_rate\tun_snp_in_bar_rate\n")
        for chr_pos, barcode_map in locate_hash.items():
            bar_num = len(barcode_map)
            total_reads = 0
            false_dup = 0
            un_snp_num = 0

            for read_map in barcode_map.values():
                reads_num = len(read_map)
                total_reads += reads_num
                if reads_num > 1:
                    dup_reads += reads_num
                    seq_counts = defaultdict(int)
                    for seq in read_map.values():
                        seq_counts[seq] += 1
                    un_snp_num += sum(count for count in seq_counts.values() if count > 1)
                else:
                    false_dup += 1

            if total_reads == 0:
                continue

            bar_reads = total_reads - false_dup
            bar_in_reads_rate = bar_num / total_reads
            dupbar_in_reads_rate = 1 - (false_dup / total_reads)
            un_snp_in_bar_rate = (un_snp_num / bar_reads) if bar_reads else un_snp_num

            un_snp_total += un_snp_num
            bar_total += bar_num
            reads_total += total_reads
            bar_reads_total += bar_reads

            info_out.write(
                f"{chr_pos}\t{bar_in_reads_rate}\t"
                f"{dupbar_in_reads_rate}\t{un_snp_in_bar_rate}\n"
            )

    dup_bar_rate = (bar_total / reads_total) if reads_total else 0
    with open(duplicate_rate_path, "w") as rate_out:
        rate_out.write(f"duplicate reads number:{reads_total}\n")
        rate_out.write(f"duplicate barcode number:{bar_total}\n")
        rate_out.write(f"bar_rate in dup reads:{dup_bar_rate}\n")


def main():
    parser = argparse.ArgumentParser(description="Barcode-aware duplicate analysis.")
    parser.add_argument("bam")
    parser.add_argument("outdir")
    args = parser.parse_args()
    analyze_duplicates(args.bam, args.outdir)


if __name__ == "__main__":
    main()
