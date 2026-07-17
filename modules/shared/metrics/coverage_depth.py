#!/usr/bin/env python3
"""Coverage/depth summary with duplicate reads removed."""

import argparse
import os
import subprocess
from collections import Counter


def count_non_n_bases(fasta_path):
    total = 0
    with open(fasta_path) as handle:
        for line in handle:
            if line.startswith(">"):
                continue
            seq = line.strip().upper()
            total += sum(1 for base in seq if base not in {"N", "n"})
    return total


def depth_histogram(bam_path, min_mapq, min_baseq, samtools):
    hist = Counter()
    exclude_flags = str(256 + 512 + 1024 + 2048)
    command = [
        samtools,
        "depth",
        "-q", str(min_baseq),
        "-Q", str(min_mapq),
        "-G", exclude_flags,
        bam_path,
    ]
    proc = subprocess.Popen(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    assert proc.stdout is not None
    for line in proc.stdout:
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 3:
            continue
        depth = int(fields[2])
        if depth > 0:
            hist[depth] += 1

    _, stderr = proc.communicate()
    if proc.returncode != 0:
        raise SystemExit(
            "samtools depth failed with exit code {}:\n{}".format(
                proc.returncode, stderr
            )
        )
    return hist


def write_tables(hist, total_bases, outdir):
    os.makedirs(outdir, exist_ok=True)
    depths = sorted(hist)

    with open(os.path.join(outdir, "depth_frequency.txt"), "w") as out:
        for depth in depths:
            out.write(f"{depth}\t{hist[depth] / total_bases}\n")

    suffix_count = 0
    cumulative = {}
    for depth in reversed(depths):
        suffix_count += hist[depth]
        cumulative[depth] = suffix_count / total_bases

    with open(os.path.join(outdir, "cumu.txt"), "w") as out:
        for depth in depths:
            out.write(f"{depth}\t{cumulative[depth]}\n")


def summarize(hist, total_bases):
    total_covered_depth = sum(depth * count for depth, count in hist.items())
    covered_bases = sum(hist.values())

    def threshold(min_depth):
        depth_sum = sum(depth * count for depth, count in hist.items() if depth >= min_depth)
        base_sum = sum(count for depth, count in hist.items() if depth >= min_depth)
        return depth_sum / total_bases, base_sum / total_bases

    avg_depth = total_covered_depth / total_bases
    coverage = covered_bases / total_bases
    _, cov4 = threshold(4)
    _, cov10 = threshold(10)
    _, cov20 = threshold(20)

    return avg_depth, coverage, cov4, cov10, cov20


def main():
    parser = argparse.ArgumentParser(description="Coverage/depth summary excluding duplicate reads.")
    parser.add_argument("--bam", required=True)
    parser.add_argument("--ref", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--mapq", type=int, default=0)
    parser.add_argument("--baseq", type=int, default=0)
    parser.add_argument("--samtools", default="samtools")
    args = parser.parse_args()

    total_bases = count_non_n_bases(args.ref)
    hist = depth_histogram(args.bam, args.mapq, args.baseq, args.samtools)
    write_tables(hist, total_bases, args.outdir)
    avg_depth, coverage, cov4, cov10, cov20 = summarize(hist, total_bases)

    print(f"Average sequencing depth\t{avg_depth:.2f}")
    print(f"Coverage\t{100 * coverage:.2f}%")
    print(f"Coverage at least 4X\t{100 * cov4:.2f}%")
    print(f"Coverage at least 10X\t{100 * cov10:.2f}%")
    print(f"Coverage at least 20X\t{100 * cov20:.2f}%")


if __name__ == "__main__":
    main()
