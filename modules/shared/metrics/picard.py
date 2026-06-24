#!/usr/bin/env python
"""Calculate legacy alignment metrics from a BAM using samtools view.

This is a Python port of the legacy Picard metrics script.  It preserves the
same positional CLI and output labels so existing summary parsing keeps working.
"""

import argparse
import os
import re
import subprocess
import sys


def percent(numerator, denominator):
    if denominator == 0:
        return "0.00%"
    return "{:.2f}%".format(100 * numerator / denominator)


def parse_optional_int(pattern, line):
    match = re.search(pattern, line)
    if match:
        return int(match.group(1))
    return None


def collect_metrics(bam, samtools):
    metrics = {
        "clean_reads": 0,
        "clean_bases": 0,
        "mapped_reads": 0,
        "mapped_bases": 0,
        "dup_reads": 0,
        "mismatch_bases": 0,
        "uniq_reads": 0,
        "uniq_bases": 0,
    }

    command = [samtools, "view", "-X", bam]
    with subprocess.Popen(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True,
    ) as proc:
        assert proc.stdout is not None
        for line in proc.stdout:
            line = line.rstrip("\n")
            if not line or line.startswith("@"):
                continue

            fields = line.split("\t")
            if len(fields) < 11:
                continue

            flag = fields[1]
            mapq = int(fields[4])
            sequence_length = len(fields[9])

            metrics["clean_reads"] += 1
            metrics["clean_bases"] += sequence_length

            if "u" in flag:
                continue

            aligned_bases = parse_optional_int(r"AS:i:(\d+)", line)
            if aligned_bases is None:
                aligned_bases = sequence_length

            if mapq >= 1:
                metrics["uniq_reads"] += 1
                metrics["uniq_bases"] += aligned_bases

            metrics["mapped_reads"] += 1
            metrics["mapped_bases"] += aligned_bases

            if "d" in flag:
                metrics["dup_reads"] += 1

            mismatches = parse_optional_int(r"NM:i:(\d+)", line)
            if mismatches is not None:
                metrics["mismatch_bases"] += mismatches

        stderr = proc.stderr.read() if proc.stderr is not None else ""
        return_code = proc.wait()
        if return_code != 0:
            sys.stderr.write(stderr)
            raise SystemExit(return_code)

    return metrics


def sample_name(bam):
    return os.path.basename(bam).split(".")[0]


def main():
    parser = argparse.ArgumentParser(
        usage="%(prog)s <in.picard.bam> <samtools>",
        description="Calculate legacy alignment metrics from a BAM.",
    )
    parser.add_argument("bam")
    parser.add_argument("samtools")
    args = parser.parse_args()

    metrics = collect_metrics(args.bam, args.samtools)
    clean_reads = metrics["clean_reads"]
    mapped_reads = metrics["mapped_reads"]
    mapped_bases = metrics["mapped_bases"]

    print("Sample\t{}".format(sample_name(args.bam)))
    print("Clean reads\t{}".format(clean_reads))
    print("Clean bases(bp)\t{}".format(metrics["clean_bases"]))
    print("Mapped reads\t{}".format(mapped_reads))
    print("Mapped bases(bp)\t{}".format(mapped_bases))
    print("Mapping rate\t{}".format(percent(mapped_reads, clean_reads)))
    print("Uniq reads\t{}".format(metrics["uniq_reads"]))
    print("Uniq bases(bp)\t{}".format(metrics["uniq_bases"]))
    print("Unique rate\t{}".format(percent(metrics["uniq_reads"], mapped_reads)))
    print("Duplicate reads\t{}".format(metrics["dup_reads"]))
    print("Duplicate rate\t{}".format(percent(metrics["dup_reads"], mapped_reads)))
    print("Mismatch bases(bp)\t{}".format(metrics["mismatch_bases"]))
    print("Mismatch rate\t{}".format(percent(metrics["mismatch_bases"], mapped_bases)))


if __name__ == "__main__":
    main()
