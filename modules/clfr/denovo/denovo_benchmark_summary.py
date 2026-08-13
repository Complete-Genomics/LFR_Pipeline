#!/usr/bin/env python3
"""Collapse this run's Benchmarks/denovo.*.txt into one table.

Snakemake's `benchmark:` directive already records wall-clock, max RSS, and
mean CPU load per rule (one TSV per rule, via psutil) -- this just prints them
side by side, the same at-a-glance view Nextflow's `-with-report` gives across
all processes in one file. No new tracking mechanism: every number here comes
from a `benchmark:` file some rule in denovo_preprocess.smk/denovo_olc.smk
already wrote.

Usage:
    python3 denovo_benchmark_summary.py [Benchmarks/]
"""
import argparse
import csv
import glob
import os


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("benchmarks_dir", nargs="?", default="Benchmarks")
    ap.add_argument("--out", help="also write the table here (TSV)")
    args = ap.parse_args()

    rows = []
    for path in sorted(glob.glob(os.path.join(args.benchmarks_dir, "denovo.*.txt"))):
        rule = os.path.basename(path)[len("denovo."):-len(".txt")]
        with open(path) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for line in reader:
                rows.append((rule, line))

    if not rows:
        print(f"no denovo.*.txt benchmark files found under {args.benchmarks_dir}/ "
              f"-- did the run include rules with a benchmark: directive?")
        return

    fields = ["rule", "s", "h:m:s", "max_rss", "mean_load", "cpu_time", "io_in", "io_out"]
    widths = {f: max(len(f), max(len(rule if f == "rule" else line.get(f, "")) for rule, line in rows))
              for f in fields}

    def fmt_row(values):
        return "  ".join(str(v).ljust(widths[f]) for f, v in zip(fields, values))

    print(fmt_row(fields))
    out_fh = open(args.out, "w", newline="") if args.out else None
    writer = csv.writer(out_fh, delimiter="\t") if out_fh else None
    if writer:
        writer.writerow(fields)
    for rule, line in rows:
        values = [rule] + [line.get(f, "") for f in fields[1:]]
        print(fmt_row(values))
        if writer:
            writer.writerow(values)
    if out_fh:
        out_fh.close()
        print(f"\nwritten to {args.out}")


if __name__ == "__main__":
    main()
