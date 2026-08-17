#!/usr/bin/env python3
"""Apply the hs1 salvage filter to an existing barcode-grouped TSV.

This is a small explicit-use wrapper around the production salvage primitive in
``denovo_noisy_preprocess_qc.py``.  It does not run the sample anomaly gate and
it never overwrites the input.  The input must already be grouped by barcode,
as ``denovo/data_R2_sorted.tsv`` is.

Default behavior, matching the hs1 salvage decision:

* discard reads shorter than 300 bp;
* retain at most the first 300 eligible reads per UMI in input order.

Example::

    python tmp_filter.py \
        --input hs1/denovo/data_R2_sorted.tsv \
        --output hs1/denovo/data_R2_tmp_filtered.tsv
"""

import argparse
import os

import denovo_noisy_preprocess_qc as qc


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True,
                        help="barcode-grouped data_R2_sorted.tsv")
    parser.add_argument("--output", required=True,
                        help="filtered TSV; must differ from --input")
    parser.add_argument("--min-read-length", type=int, default=300,
                        help="minimum retained read length [300]")
    parser.add_argument("--max-reads-per-umi", type=int, default=300,
                        help="maximum retained eligible reads per UMI [300]")
    args = parser.parse_args(argv)

    source = os.path.abspath(args.input)
    destination = os.path.abspath(args.output)
    if source == destination:
        parser.error("--output must differ from --input; original reads are preserved")
    if args.min_read_length < 1:
        parser.error("--min-read-length must be positive")
    if args.max_reads_per_umi < 1:
        parser.error("--max-reads-per-umi must be positive")
    if not os.path.isfile(source):
        parser.error("input does not exist: {}".format(source))

    parent = os.path.dirname(destination)
    if parent:
        os.makedirs(parent, exist_ok=True)
    kept = qc.filter_tsv(source, destination, args.min_read_length,
                         args.max_reads_per_umi)
    print("kept_reads\t{}".format(kept))
    print("output\t{}".format(destination))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
