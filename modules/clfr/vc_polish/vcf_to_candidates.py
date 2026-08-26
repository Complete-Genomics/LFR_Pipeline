#!/usr/bin/env python3
"""Production adapter: turn a called VCF's SNP records into the candidates TSV
that 02_extract_features.py expects (chrom/pos/ref/alt/label).

01_make_candidates.py (from cLFR_eval) is a TRAINING-data tool: it requires
--truth-vcf + --confident-bed to LABEL each site true/error against GIAB. At
inference time there is no GIAB truth for the sample being polished -- the
candidate positions instead come from whatever caller already ran
(HaplotypeCaller in this pipeline's make_vcf.smk). This script bridges that
gap: it emits the same 5-column header 01_make_candidates.py produces, with
label filled as a placeholder (-1, never used downstream -- 02 only carries
it through unread, and 04_apply_rescore.py never reads it either).

Only bi-allelic SNVs are emitted (ref/alt both single-base); indels and
multi-allelic sites are skipped -- the shipped model was trained on SNVs only.

Usage
-----
  python vcf_to_candidates.py --vcf sample.snp.vcf.gz --out candidates.tsv
"""
import argparse
import sys

import pysam

BASES = {"A", "C", "G", "T"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vcf", required=True, help="called VCF (bgzip+tabix or plain)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--regions", nargs="*", default=None,
                     help="restrict to these chroms; default = whole VCF")
    args = ap.parse_args()

    vf = pysam.VariantFile(args.vcf)
    n_in = n_out = 0
    with open(args.out, "w") as out:
        out.write("chrom\tpos\tref\talt\tlabel\n")
        chroms = args.regions if args.regions else vf.header.contigs.keys()
        for chrom in chroms:
            try:
                records = vf.fetch(chrom)
            except ValueError:
                continue  # contig not present in this VCF
            for rec in records:
                n_in += 1
                ref = rec.ref
                if ref is None or len(ref) != 1 or ref.upper() not in BASES:
                    continue
                for alt in (rec.alts or ()):
                    if len(alt) != 1 or alt.upper() not in BASES:
                        continue
                    out.write(f"{chrom}\t{rec.pos}\t{ref.upper()}\t{alt.upper()}\t-1\n")
                    n_out += 1
    sys.stderr.write(f"[vcf_to_candidates] {n_in} VCF records -> {n_out} SNV candidates -> {args.out}\n")


if __name__ == "__main__":
    main()
