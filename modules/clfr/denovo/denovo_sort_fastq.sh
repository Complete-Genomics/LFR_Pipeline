#!/usr/bin/env bash
# Barcode-sort a split/trimmed FASTQ, preserving quality.
#
# denovo_preprocess.smk's reformat_fasta2 already sorts by barcode, but it
# drops the quality string on the way (it emits a 2-column TSV of header+seq),
# so it cannot feed anything that needs quality -- specifically
# denovo_read_features.py, whose model uses five quality-derived features.
# This produces the same barcode grouping while keeping all four FASTQ lines.
#
# Implemented with samtools rather than sort(1), because samtools already
# sorts by an arbitrary aux tag and the split-read header carries the barcode
# as exactly that (BX:Z:). Measured on hs6 (132k reads): 0.8-1.1s versus 7.0s
# for an equivalent awk+sort(1) pipeline -- ~7x, plus it is multi-threaded and
# spills to disk, so it does not degrade at 3M-UMI scale. samtools is already
# a pipeline dependency (config params.samtools).
#
# NOT usable here: `seqkit sort`. It sorts on the whole sequence ID, and the
# barcode sits in the MIDDLE of these ids ("@<runid>#<barcode>/2"), so sorting
# by id groups by run/tile, not by barcode. It also loads the whole FASTQ into
# memory (its low-memory two-pass mode is FASTA-only).
#
# Round-trip is byte-exact on real data (132,323 records, 0 differences in
# header/sequence/quality) and deterministic across runs. Two details worth
# knowing:
#   - `samtools import` drops the "/2" suffix into the pairing flag, and
#     `samtools fastq` restores it; the intermediate SAM read name is bare.
#   - Read order WITHIN a barcode is not preserved from the input (it is
#     deterministic, just different). That matters because
#     denovo_read_filter.py examines only the first --max-reads of each pool,
#     so which reads it inspects can shift versus an unsorted-input run.
#
# Usage: denovo_sort_fastq.sh <in.fastq[.gz]> <out.fastq[.gz]> [threads] [samtools]
set -euo pipefail

in=${1:?usage: denovo_sort_fastq.sh <in.fastq[.gz]> <out.fastq[.gz]> [threads] [samtools]}
out=${2:?usage: denovo_sort_fastq.sh <in.fastq[.gz]> <out.fastq[.gz]> [threads] [samtools]}
threads=${3:-4}
samtools=${4:-samtools}

tmp_prefix="$(dirname "$out")/.denovo_sort_fastq.$$"

# -T '*' carries every header comment field through as an aux tag, which is
# what puts BX:Z: (and anything else upstream added) where `sort -t BX` can
# see it; `samtools fastq -T BX` writes it back out as a header comment.
"$samtools" import -T '*' -@ "$threads" "$in" \
  | "$samtools" sort -t BX -@ "$threads" -O sam -T "$tmp_prefix" \
  | "$samtools" fastq -T BX -@ "$threads" - \
  | { case "$out" in
          *.gz) gzip -c > "$out" ;;
          *)    cat > "$out" ;;
      esac; }
