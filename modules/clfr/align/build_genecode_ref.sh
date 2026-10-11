#!/usr/bin/env bash
# Build the combined GENCODE reference for genecode.smk (run ONCE per GENCODE version, not per sample):
#   transcripts + rRNA + genome decoy -> one gz fasta, streamed (no decompressed intermediate on disk).
# Headers are rewritten to the GENCODE style  id|id|-|-|id|id|0|<biotype>|  with biotype rRNA / genome_decoy
# for the custom records; genecode_classify.py reads only the id and the biotype field.
#
# Usage:
#   build_genecode_ref.sh -t gencode.v50.transcripts.fa.gz -r human_rRNA.fa -g GRCh38.primary_assembly.genome.fa.gz [-o out.fa.gz]
# Default output (what genecode.smk expects): the -t path with .transcripts.fa.gz -> .transcripts_rRNA_decoy.fa.gz
# -r / -g may be omitted (rRNA reads then cannot be counted / genomic reads are forced onto transcripts).
# Env: GZIP (default pigz if installed, else gzip) for the output compression.

set -euo pipefail

TX=""; RRNA=""; GENOME=""; OUT=""
while getopts "t:r:g:o:h" opt; do
    case "$opt" in
        t) TX="$OPTARG" ;;
        r) RRNA="$OPTARG" ;;
        g) GENOME="$OPTARG" ;;
        o) OUT="$OPTARG" ;;
        *) sed -n '2,12p' "$0"; exit 1 ;;
    esac
done
[[ -n "$TX" && -s "$TX" ]] || { echo "missing or empty transcripts fasta (-t)" >&2; exit 1; }
[[ -z "$RRNA" || -s "$RRNA" ]] || { echo "rRNA fasta not found: $RRNA" >&2; exit 1; }
[[ -z "$GENOME" || -s "$GENOME" ]] || { echo "genome fasta not found: $GENOME" >&2; exit 1; }
if [[ -z "$OUT" ]]; then
    case "$TX" in
        *.transcripts.fa.gz) OUT="${TX/.transcripts.fa.gz/.transcripts_rRNA_decoy.fa.gz}" ;;
        *) OUT="$TX.rRNA_decoy.fa.gz" ;;
    esac
fi
GZ="${GZIP:-$(command -v pigz || command -v gzip)}"

# cat a plain or gzipped fasta
fa() { case "$1" in *.gz) gzip -dc "$1" ;; *) cat "$1" ;; esac; }
# keep only the first word of each header and rewrite it; sequence lines are untouched
hdr() { sed "s/^>\([^ ]*\).*/>\1|\1|-|-|\1|\1|0|$1|/"; }

# write to .tmp then mv: a half-written reference is never visible under the final name
TMP="$OUT.tmp.$$"
trap 'rm -f "$TMP"' EXIT
{
    fa "$TX"
    if [[ -n "$RRNA" ]]; then fa "$RRNA" | hdr rRNA; else echo "WARNING: no rRNA fasta (-r); rRNA reads will not be countable" >&2; fi
    if [[ -n "$GENOME" ]]; then fa "$GENOME" | hdr genome_decoy; else echo "WARNING: no genome fasta (-g); no genome decoy" >&2; fi
} | "$GZ" > "$TMP"
gzip -t "$TMP"
mv "$TMP" "$OUT"
trap - EXIT
echo "built $OUT ($(du -h "$OUT" | cut -f1))"
