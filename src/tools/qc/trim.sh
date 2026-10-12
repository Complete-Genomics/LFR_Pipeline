#!/usr/bin/env bash
# Per-batch trim parameters, from the RAW reads, in one command. Prints the trim_* / adapter_ref_* values to put in
# the run's config.yaml (trim_reads then applies them); replaces the hand-written trim.sh of each batch.
#
#   cd <the run dir>; bash trim.sh          # uses ./config.yaml and ./data (data/read_2.fq.gz; PE also data/read_1.fq.gz)
# Optional overrides, all positional and in this order:
#   bash trim.sh [config.yaml] [data_dir_or_raw.fq.gz] [outdir=trim_scan] [raw_r2.fq.gz]
#   (the old order  raw.fq.gz config.yaml  also works: the argument ending in .yaml is the config)
#
# Needs: conda activate lfr_pipeline (snakemake, bbduk.sh, numpy, matplotlib, pyyaml). Nothing is written to config.yaml.
#
# What it does (the first N_READS raw reads, default 100000; the pipeline's own rules do the split and the trim):
#   A  split the sample                                  -> outdir/test_split_read.{1,2}.fq.gz   (QC copy of the split reads;
#                                                           test_ prefix = QC output, not the full-data split_read.*)
#   B  adapter evidence on the split reads (qc_adapter_scan.sh report + qc_adapter_decide.py rules)
#        adapter_ref_5p / trim_5p_restrict_left (adapter183 anchored at the read start), adapter_ref_3p (linker)
#   C  trim with those keys, then look at what is left   -> outdir/test_split_read_2_trimmed.fastq.gz
#        trim_mgi_adapter (+ adapter_ref_k 18) if MGI fragments remain, trim_polya if polyA tails remain
#   D  quality scan with the decided adapter keys (qc_quality_scan.py): trim_q / trim_maq (the 3' low-quality tail is
#      largely adapter cassette, so the quality keys are decided last)
# Prints ALL trim keys: the ones decided from data, and trim_qtrim / trim_minlen / trim_maxns (not data-decidable) with
# the value from your config.yaml (or the pipeline default); each line says which it is.
# Thresholds are rules of thumb (see qc_adapter_decide.py); the evidence is in outdir/trim_report.txt, the plots are
# outdir/quality_scan.png and outdir/low_quality_distribution.png. Quality is one
# side only: check a changed trim on the aligned data (low_cov fraction, soft-clip length).

set -euo pipefail

USAGE="usage: cd <run dir>; bash $0     (needs ./config.yaml and ./data/;  optional: $0 [config.yaml] [data_dir_or_raw.fq.gz] [outdir] [raw_r2.fq.gz])"
A="${1:-config.yaml}"; B="${2:-data}"
if [[ "$A" == *.yaml || "$A" == *.yml ]]; then CONFIG="$A"; SRC="$B"; else CONFIG="$B"; SRC="$A"; fi
OUT="${3:-trim_scan}"
RAW_R2="${4:-}"
N_READS="${N_READS:-100000}"      # raw reads taken for the scans
ADAPTER_SAMPLE="${ADAPTER_SAMPLE:-5000}"

TOOLS="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ADAPTERS="${ADAPTERS:-$TOOLS/../../../config/adapters}"
[[ -f "$CONFIG" ]] || { echo "config not found (or not a file): $CONFIG (run from the run dir that holds config.yaml)" >&2; echo "$USAGE" >&2; exit 1; }
if [[ -d "$SRC" ]]; then          # the run's data/ dir: SE reads live in read_2, PE has read_1 + read_2
    SEQ="$(python -c "import sys,yaml; print(str(yaml.safe_load(open(sys.argv[1]))['params'].get('sequence_type','se')).lower())" "$CONFIG")"
    if [[ "$SEQ" == "se" ]]; then RAW="$SRC/read_2.fq.gz"; RAW_R2=""; else RAW="$SRC/read_1.fq.gz"; RAW_R2="$SRC/read_2.fq.gz"; fi
else
    RAW="$SRC"
fi
[[ -f "$RAW" ]] || { echo "raw reads not found: $RAW" >&2; echo "$USAGE" >&2; exit 1; }
[[ -z "$RAW_R2" || -f "$RAW_R2" ]] || { echo "raw R2 not found: $RAW_R2" >&2; exit 1; }
R2_ARG=(); [[ -n "$RAW_R2" ]] && R2_ARG=(--raw_r2 "$RAW_R2")
# preflight: the split needs no bbduk, the trim does -- fail now, not in step C
BBDUK_CMD="$(python -c "import sys,yaml; print(str(yaml.safe_load(open(sys.argv[1]))['params'].get('bbduk','bbduk.sh') or 'bbduk.sh').split()[0])" "$CONFIG")"
for tool in snakemake; do
    command -v "$tool" >/dev/null 2>&1 || { echo "$tool not found in PATH: conda activate lfr_pipeline" >&2; exit 1; }
done
if ! command -v "$BBDUK_CMD" >/dev/null 2>&1 && ! command -v bbduk.sh >/dev/null 2>&1; then
    echo "bbduk not found (config params.bbduk = $BBDUK_CMD, and bbduk.sh is not in PATH): conda activate lfr_pipeline" >&2; exit 1
fi
mkdir -p "$OUT"
REPORT="$OUT/trim_report.txt"; : > "$REPORT"

scan() {   # scan <config> <outdir> [extra qc_quality_scan.py args]: log in <outdir>.log
    local cfg="$1" o="$2"; shift 2
    python "$TOOLS/qc_quality_scan.py" "$RAW" --raw --config "$cfg" -o "$o" -n "$N_READS" \
        ${R2_ARG[@]+"${R2_ARG[@]}"} "$@" > "$o.log" 2>&1 || { echo "FAILED: see $o.log (last lines:)" >&2; tail -15 "$o.log" >&2; exit 1; }
}
decide() { # decide <subcommand> <reads> <kv-out>: evidence -> report, decided keys -> <kv-out>
    python "$TOOLS/qc_adapter_decide.py" "$1" --reads "$2" --yaml_out "$3" | tee -a "$REPORT" | sed 's/^/   /'
}
kv() { sed 's/=/: /' "$1"; }

echo "config $CONFIG | raw $RAW ${RAW_R2:+| raw R2 $RAW_R2} | outdir $OUT | $N_READS reads"

echo "A/ split the sample"
scan "$CONFIG" "$OUT/qs_split" --split_only
cp "$OUT/qs_split/test_split_read.2.fq.gz" "$OUT/test_split_read.2.fq.gz"
[[ -n "$RAW_R2" ]] && cp "$OUT/qs_split/test_split_read.1.fq.gz" "$OUT/test_split_read.1.fq.gz"
SPLIT2="$OUT/test_split_read.2.fq.gz"; SPLIT1="${RAW_R2:+$OUT/test_split_read.1.fq.gz}"; SPLIT1="${SPLIT1:-$SPLIT2}"   # SE: same file twice

echo "B/ adapter evidence on the split reads"
bash "$TOOLS/qc_adapter_scan.sh" "$SPLIT1" "$SPLIT2" "$OUT/qc_adapter" "$ADAPTER_SAMPLE" 40 20 \
    "$ADAPTERS/mgi_dnbseq_adapters.fa" "$ADAPTERS/mgi_dnbseq_5p_contam.fa" "$ADAPTERS/clfr_se600_3p_linker.fa" \
    > "$OUT/qc_adapter.log" 2>&1 || true            # full scan (also finds adapters not in the catalogue): qc_adapter.log
decide split "$SPLIT2" "$OUT/d_split.kv"
SET1=($(cat "$OUT/d_split.kv"))

echo "C/ trim with those keys and look at what is left"
python "$TOOLS/qc_adapter_decide.py" setcfg --config "$CONFIG" --out "$OUT/config.c1.yaml" "${SET1[@]}" trim_mgi_adapter=false
scan "$OUT/config.c1.yaml" "$OUT/c1" --candidates 6:0 --no_plot
cp "$OUT/c1/test_split_read_2_trimmed.fastq.gz" "$OUT/test_split_read_2_trimmed.fastq.gz"
decide mgi "$OUT/test_split_read_2_trimmed.fastq.gz" "$OUT/d_mgi.kv"
SET2=($(cat "$OUT/d_mgi.kv"))
python "$TOOLS/qc_adapter_decide.py" setcfg --config "$OUT/config.c1.yaml" --out "$OUT/config.c2.yaml" "${SET2[@]}"
if grep -q "trim_mgi_adapter=true" "$OUT/d_mgi.kv"; then      # polyA only shows once the MGI cassette behind it is gone
    scan "$OUT/config.c2.yaml" "$OUT/c2" --candidates 6:0 --no_plot
    cp "$OUT/c2/test_split_read_2_trimmed.fastq.gz" "$OUT/test_split_read_2_trimmed.fastq.gz"
fi
decide polya "$OUT/test_split_read_2_trimmed.fastq.gz" "$OUT/d_polya.kv"
SET3=($(cat "$OUT/d_polya.kv"))
python "$TOOLS/qc_adapter_decide.py" setcfg --config "$OUT/config.c2.yaml" --out "$OUT/config.final.yaml" "${SET3[@]}"

echo "D/ quality scan with the decided adapter keys"
scan "$OUT/config.final.yaml" "$OUT/qs_final"
cat "$OUT/qs_final/report.txt" >> "$REPORT"
for png in quality_scan low_quality_distribution; do   # the plots, next to trim_params.yaml
    [[ -f "$OUT/qs_final/$png.png" ]] && cp "$OUT/qs_final/$png.png" "$OUT/$png.png"
done
TQ="$(awk '$1=="trim_q:"{print $2}' "$OUT/qs_final/suggested_config.yaml")"
TM="$(awk '$1=="trim_maq:"{print $2}' "$OUT/qs_final/suggested_config.yaml")"
REASON="$(grep '^recommendation' "$OUT/qs_final/report.txt" | sed 's/^recommendation: [^(]*//')"

python "$TOOLS/qc_adapter_decide.py" params --config "$CONFIG" --kv "$OUT/d_split.kv" "$OUT/d_mgi.kv" "$OUT/d_polya.kv" \
    --trim_q "$TQ" --trim_maq "$TM" --reason "$REASON" --label "from $N_READS raw reads of $(basename "$RAW"); evidence: $REPORT" \
    --out "$OUT/trim_params.yaml"

echo
echo "================ trim parameters for config.yaml (also in $OUT/trim_params.yaml) ================"
cat "$OUT/trim_params.yaml"
echo
echo "plots: $OUT/quality_scan.png (quality by position, per-read mean Q, trim trade-off)  $OUT/low_quality_distribution.png"
echo "report: $REPORT"
