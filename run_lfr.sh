#!/usr/bin/env bash
# run_lfr.sh — execute the LFR pipeline via Snakemake
#
# Usage:
#   cd <analysis_dir>          # must contain config.yaml
#   bash /path/to/run_lfr.sh              # stLFR, 20 cores
#   bash /path/to/run_lfr.sh clfr 40      # cLFR, 40 cores
#
# Prerequisites:
#   - Snakemake installed and on PATH (conda activate snakemake_env)
#   - config.yaml copied from config/stlfr.yaml or config/clfr.yaml and filled in

set -euo pipefail

PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MODE="stlfr"
THREADS=20

if [[ "${1:-}" == "stlfr" || "${1:-}" == "clfr" ]]; then
    MODE="$1"
    shift
fi

if [[ "${1:-}" =~ ^[0-9]+$ ]]; then
    THREADS="$1"
    shift
fi

if [[ "$MODE" == "clfr" ]]; then
    SNAKEFILE="$PIPELINE_DIR/workflows/clfr.smk"
else
    SNAKEFILE="$PIPELINE_DIR/workflows/stlfr.smk"
fi

echo "Running $MODE pipeline"
echo "Snakefile : $SNAKEFILE"
echo "Threads   : $THREADS"
echo "Config    : $(pwd)/config.yaml"
echo ""

snakemake \
    --snakefile "$SNAKEFILE" \
    --cores "$THREADS" \
    --rerun-incomplete \
    --latency-wait 60 \
    --printshellcmds \
    "$@"
