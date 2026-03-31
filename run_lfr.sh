#!/usr/bin/env bash
# run_lfr.sh — execute the LFR pipeline via Snakemake
#
# Usage:
#   cd <analysis_dir>          # must contain config.yaml
#   bash /path/to/run_lfr.sh   # stLFR
#   bash /path/to/run_lfr.sh clfr   # cLFR
#
# Prerequisites:
#   - Snakemake installed and on PATH (conda activate snakemake_env)
#   - config.yaml copied from config/stlfr.yaml or config/clfr.yaml and filled in
#   - src_dir in config.yaml points to this repository root

set -euo pipefail

MODE="${1:-stlfr}"
THREADS="${2:-20}"
PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

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
