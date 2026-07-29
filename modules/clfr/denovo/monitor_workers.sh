#!/usr/bin/env bash
# Samples denovo_seed_olc.py worker CPU utilization every INTERVAL seconds
# and logs how many workers are actually busy (>50% CPU) vs idle at each
# sample. Run this alongside a real production job to see whether the
# active-worker count stays near num_processes throughout, or dips/
# fluctuates -- without needing to guess from occasional manual `ps` checks.
#
# Usage:
#   ./monitor_workers.sh [interval_seconds] > worker_monitor.log 2>&1 &
#   # ... let the real denovo_seed_olc.py job run ...
#   # afterwards, summarize:
#   grep busy_workers worker_monitor.log | awk -F= '{print $2}' | sort -n | uniq -c

INTERVAL="${1:-10}"

while true; do
    ts=$(date '+%Y-%m-%d %H:%M:%S')
    snapshot=$(ps aux | grep '[d]enovo_seed_olc.py')
    n_total=$(echo "$snapshot" | grep -c .)
    n_busy=$(echo "$snapshot" | awk '{ if ($3+0 > 50) c++ } END { print c+0 }')
    echo "ts=${ts} total_workers=${n_total} busy_workers=${n_busy}"
    sleep "$INTERVAL"
done
