#!/usr/bin/env bash
# qc_adapter_scan.sh -- standalone adapter QC, no pipeline wiring, no adapter
# reference file required. Run this on a batch's raw fastq BEFORE deciding
# whether to enable trim_reads' adapter_ref_5p switch (see
# LFR_Pipeline/memory/project_lfr_chr18_troubleshoot.md) -- confirm with
# wet-lab first, then flip the config switch; this script does not touch
# the workflow.
#
# Checks two independent, unrelated failure modes:
#   1) Classic 3' read-through adapter (insert shorter than read length) --
#      auto-detected via bbmerge's overlap-based adapter finder. No need to
#      know or guess the adapter sequence in advance.
#   2) 5'-anchored fixed-position contamination (like "adapter183" found in
#      this project: GAGACGTTCTCGACTCAGCAGAGGG at R1 pos1-25, present
#      regardless of insert length) -- detected by scanning per-position
#      base composition on a true random sample of reads and flagging
#      positions where one base dominates far above the ~25% random
#      background. bbmerge/tbo/tpe CANNOT catch this class: it has nothing
#      to do with R1/R2 overlap, it's baked into the read from position 1.
#
# Usage:
#   ./qc_adapter_scan.sh R1.fq.gz R2.fq.gz [output_dir] [sample_size] [scan_len] [known_adapter.fa ...]
#
# Example:
#   ./qc_adapter_scan.sh data/split_read.1.fq.gz data/split_read.2.fq.gz \
#       qc_out 5000 40 \
#       /path/to/LFR_Pipeline/config/adapters/mgi_dnbseq_adapters.fa \
#       /path/to/LFR_Pipeline/config/adapters/mgi_dnbseq_5p_contam.fa

set -euo pipefail

if [[ $# -lt 2 ]]; then
    echo "Usage: $0 R1.fq.gz R2.fq.gz [output_dir] [sample_size] [scan_len] [known_adapter.fa ...]" >&2
    exit 1
fi

R1="$1"
R2="$2"
OUTDIR="${3:-qc_adapter_out}"
SAMPLE_N="${4:-5000}"
SCAN_LEN="${5:-40}"
shift $(( $# < 5 ? $# : 5 ))
KNOWN_ADAPTERS=("$@")

mkdir -p "$OUTDIR"

echo "=================================================================="
echo "[1/3] 3' read-through adapter auto-detection (bbmerge overlap)"
echo "=================================================================="
if command -v bbmerge.sh >/dev/null 2>&1; then
    bbmerge.sh in1="$R1" in2="$R2" outa="$OUTDIR/detected_3p_adapters.fa" \
        > "$OUTDIR/bbmerge.log" 2>&1 || true
    if [[ -s "$OUTDIR/detected_3p_adapters.fa" ]]; then
        echo "Detected 3' adapter(s) -- see $OUTDIR/detected_3p_adapters.fa:"
        cat "$OUTDIR/detected_3p_adapters.fa"
    else
        echo "No adapter sequence recovered by overlap detection."
        echo "(Either the data is clean, or inserts are consistently longer"
        echo " than the read length so R1/R2 never overlap into adapter --"
        echo " that does NOT rule out a 5'-anchored contaminant; see step 3.)"
    fi
else
    echo "SKIPPED: bbmerge.sh not found in PATH."
fi
echo ""

echo "=================================================================="
echo "[2/3] True random sampling ($SAMPLE_N reads each from R1/R2)"
echo "=================================================================="
# Deliberately NOT using `head` -- a coordinate-sorted-bam-style bias trap
# already documented in this project's memory (head only ever samples
# whichever records happen to sit first in the file, which for raw fastq
# is fine, but for consistency with the bam-based version of this check
# elsewhere in the debug history, and to avoid any accidental lane/tile
# ordering bias, sample uniformly at random across the whole file).
sample_fastq() {
    local infile="$1" outfile="$2" n="$3"
    if command -v seqtk >/dev/null 2>&1; then
        seqtk sample -s100 "$infile" "$n" > "$outfile"
    else
        # `sort` emits every record before `head -n` stops reading, so once
        # head has enough lines it closes the pipe and SIGPIPEs sort/paste/
        # zcat -- that's expected here, not a real failure, so pipefail must
        # be off for just this one pipeline (same gotcha as
        # merge_14bp_umi.sh's head -1 probe elsewhere in this project).
        set +o pipefail
        zcat -f "$infile" \
            | paste - - - - \
            | awk 'BEGIN{srand()}{print rand()"\t"$0}' \
            | sort -n -k1,1 \
            | cut -f2- \
            | head -n "$n" \
            | tr '\t' '\n' \
            > "$outfile"
        set -o pipefail
    fi
}
sample_fastq "$R1" "$OUTDIR/r1_sample.fq" "$SAMPLE_N"
sample_fastq "$R2" "$OUTDIR/r2_sample.fq" "$SAMPLE_N"
echo "R1 sampled: $(( $(wc -l < "$OUTDIR/r1_sample.fq") / 4 )) reads"
echo "R2 sampled: $(( $(wc -l < "$OUTDIR/r2_sample.fq") / 4 )) reads"
echo ""

echo "=================================================================="
echo "[3/3] Per-position base composition scan (pos 1-$SCAN_LEN)"
echo "=================================================================="
python3 - "$OUTDIR/r1_sample.fq" "$OUTDIR/r2_sample.fq" "$SCAN_LEN" "${KNOWN_ADAPTERS[@]}" <<'PYEOF'
import sys
from collections import Counter

r1_path, r2_path, scan_len = sys.argv[1], sys.argv[2], int(sys.argv[3])
known_adapter_paths = sys.argv[4:]

# Flag threshold: random background for one-of-four bases is 25%. A real
# fixed-position contaminant/adapter sits far above that (this project's
# adapter183 measured 70-86% at its anchored positions). 40% is a
# conservative trip-wire that won't false-fire on ordinary base-composition
# skew (GC bias etc rarely pushes a single base past ~35% at one position).
FLAG_THRESHOLD = 0.40
MIN_RUN_LEN = 5  # a single flagged position is likely noise; a contiguous
                 # run of >=5 is the signature of a real fixed sequence.

def read_seqs(path):
    seqs = []
    with open(path) as f:
        for i, line in enumerate(f):
            if i % 4 == 1:
                seqs.append(line.strip())
    return seqs

def scan(seqs, label):
    n = len(seqs)
    if n == 0:
        print(f"  {label}: no reads sampled, skipping")
        return []
    dominant = []
    for pos in range(scan_len):
        c = Counter(s[pos] for s in seqs if len(s) > pos)
        total = sum(c.values())
        if total == 0:
            dominant.append((None, 0.0))
            continue
        base, cnt = c.most_common(1)[0]
        dominant.append((base, cnt / total))
    print(f"  {label} (n={n}):")
    line = "    pos:  " + " ".join(f"{p+1:>4}" for p in range(scan_len))
    print(line)
    line = "    base: " + " ".join(f"{(b or '-'):>4}" for b, _ in dominant)
    print(line)
    line = "    pct:  " + " ".join(f"{f*100:>3.0f}%" for _, f in dominant)
    print(line)
    # find contiguous flagged runs
    runs = []
    run_start = None
    for pos, (base, frac) in enumerate(dominant):
        if frac >= FLAG_THRESHOLD:
            if run_start is None:
                run_start = pos
        else:
            if run_start is not None and pos - run_start >= MIN_RUN_LEN:
                runs.append((run_start, pos))
            run_start = None
    if run_start is not None and scan_len - run_start >= MIN_RUN_LEN:
        runs.append((run_start, scan_len))
    return runs

def consensus_seq(seqs, start, end):
    out = []
    for pos in range(start, end):
        c = Counter(s[pos] for s in seqs if len(s) > pos)
        if not c:
            out.append("N")
            continue
        out.append(c.most_common(1)[0][0])
    return "".join(out)

def load_known_adapters(paths):
    adapters = {}
    for p in paths:
        try:
            with open(p) as f:
                name, seq = None, []
                for line in f:
                    line = line.strip()
                    if not line:
                        continue
                    if line.startswith(">"):
                        if name:
                            adapters[name] = "".join(seq)
                        name, seq = line[1:], []
                    else:
                        seq.append(line)
                if name:
                    adapters[name] = "".join(seq)
        except OSError as e:
            print(f"  WARNING: could not read {p}: {e}")
    return adapters

def best_match(query, adapters):
    best_name, best_len = None, 0
    rc = query.translate(str.maketrans("ACGT", "TGCA"))[::-1]
    for name, seq in adapters.items():
        for q in (query, rc):
            for a, b in ((q, seq), (seq, q)):
                for i in range(len(a)):
                    for j in range(len(b)):
                        k = 0
                        while i + k < len(a) and j + k < len(b) and a[i + k] == b[j + k]:
                            k += 1
                        if k > best_len:
                            best_len, best_name = k, name
    return best_name, best_len

r1_seqs = read_seqs(r1_path)
r2_seqs = read_seqs(r2_path)

print("R1:")
r1_runs = scan(r1_seqs, "R1")
print("R2:")
r2_runs = scan(r2_seqs, "R2")

print("")
print("Flagged contiguous fixed-position runs (>= {}bp at >= {:.0f}% dominance):".format(
    MIN_RUN_LEN, FLAG_THRESHOLD * 100))

known = load_known_adapters(known_adapter_paths) if known_adapter_paths else {}

any_flagged = False
for label, seqs, runs in (("R1", r1_seqs, r1_runs), ("R2", r2_seqs, r2_runs)):
    for start, end in runs:
        any_flagged = True
        cons = consensus_seq(seqs, start, end)
        print(f"  {label} pos{start+1}-{end}: {cons}")
        if known:
            name, length = best_match(cons, known)
            if name and length >= max(10, (end - start) // 2):
                print(f"    -> matches known adapter '{name}' (longest common substring {length}bp) -- already cataloged, check trim direction (ktrim=l for 5'-anchored, ktrim=r for 3' read-through)")
            else:
                print(f"    -> UNKNOWN (best match against provided reference(s): {name or 'none'}, {length}bp overlap) -- confirm with wet-lab before wiring into the pipeline")
        else:
            print(f"    -> no reference adapter file(s) provided to compare against")

if not any_flagged:
    print("  none -- no fixed-position contamination signature detected at this threshold")

print("")
print("Note: this scan only catches FIXED-POSITION signals (same base(s) at the")
print("same coordinate across reads). Ordinary 3' read-through adapters at variable")
print("insert-length-dependent positions won't show up here -- that's what step 1")
print("(bbmerge) is for. The two checks are complementary, not redundant.")
PYEOF

echo ""
echo "Done. Outputs in $OUTDIR/"
