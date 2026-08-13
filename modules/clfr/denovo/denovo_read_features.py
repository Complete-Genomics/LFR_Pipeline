#!/usr/bin/env python3
"""Per-read feature extraction for the learned read-quality score.

Feeds the model whose prediction denovo_read_filter.py's conflict-graph uses
as a tie-break (denovo.md sec 51/59/62: four samples, +0.03..+0.07 identity
against a greengenes local truth, all significant). Feature set is fixed by
the trained model -- changing or reordering these columns silently invalidates
it, so treat the column list as part of the model artifact, not as tunable.

WHY THIS IS ITS OWN PASS: the features need BOTH quality strings (present only
in the FASTQ, dropped by denovo_preprocess.smk's reformat_fasta2) and the
whole barcode's read pool (pool_kmer_popular_frac). Nothing else in the
pipeline holds both at once.

INPUT MUST BE BARCODE-SORTED. This streams one barcode group at a time so
memory stays bounded by --batch-barcodes rather than by the run: at 3M UMI a
load-everything-then-group approach (which is what the original prototype did)
would need tens of GB of sequence resident. Same requirement, and same reason,
as denovo_read_filter.py's iter_barcode_groups. Note that
data/split_read_2_trimmed.fastq.gz is NOT barcode-sorted -- it comes off the
pre-sort stage -- so it must be sorted before being fed here.

Runtime, measured on hs6 (3000 barcodes / 132k reads, denovo.md sec 63): the
original single-process prototype took 75.5s, which extrapolates to ~23h at
3M UMI and was the single largest cost of adopting the learned score at all.
Two fixes, both here: an O(n) sliding-window minimizer (was O(n*w) -- it alone
was 69% of profile time) and the same process-pool scaffolding the rest of the
denovo modules use.
"""
import argparse
import gzip
import os
import sys
from collections import Counter, defaultdict, deque

K_MINI = 15
W_MINI = 10
K_POOL = 21

# Fixed by the trained model -- see module docstring.
FIELDS = ["read_id", "barcode", "length", "qual_mean", "qual_min",
          "qual_head", "qual_tail", "qual_trend", "hp_frac", "hp_max_run",
          "mini_gap_mean", "mini_gap_max", "mini_gap_var",
          "pool_size", "pool_kmer_popular_frac"]


def minimizer_gaps(seq, k=K_MINI, w=W_MINI):
    """Gaps between successive distinct minimizer positions.

    Monotonic-deque sliding-window minimum: each k-mer is pushed and popped at
    most once, so this is O(n) instead of the O(n*w) that re-scanning every
    window with min() costs. Verified to produce byte-identical output to that
    version on random sequence; 2.4x faster in isolation, and it was 69% of the
    profile, so it dominates this module's total.
    """
    n = len(seq)
    if n < k:
        return []
    m = n - k + 1
    if m < w:
        return []
    kmers = [seq[i:i + k] for i in range(m)]
    dq = deque()
    positions = []
    last = None
    for i in range(m):
        while dq and kmers[dq[-1]] > kmers[i]:
            dq.pop()
        dq.append(i)
        if dq[0] <= i - w:
            dq.popleft()
        if i >= w - 1:
            mpos = dq[0]
            if mpos != last:
                positions.append(mpos)
                last = mpos
    positions = sorted(set(positions))
    return [positions[i + 1] - positions[i] for i in range(len(positions) - 1)]


def homopolymer_stats(seq, min_run=3):
    n = len(seq)
    i = 0
    hp_bases = 0
    max_run = 0
    while i < n:
        c = seq[i]
        j = i + 1
        while j < n and seq[j] == c:
            j += 1
        run = j - i
        if run >= min_run:
            hp_bases += run
            max_run = max(max_run, run)
        i = j
    return (hp_bases / n if n else 0.0), max_run


def qual_stats(qual_str):
    q = [ord(c) - 33 for c in qual_str]
    n = len(q)
    if n == 0:
        return 0, 0, 0, 0
    head = q[: max(1, n // 5)]
    tail = q[-max(1, n // 5):]
    return (sum(q) / n, min(q), sum(head) / len(head), sum(tail) / len(tail))


def kmer_set(seq, k=K_POOL):
    return {seq[i:i + k] for i in range(len(seq) - k + 1)}


def features_for_pool(barcode, reads):
    """(barcode, [(read_id, seq, qual), ...]) -> list of formatted TSV rows.

    Pure function of its input so it is safe to farm out to a process pool.
    """
    pool_size = len(reads)
    kmer_read_count = Counter()
    per_read_kmers = {}
    for rid, seq, _q in reads:
        ks = kmer_set(seq)
        per_read_kmers[rid] = ks
        for km in ks:
            kmer_read_count[km] += 1

    rows = []
    for rid, seq, qual in reads:
        qm, qmin, qh, qt = qual_stats(qual)
        hp_frac, hp_max = homopolymer_stats(seq)
        gaps = minimizer_gaps(seq)
        if gaps:
            gmean = sum(gaps) / len(gaps)
            gmax = max(gaps)
            gvar = sum((g - gmean) ** 2 for g in gaps) / len(gaps)
        else:
            gmean = gmax = gvar = 0.0

        ks = per_read_kmers[rid]
        # "popular" == this k-mer also appears in at least one OTHER read of
        # the same pool; the single most important feature in the trained
        # model (denovo.md sec 50).
        pop_frac = (sum(1 for km in ks if kmer_read_count[km] >= 2) / len(ks)
                    if ks else 0.0)

        rows.append("\t".join(str(x) for x in [
            rid, barcode, len(seq), f"{qm:.2f}", qmin, f"{qh:.2f}", f"{qt:.2f}",
            f"{qt - qh:.2f}", f"{hp_frac:.4f}", hp_max,
            f"{gmean:.3f}", gmax, f"{gvar:.3f}",
            pool_size, f"{pop_frac:.4f}"]))
    return rows


def _features_one(job):
    barcode, reads = job
    return features_for_pool(barcode, reads)


def iter_fastq_barcode_groups(path):
    """Yield (barcode, [(read_id, seq, qual), ...]) for a BARCODE-SORTED FASTQ.

    Read id is the header up to the first tab -- the sgrep-style header carries
    a trailing "\\tBX:Z:<barcode>" field, and taking the whole line as the id
    silently corrupts every downstream column (a real bug in the prototype).
    Barcode is parsed from the id itself ("...#<barcode>/2") so this does not
    depend on that second field being present.
    """
    opener = gzip.open if path.endswith(".gz") else open
    cur_bc, group = None, []
    with opener(path, "rt") as fh:
        while True:
            header = fh.readline()
            if not header:
                break
            seq = fh.readline().rstrip("\n")
            fh.readline()
            qual = fh.readline().rstrip("\n")
            rid = header.rstrip("\n")[1:].split("\t")[0]
            try:
                bc = rid.split("#")[1].split("/")[0]
            except IndexError:
                continue
            if bc != cur_bc:
                if group:
                    yield cur_bc, group
                cur_bc, group = bc, []
            group.append((rid, seq, qual))
    if group:
        yield cur_bc, group


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--fastq", required=True,
                    help="BARCODE-SORTED FASTQ (.gz ok); see module docstring")
    ap.add_argument("--out", required=True)
    ap.add_argument("--num_processes", type=int, default=1)
    ap.add_argument("--batch-barcodes", type=int, default=20000,
                    help="barcodes held in memory per parallel batch; bounds "
                         "peak memory independently of run size")
    args = ap.parse_args()

    pool = None
    if args.num_processes > 1:
        import multiprocessing as mp
        pool = mp.Pool(args.num_processes)

    n_bc = n_reads = 0
    with open(args.out, "w") as out:
        out.write("\t".join(FIELDS) + "\n")

        def handle_batch(batch):
            nonlocal n_bc, n_reads
            if not batch:
                return
            if pool is not None:
                # chunksize=1: per-barcode cost scales with pool depth, which
                # is heavy-tailed on real data, so pre-assigned chunks stall on
                # whichever worker draws the deep barcodes (same reason
                # denovo_read_filter.py and denovo_junction_qc.py pin it too).
                results = pool.map(_features_one, batch, chunksize=1)
            else:
                results = [_features_one(j) for j in batch]
            for rows in results:
                n_bc += 1
                n_reads += len(rows)
                out.write("\n".join(rows))
                out.write("\n")

        batch = []
        for group in iter_fastq_barcode_groups(args.fastq):
            batch.append(group)
            if len(batch) >= args.batch_barcodes:
                handle_batch(batch)
                batch = []
        handle_batch(batch)

    if pool is not None:
        pool.close()
        pool.join()

    print(f"barcodes={n_bc}", file=sys.stderr)
    print(f"reads={n_reads}", file=sys.stderr)


if __name__ == "__main__":
    main()
