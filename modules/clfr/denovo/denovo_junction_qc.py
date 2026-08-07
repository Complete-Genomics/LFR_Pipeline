#!/usr/bin/env python3
"""Reference-free chimera detection via spanning-read depth discontinuity.

Why this works where read-back QC does not (denovo.md sec 30): a chimera
built from a genuinely mixed read pool IS well covered by that same pool --
every read belongs to the barcode, so plain coverage/breadth looks fine.
But coverage is not the same as *spanning* support. At a chimera junction,
reads to the left come from source molecule A and reads to the right from
source molecule B, and no physical molecule contains both sides, so no read
can genuinely span the junction with real margin on both sides. The
assembler joined the two halves through a short conserved-region overlap
(denovo.md sec 28), which is exactly why the join exists at all -- but that
short overlap cannot manufacture deeply-spanning reads.

So: place each contig's own reads back onto it, and for every interior
position count reads that cross it with >= margin bp on BOTH sides. A true
single-molecule contig has a roughly uniform spanning profile; a chimera has
a sharp dip (often to zero) at the junction while ordinary coverage stays
healthy on both flanks.
"""
import argparse
import csv
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + "/test")
from analyze_olc_shortfall import kmer_index, place, rc  # noqa: E402
from compare_olc_readback import load_fasta  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from denovo_read_filter import iter_barcode_groups  # noqa: E402


def place_oriented(read, contig, index):
    """Like analyze_olc_shortfall.place_either, but also returns the oriented
    read, since scoring identity requires knowing which strand was placed."""
    hits = []
    for oriented in (read, rc(read)):
        result = place(oriented, contig, index)
        if result is not None:
            hits.append((result, oriented))
    if len(hits) != 1:
        return None
    (start, end, _, _), oriented = hits[0]
    return start, end, oriented


def _match_profile(read, contig, start):
    """Per-contig-position match flags for an ungapped placement at `start`.
    Returns (lo, hi, matches) where matches[i] is 1/0 for contig[lo+i]."""
    n = len(contig)
    lo, hi = max(0, start), min(n, start + len(read))
    matches = []
    for pos in range(lo, hi):
        ri = pos - start
        matches.append(1 if read[ri] == contig[pos] else 0)
    return lo, hi, matches


def spanning_profile(contig, reads, margin, min_side_identity, window):
    """Return (spanning_depth[], plain_depth[], placed).

    A read only counts as spanning position p if it extends >= margin bp past
    p on both sides AND actually matches the contig on both sides at
    >= min_side_identity. That identity check is the whole point: placement
    alone is far too permissive here, because 16S conserved regions let a read
    from species A anchor onto (and appear to span) a stretch that actually
    belongs to species B (denovo.md sec 31). Requiring real two-sided identity
    is what makes a chimera junction show up as a spanning-depth collapse
    while ordinary coverage stays healthy.
    """
    index = kmer_index(contig, 17)
    n = len(contig)
    span_diff = [0] * (n + 1)
    plain_diff = [0] * (n + 1)
    placed = 0
    for read in set(reads):
        hit = place_oriented(read, contig, index)
        if hit is None:
            continue
        start, _end, oriented = hit
        lo, hi, matches = _match_profile(oriented, contig, start)
        if lo >= hi:
            continue
        placed += 1
        plain_diff[lo] += 1
        plain_diff[hi] -= 1

        s_lo, s_hi = lo + margin, hi - margin
        if s_lo >= s_hi:
            continue
        # prefix[i] = matches in matches[:i]; lets us score any window in O(1)
        prefix = [0] * (len(matches) + 1)
        for i, m in enumerate(matches):
            prefix[i + 1] = prefix[i] + m
        m_len = len(matches)
        for p in range(s_lo, s_hi):
            k = p - lo
            # Score a local window on each side rather than the whole side:
            # placement here is ungapped (single k-mer offset), so one indel
            # anywhere would frameshift and wreck whole-side identity even for
            # a perfectly good read. A window local to p keeps the test about
            # "does this read really match across THIS position".
            l_start = max(0, k - window)
            r_end = min(m_len, k + window)
            l_len, r_len = k - l_start, r_end - k
            if l_len < margin or r_len < margin:
                continue
            left_id = (prefix[k] - prefix[l_start]) / l_len
            right_id = (prefix[r_end] - prefix[k]) / r_len
            if left_id >= min_side_identity and right_id >= min_side_identity:
                span_diff[p] += 1
                span_diff[p + 1] -= 1

    span, plain = [], []
    sv = pv = 0
    for i in range(n):
        sv += span_diff[i]
        pv += plain_diff[i]
        span.append(sv)
        plain.append(pv)
    return span, plain, placed


def analyze(contig, reads, margin, edge_skip, min_side_identity, window, smooth):
    span, plain, placed = spanning_profile(contig, reads, margin, min_side_identity,
                                            window)
    n = len(contig)
    lo, hi = edge_skip, n - edge_skip
    if hi <= lo or placed == 0:
        return None

    interior_span = span[lo:hi]
    interior_plain = plain[lo:hi]
    min_span = min(interior_span)
    min_pos = lo + interior_span.index(min_span)
    # plain coverage at the weakest spanning point -- a real junction shows
    # healthy ordinary coverage but no spanning support
    plain_at_min = plain[min_pos]
    median_span = sorted(interior_span)[len(interior_span) // 2]
    median_plain = sorted(interior_plain)[len(interior_plain) // 2]

    # Normalising spanning depth by local ordinary coverage is what makes the
    # metric comparable across contigs: a low-depth UMI has low spanning depth
    # everywhere for innocent reasons, whereas a junction is specifically a
    # place where coverage is fine but verified spanning support is not.
    # Smoothing first keeps a single noisy position from dominating (the raw
    # per-position minimum turned out to be far too noisy to threshold on).
    half = max(1, smooth // 2)
    m = len(interior_span)
    ps_span = [0] * (m + 1)
    ps_plain = [0] * (m + 1)
    for i in range(m):
        ps_span[i + 1] = ps_span[i] + interior_span[i]
        ps_plain[i + 1] = ps_plain[i] + interior_plain[i]
    min_ratio, min_ratio_pos = 1.0, -1
    found = False
    for i in range(m):
        a, b = max(0, i - half), min(m, i + half + 1)
        p_sum = ps_plain[b] - ps_plain[a]
        if p_sum <= 0:
            continue
        r = (ps_span[b] - ps_span[a]) / p_sum
        if not found or r < min_ratio:
            min_ratio, min_ratio_pos, found = r, lo + i, True

    return {
        "contig_len": n,
        "placed_reads": placed,
        "min_spanning_depth": min_span,
        "min_spanning_pos": min_pos,
        "plain_depth_at_min": plain_at_min,
        "median_spanning_depth": median_span,
        "median_plain_depth": median_plain,
        "span_cov_ratio": (median_span / median_plain) if median_plain else 1.0,
        "min_local_span_ratio": round(min_ratio, 4),
        "min_local_span_pos": min_ratio_pos,
    }


_CFG = {}
_CONTIGS = {}


def _init_worker(cfg, contigs):
    _CFG.update(cfg)
    _CONTIGS.update(contigs)


def _analyze_one(job):
    """(barcode, reads) -> (barcode, result dict or None). Contigs come from
    the worker-global map so the (large) contig set is shipped once at pool
    startup rather than once per barcode."""
    barcode, reads = job
    seqs = _CONTIGS.get(barcode)
    if not seqs:
        return barcode, None
    contig = max(seqs, key=len)
    return barcode, analyze(contig, reads, _CFG["margin"], _CFG["edge_skip"],
                             _CFG["min_side_identity"], _CFG["window"],
                             _CFG["smooth"])


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--r2", required=True)
    ap.add_argument("--contigs", required=True)
    ap.add_argument("--n", type=int, default=100000000)
    ap.add_argument("--margin", type=int, default=30,
                     help="bp a read must extend past a position on both sides "
                          "to count as spanning it")
    ap.add_argument("--edge-skip", type=int, default=100,
                     help="ignore this many bp at each contig end (spanning "
                          "depth is legitimately low there)")
    ap.add_argument("--min-side-identity", type=float, default=0.95,
                     help="a read must match the contig at least this well on "
                          "BOTH sides of a position to count as spanning it")
    ap.add_argument("--window", type=int, default=100,
                     help="bp of local context scored on each side of a position")
    ap.add_argument("--smooth", type=int, default=101,
                     help="bp smoothing window for the local span/coverage ratio")
    ap.add_argument("--max-span-ratio", type=float, default=0.25,
                     help="flag contig if median verified-spanning depth divided "
                          "by median ordinary coverage falls below this. Validated "
                          "on the ZymoBIOMICS control (denovo.md sec 31): AUC 0.827; "
                          "at 0.25 this recovers ~54%% of confident chimeras at ~3x "
                          "enrichment, dropping ~14%% of clean contigs with it")
    ap.add_argument("--min-local-span-ratio", type=float, default=0.0,
                     help="optional high-purity gate: also require the weakest "
                          "smoothed local span/coverage ratio to reach this. "
                          "0 disables. Ranks worse than span_cov_ratio overall "
                          "(AUC 0.66 vs 0.83) but its extreme tail is very pure: "
                          "on the Zymo control, 0.20 keeps 10%% of contigs at 97.8%% "
                          "mean identity vs 95.4%% for span_cov_ratio alone")
    ap.add_argument("--num_processes", type=int, default=1)
    ap.add_argument("--batch-barcodes", type=int, default=20000,
                     help="barcodes held in memory per parallel batch; the read "
                          "set is streamed rather than fully loaded, which at "
                          "1.5-3M UMI is the difference between a few GB and "
                          "tens of GB of RAM")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    contigs = load_fasta(args.contigs)
    cfg = {"margin": args.margin, "edge_skip": args.edge_skip,
           "min_side_identity": args.min_side_identity,
           "window": args.window, "smooth": args.smooth}

    fields = ["barcode", "contig_len", "placed_reads", "min_spanning_depth",
              "min_spanning_pos", "plain_depth_at_min", "median_spanning_depth",
              "median_plain_depth", "span_cov_ratio", "min_local_span_ratio",
              "min_local_span_pos", "junction_suspect"]
    pool = None
    if args.num_processes > 1:
        import multiprocessing as mp
        pool = mp.Pool(args.num_processes, initializer=_init_worker,
                        initargs=(cfg, contigs))
    else:
        _init_worker(cfg, contigs)

    n_flag = 0
    n_total = 0
    with open(args.out, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        writer.writeheader()

        def run_batch(batch):
            nonlocal n_flag, n_total
            if not batch:
                return
            # chunksize=1: cost per barcode scales with contig length x read
            # count, so a deep barcode would stall a whole pre-assigned chunk.
            if pool is not None:
                results = pool.map(_analyze_one, batch, chunksize=1)
            else:
                results = [_analyze_one(j) for j in batch]
            for barcode, res in results:
                if res is None:
                    continue
                n_total += 1
                suspect = int(res["span_cov_ratio"] < args.max_span_ratio
                                or res["min_local_span_ratio"] < args.min_local_span_ratio)
                n_flag += suspect
                row = {"barcode": barcode, "junction_suspect": suspect}
                row.update(res)
                writer.writerow(row)

        batch = []
        seen = 0
        for barcode, _lines, seqs in iter_barcode_groups(args.r2):
            seen += 1
            if seen > args.n:
                break
            batch.append((barcode, seqs))
            if len(batch) >= args.batch_barcodes:
                run_batch(batch)
                batch = []
        run_batch(batch)

    if pool is not None:
        pool.close()
        pool.join()

    print(f"contigs_checked={n_total}")
    if n_total:
        print(f"junction_suspect={n_flag} ({100*n_flag/n_total:.2f}%)")


if __name__ == "__main__":
    main()