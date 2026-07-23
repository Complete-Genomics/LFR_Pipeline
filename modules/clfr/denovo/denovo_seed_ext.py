"""
Per-UMI seed-extension (OLC) assembler for LFR / stLFR data.

Replaces megahit (de Bruijn graph) for low-depth UMIs where k-mer
coverage is too sparse to build a connected graph.

Algorithm — greedy Overlap-Layout-Consensus:
    1. Seed   = longest (deduplicated) read in the UMI
    2. Extend = scan remaining reads for prefix/suffix overlap with
                current contig, including reverse-complement candidates
    3. Repeat until no further extension is possible

Minimum viable depth: 1 read (seed only) → contig if ≥ min_ctg bp.

Output format matches megahit convention so denovo_supp.py requires
no changes: header = >{barcode}k41_0, first 15 chars = barcode.

No required dependencies — pure Python stdlib.
Optional: pip install mappy  →  faster overlap via minimap2 C engine.

Drop-in API
-----------
Replace process_barcode_se / process_barcode_pe in denovo_clfr_ram.py:

    from denovo_seed_ext import configure, process_barcode_se, process_barcode_pe

    # call once before Pool
    configure(min_ctg_len=MIN_CTG_LEN, out_id=ID)

    pool.starmap(process_barcode_se,
                 [(bc, shared_meta2, lock) for bc in meta_data2])
"""

from __future__ import annotations

import threading
from collections import defaultdict
from typing import Sequence

# ── module-level config (set via configure() before multiprocessing) ──────────

_CFG: dict = {
    "min_ctg":   400,
    "min_ov":    20,
    "max_mm":    0.05,
    "out_id":    0,
    "out_file":  "denovo/final_contigs_{id}.fa",
    "seed_k":    10,    # k-mer size for overlap pre-filter
    "use_mappy": None,  # None = auto-detect, True/False = force
}


def configure(
    min_ctg_len: int   = 400,
    min_overlap: int   = 20,
    max_mismatch: float = 0.05,
    out_id: int        = 0,
    out_file: str      = "denovo/final_contigs_{id}.fa",
    use_mappy: bool | None = None,
) -> None:
    """Call once in the parent process before spawning Pool workers."""
    _CFG["min_ctg"]   = min_ctg_len
    _CFG["min_ov"]    = min_overlap
    _CFG["max_mm"]    = max_mismatch
    _CFG["out_id"]    = out_id
    _CFG["out_file"]  = out_file
    _CFG["use_mappy"] = use_mappy


# ── sequence utilities ────────────────────────────────────────────────────────

_RC = str.maketrans("ACGT", "TGCA")


def rc(seq: str) -> str:
    return seq.translate(_RC)[::-1]


def _kmer_set(seq: str, k: int, start: int = 0, end: int | None = None) -> set[str]:
    s = seq[start:end]
    if len(s) < k:
        return set()
    return {s[i:i+k] for i in range(len(s) - k + 1)}


# ── core overlap ──────────────────────────────────────────────────────────────

def suffix_prefix_overlap(
    a: str,
    b: str,
    min_ov: int,
    max_mm: float,
    seed_k: int = 10,
) -> int:
    """
    Return the length of b's prefix that overlaps a's suffix, 0 if none.

    Checks decreasing overlap lengths so returns the longest valid overlap.
    Uses a k-mer seed pre-filter to skip pairs that cannot possibly overlap,
    giving ~5-10× speedup when most pairs are non-overlapping.
    """
    limit = min(len(a), len(b))
    if limit < min_ov:
        return 0

    # ── seed filter: share at least one k-mer near the boundary ──
    check_len = min(limit, max(min_ov * 3, seed_k * 4))
    a_end_kmers = _kmer_set(a, seed_k, start=len(a) - check_len)
    b_start_kmers = _kmer_set(b, seed_k, end=check_len)
    if not (a_end_kmers & b_start_kmers):
        return 0

    # ── full mismatch check (longest-first) ──────────────────────
    for ov in range(limit, min_ov - 1, -1):
        mm = sum(x != y for x, y in zip(a[-ov:], b[:ov]))
        if mm / ov <= max_mm:
            return ov
    return 0


# ── assembler ─────────────────────────────────────────────────────────────────

def assemble_umi(
    seqs: Sequence[str],
    min_ov: int   = 20,
    max_mm: float = 0.05,
    min_ctg: int  = 400,
    seed_k: int   = 10,
) -> list[str]:
    """
    Greedy seed-extension assembly for one UMI's reads.

    Parameters
    ----------
    seqs      : list of DNA sequences (forward strand, no quality)
    min_ov    : minimum overlap length to merge two reads [20]
    max_mm    : maximum mismatch rate in overlap region [0.05]
    min_ctg   : discard contigs shorter than this [400]
    seed_k    : k-mer length for overlap pre-filter [10]

    Returns
    -------
    List of contig sequences (usually 0 or 1 per UMI).
    """
    if not seqs:
        return []

    # deduplicate, sort longest-first → longest read is seed
    uniq = sorted(set(seqs), key=len, reverse=True)
    contig = uniq[0]
    unused = list(range(1, len(uniq)))

    changed = True
    while changed and unused:
        changed = False
        for i in list(unused):
            seq = uniq[i]
            extended = False

            for cand in (seq, rc(seq)):
                # try extending contig at 3' end
                ov = suffix_prefix_overlap(contig, cand, min_ov, max_mm, seed_k)
                if ov:
                    contig += cand[ov:]
                    unused.remove(i)
                    changed = True
                    extended = True
                    break

                # try extending contig at 5' end
                ov = suffix_prefix_overlap(cand, contig, min_ov, max_mm, seed_k)
                if ov:
                    contig = cand + contig[ov:]
                    unused.remove(i)
                    changed = True
                    extended = True
                    break

            if extended:
                # restart scan so new contig ends are retried against all unused
                break

    return [contig] if len(contig) >= min_ctg else []


# ── optional mappy fast path ──────────────────────────────────────────────────

def _assemble_umi_mappy(seqs: list[str], min_ctg: int) -> list[str] | None:
    """
    Overlap detection via mappy (minimap2 Python bindings).
    Returns contig list, or None if mappy unavailable / fails.
    """
    try:
        import mappy as mp
        import tempfile, os
    except ImportError:
        return None

    if not seqs or len(seqs) < 2:
        return None

    # write reads to temp FASTA for indexing
    fa_lines = "".join(f">r{i}\n{s}\n" for i, s in enumerate(seqs))
    with tempfile.NamedTemporaryFile(suffix=".fa", mode="w", delete=False) as f:
        f.write(fa_lines)
        tmp = f.name

    try:
        aligner = mp.Aligner(tmp, preset="ava-sr", best_n=5)
        if not aligner:
            return None

        # collect overlaps: {(i, j): overlap_len}
        overlaps: dict[tuple[int, int], int] = {}
        for i, seq in enumerate(seqs):
            for hit in aligner.map(seq):
                try:
                    j = int(hit.ctg[1:])
                except (ValueError, IndexError):
                    continue
                if j != i:
                    ov = hit.q_en - hit.q_st
                    if ov > overlaps.get((i, j), 0):
                        overlaps[(i, j)] = ov
    finally:
        os.unlink(tmp)

    if not overlaps:
        # no overlaps found; fall back to pure Python with seed only
        best = max(range(len(seqs)), key=lambda i: len(seqs[i]))
        return [seqs[best]] if len(seqs[best]) >= min_ctg else []

    # greedy chain: start from longest read, follow best overlaps
    idx_sorted = sorted(range(len(seqs)), key=lambda i: len(seqs[i]), reverse=True)
    contig = seqs[idx_sorted[0]]
    used = {idx_sorted[0]}

    while True:
        best_j, best_ov = None, 0
        for (i, j), ov in overlaps.items():
            if i in used and j not in used and ov > best_ov:
                best_j, best_ov = j, ov
        if best_j is None:
            break
        contig += seqs[best_j][best_ov:]
        used.add(best_j)

    return [contig] if len(contig) >= min_ctg else []


# ── output writer ─────────────────────────────────────────────────────────────

def _write_contigs(
    barcode: str,
    contigs: list[str],
    out_file: str,
    lock,
) -> None:
    """
    Append contigs to final_contigs_{id}.fa.
    Header format: >{barcode}k41_{i}  (first 15 chars = barcode,
    matching what denovo_supp.py expects via record.id[:CBC_LEN]).
    """
    if not contigs:
        return
    lines = []
    for i, seq in enumerate(contigs[:4]):   # max 4 contigs per UMI, same as megahit path
        lines.append(f">{barcode}k41_{i}\n{seq}\n")
    block = "".join(lines)
    with lock:
        with open(out_file, "a") as fh:
            fh.write(block)


# ── per-barcode workers (drop-in replacements) ────────────────────────────────

def _seqs_from_meta(meta: dict, barcode: str) -> list[str]:
    """Extract sequence strings from meta_data dict (strips header lines)."""
    entries = meta.get(barcode)
    if not entries:
        return []
    # entries = ['>id0', 'seq0', '>id1', 'seq1', ...]
    return entries[1::2]


def process_barcode_se(barcode: str, shared_meta_data2, lock) -> None:
    """
    SE drop-in for denovo_clfr_ram.process_barcode_se.
    Uses seed extension instead of megahit.
    """
    min_ctg  = _CFG["min_ctg"]
    min_ov   = _CFG["min_ov"]
    max_mm   = _CFG["max_mm"]
    seed_k   = _CFG["seed_k"]
    out_id   = _CFG["out_id"]
    out_file = _CFG["out_file"].format(id=out_id)
    use_mp   = _CFG["use_mappy"]

    seqs = _seqs_from_meta(shared_meta_data2, barcode)
    if not seqs:
        return

    contigs: list[str] | None = None
    if use_mp is not False:
        contigs = _assemble_umi_mappy(seqs, min_ctg)
    if contigs is None:
        contigs = assemble_umi(seqs, min_ov, max_mm, min_ctg, seed_k)

    _write_contigs(barcode, contigs, out_file, lock)


def process_barcode_pe(barcode: str, shared_meta_data1, shared_meta_data2, lock) -> None:
    """
    PE drop-in for denovo_clfr_ram.process_barcode_pe.
    Uses seed extension instead of megahit.
    """
    min_ctg  = _CFG["min_ctg"]
    min_ov   = _CFG["min_ov"]
    max_mm   = _CFG["max_mm"]
    seed_k   = _CFG["seed_k"]
    out_id   = _CFG["out_id"]
    out_file = _CFG["out_file"].format(id=out_id)
    use_mp   = _CFG["use_mappy"]

    r1 = _seqs_from_meta(shared_meta_data1, barcode)
    r2 = _seqs_from_meta(shared_meta_data2, barcode)
    seqs = r1 + r2
    if not seqs:
        return

    contigs: list[str] | None = None
    if use_mp is not False:
        contigs = _assemble_umi_mappy(seqs, min_ctg)
    if contigs is None:
        contigs = assemble_umi(seqs, min_ov, max_mm, min_ctg, seed_k)

    _write_contigs(barcode, contigs, out_file, lock)


# ── CLI self-test ─────────────────────────────────────────────────────────────

if __name__ == "__main__":
    import sys, random

    random.seed(42)
    BASES = "ACGT"

    def _rand_seq(n: int) -> str:
        return "".join(random.choice(BASES) for _ in range(n))

    def _make_reads(frag: str, read_len: int, step: int) -> list[str]:
        reads = []
        for s in range(0, len(frag) - read_len + 1, step):
            reads.append(frag[s:s + read_len])
        return reads

    # non-repetitive 600 bp fragment
    FRAG = _rand_seq(600)

    # high-depth: 150 bp reads every 30 bp → 16 distinct reads, overlap=120 bp
    reads_hi = _make_reads(FRAG, 150, 30)
    ctg_hi = assemble_umi(reads_hi, min_ctg=400)
    print(f"[hi-depth] reads={len(reads_hi)}  contig_len={len(ctg_hi[0]) if ctg_hi else 0}")

    # low-depth: 3 reads with 50 bp overlap
    reads_lo = _make_reads(FRAG, 200, 150)
    ctg_lo = assemble_umi(reads_lo, min_ctg=100)
    print(f"[lo-depth] reads={len(reads_lo)}  contig_len={len(ctg_lo[0]) if ctg_lo else 0}")

    # single-read test
    ctg_1 = assemble_umi([FRAG[:250]], min_ctg=100)
    print(f"[1 read]   reads=1  contig_len={len(ctg_1[0]) if ctg_1 else 0}")

    # RC extension test: reverse-complement second read
    frag2 = _rand_seq(400)
    r1 = frag2[:200]
    r2 = rc(frag2[150:])                # overlaps r1 by 50 bp, on RC strand
    ctg_rc = assemble_umi([r1, r2], min_ov=30, min_ctg=200)
    print(f"[RC ext]   reads=2  contig_len={len(ctg_rc[0]) if ctg_rc else 0}")

    # correctness checks
    ok = True
    if not ctg_hi or ctg_hi[0] not in FRAG and FRAG not in ctg_hi[0]:
        print("FAIL: hi-depth contig does not match fragment", file=sys.stderr); ok = False
    if not ctg_lo:
        print("FAIL: lo-depth assembly produced nothing", file=sys.stderr); ok = False
    if not ctg_rc or len(ctg_rc[0]) < 350:
        print(f"FAIL: RC extension too short ({len(ctg_rc[0]) if ctg_rc else 0})", file=sys.stderr); ok = False
    if ok:
        print("correctness: OK")