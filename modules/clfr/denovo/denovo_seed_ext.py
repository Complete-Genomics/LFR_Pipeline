"""
Per-UMI seed-extension (OLC) assembler for LFR / stLFR data.

Replaces megahit (de Bruijn graph) for low-depth UMIs where k-mer
coverage is too sparse to build a connected graph.

Algorithm - greedy Overlap-Layout-Consensus:
    1. Seed   = longest (deduplicated) read in the UMI
    2. Extend = scan remaining reads for prefix/suffix overlap with
                current contig, including reverse-complement candidates
    3. Repeat until no further extension is possible

Minimum viable depth: 1 read (seed only) -> contig if >= min_ctg bp.

Output format matches megahit convention so denovo_supp.py requires
no changes: header = >{barcode}>k41_0, first 15 chars = barcode.

No required dependencies - pure Python stdlib (Python >= 3.6).
Optional: pip install mappy  ->  faster overlap via minimap2 C engine.

Self-test (no args)
-------------------
    python3 denovo_seed_ext.py

Standalone CLI (drop-in replacement for denovo_clfr_ram.py --module denovo_parallel)
-------------------------------------------------------------------------------------
Run directly from a Snakemake work dir containing denovo/data_R1_sgrep.tsv
and denovo/data_R2_sgrep.tsv (same layout denovo_clfr_ram.py expects).
No megahit binary, no tmp_dir, no subprocess fork per UMI.

    python3 denovo_seed_ext.py \\
        --sequence_type se \\
        --num_processes 30 \\
        --n_line_chunk 2000000 \\
        --min_ctg_len 400 \\
        --nth_of_nodes 0 \\
        --n 1000            # optional: only assemble first 1000 UMIs (config: assembly_N_umi); omit/empty = all UMIs

Writes contigs to denovo/final_contigs_{nth_of_nodes}.fa and touches
denovo/frag_denovo_done, matching denovo_clfr_ram.py's output contract
so downstream rules (map_denovo, correc_direction_denovo, ...) need no changes.

See denovo_clfr.smk rule run_denovo_parallel for the branch that invokes
this script directly when frag_de_novo.assembler == 'seedext'.

Benchmark only, no pipeline side effects (run from Snakemake work dir)
------------------------------------------------------------------------
    python3 /path/to/benchmark_seedext.py \\
        --n 1000 \\
        --r2 denovo/data_R2_sgrep.tsv \\
        --min_ctg 400

Output printed:
    per-UMI latency, throughput (UMI/s), contig yield, 1M/3M extrapolation
"""

import itertools
import os
from collections import defaultdict

# ── module-level config (set via configure() before multiprocessing) ──────────

_CFG = {
    "min_ctg":   400,
    "min_ov":    20,
    "max_mm":    0.05,
    "out_id":    0,
    "out_file":  "denovo/final_contigs_{id}.fa",
    "seed_k":    10,    # k-mer size for overlap pre-filter
    "use_mappy": None,  # None = auto-detect, True/False = force
}


def configure(min_ctg_len=400, min_overlap=20, max_mismatch=0.05,
              out_id=0, out_file="denovo/final_contigs_{id}.fa", use_mappy=None):
    """Call once in the parent process before spawning Pool workers."""
    _CFG["min_ctg"]   = min_ctg_len
    _CFG["min_ov"]    = min_overlap
    _CFG["max_mm"]    = max_mismatch
    _CFG["out_id"]    = out_id
    _CFG["out_file"]  = out_file
    _CFG["use_mappy"] = use_mappy


# ── sequence utilities ────────────────────────────────────────────────────────

_RC = str.maketrans("ACGT", "TGCA")


def rc(seq):
    return seq.translate(_RC)[::-1]


def _kmer_set(seq, k, start=0, end=None):
    s = seq[start:end]
    if len(s) < k:
        return set()
    return {s[i:i+k] for i in range(len(s) - k + 1)}


# ── core overlap ──────────────────────────────────────────────────────────────

def suffix_prefix_overlap(a, b, min_ov, max_mm, seed_k=10):
    """
    Return the length of b's prefix that overlaps a's suffix, 0 if none.

    Checks decreasing overlap lengths so returns the longest valid overlap.
    Uses a k-mer seed pre-filter to skip pairs that cannot possibly overlap,
    giving ~5-10x speedup when most pairs are non-overlapping.
    """
    limit = min(len(a), len(b))
    if limit < min_ov:
        return 0

    # seed filter: share at least one k-mer near the boundary
    check_len = min(limit, max(min_ov * 3, seed_k * 4))
    a_end_kmers = _kmer_set(a, seed_k, start=len(a) - check_len)
    b_start_kmers = _kmer_set(b, seed_k, end=check_len)
    if not (a_end_kmers & b_start_kmers):
        return 0

    # full mismatch check (longest-first)
    for ov in range(limit, min_ov - 1, -1):
        mm = sum(x != y for x, y in zip(a[-ov:], b[:ov]))
        if mm / ov <= max_mm:
            return ov
    return 0


# ── assembler ─────────────────────────────────────────────────────────────────

def assemble_umi(seqs, min_ov=20, max_mm=0.05, min_ctg=400, seed_k=10):
    """
    Greedy seed-extension assembly for one UMI's reads.

    seqs    : list of DNA sequences (forward strand, no quality)
    min_ov  : minimum overlap length to merge two reads [20]
    max_mm  : maximum mismatch rate in overlap region [0.05]
    min_ctg : discard contigs shorter than this [400]
    seed_k  : k-mer length for overlap pre-filter [10]

    Returns list of contig sequences (usually 0 or 1 per UMI).
    """
    if not seqs:
        return []

    # deduplicate, sort longest-first -> longest read is seed
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

def _assemble_umi_mappy(seqs, min_ctg):
    """
    Overlap detection via mappy (minimap2 Python bindings).
    Returns contig list, or None if mappy unavailable / fails.
    """
    try:
        import mappy as mp
        import tempfile
        import os
    except ImportError:
        return None

    if not seqs or len(seqs) < 2:
        return None

    fa_lines = "".join(">r{}\n{}\n".format(i, s) for i, s in enumerate(seqs))
    with tempfile.NamedTemporaryFile(suffix=".fa", mode="w", delete=False) as f:
        f.write(fa_lines)
        tmp = f.name

    try:
        aligner = mp.Aligner(tmp, preset="ava-sr", best_n=5)
        if not aligner:
            return None

        overlaps = {}
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
        best = max(range(len(seqs)), key=lambda i: len(seqs[i]))
        return [seqs[best]] if len(seqs[best]) >= min_ctg else []

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

def _write_contigs(barcode, contigs, out_file, lock):
    """
    Append contigs to final_contigs_{id}.fa.
    Header: >{barcode}>k41_{i}  (first 15 chars = barcode,
    matching denovo_supp.py record.id[:CBC_LEN]; second '>' marks the
    barcode/UMI boundary, matching megahit-path convention).
    """
    if not contigs:
        return
    lines = []
    for i, seq in enumerate(contigs[:4]):   # max 4 per UMI, same as megahit path
        lines.append(">{barcode}>k41_{i}\n{seq}\n".format(barcode=barcode, i=i, seq=seq))
    block = "".join(lines)
    with lock:
        with open(out_file, "a") as fh:
            fh.write(block)


# ── per-barcode workers (drop-in replacements) ────────────────────────────────

def _seqs_from_meta(meta, barcode):
    """Extract sequence strings from meta_data dict (strips header lines)."""
    entries = meta.get(barcode)
    if not entries:
        return []
    # entries = ['>id0', 'seq0', '>id1', 'seq1', ...]
    return entries[1::2]


def process_barcode_se(barcode, shared_meta_data2, lock):
    """SE drop-in for denovo_clfr_ram.process_barcode_se."""
    min_ctg  = _CFG["min_ctg"]
    min_ov   = _CFG["min_ov"]
    max_mm   = _CFG["max_mm"]
    seed_k   = _CFG["seed_k"]
    out_file = _CFG["out_file"].format(id=_CFG["out_id"])
    use_mp   = _CFG["use_mappy"]

    seqs = _seqs_from_meta(shared_meta_data2, barcode)
    if not seqs:
        return

    contigs = None
    if use_mp is not False:
        contigs = _assemble_umi_mappy(seqs, min_ctg)
    if contigs is None:
        contigs = assemble_umi(seqs, min_ov, max_mm, min_ctg, seed_k)

    _write_contigs(barcode, contigs, out_file, lock)


def process_barcode_pe(barcode, shared_meta_data1, shared_meta_data2, lock):
    """PE drop-in for denovo_clfr_ram.process_barcode_pe."""
    min_ctg  = _CFG["min_ctg"]
    min_ov   = _CFG["min_ov"]
    max_mm   = _CFG["max_mm"]
    seed_k   = _CFG["seed_k"]
    out_file = _CFG["out_file"].format(id=_CFG["out_id"])
    use_mp   = _CFG["use_mappy"]

    r1 = _seqs_from_meta(shared_meta_data1, barcode)
    r2 = _seqs_from_meta(shared_meta_data2, barcode)
    seqs = r1 + r2
    if not seqs:
        return

    contigs = None
    if use_mp is not False:
        contigs = _assemble_umi_mappy(seqs, min_ctg)
    if contigs is None:
        contigs = assemble_umi(seqs, min_ov, max_mm, min_ctg, seed_k)

    _write_contigs(barcode, contigs, out_file, lock)


# ── sgrep TSV parsing (same format as denovo_clfr_ram.add_sgrep_line) ──────────

def _add_sgrep_line(meta_data, line):
    """
    Parse one line of denovo/data_R{1,2}_sgrep.tsv into meta_data[barcode].
    Line format: <readname>\\t<seq>, readname[5:20] = 15-char barcode.
    Appends '>id' then 'seq' so meta_data[bc] = ['>id0','seq0','>id1','seq1',...],
    matching denovo_clfr_ram.py's convention exactly.
    """
    bc_len = 15
    info = line.rstrip("\n").split("\t")
    if len(info) < 2:
        return False
    bc = info[0][5:5 + bc_len]
    rid = ">" + info[0][22:]
    seq = info[1]
    meta_data[bc].append(rid)
    meta_data[bc].append(seq)
    return True


def _iter_se_chunks(r2_path, start_idx, n_line_chunk):
    with open(r2_path) as f:
        for _ in itertools.islice(f, start_idx):
            pass
        chunk_start = start_idx
        while True:
            meta_data2 = defaultdict(list)
            n_lines = 0
            for line in itertools.islice(f, n_line_chunk):
                if _add_sgrep_line(meta_data2, line):
                    n_lines += 1
            if n_lines == 0:
                break
            yield chunk_start, meta_data2
            chunk_start += n_lines


def _iter_pe_chunks(r1_path, r2_path, start_idx, n_line_chunk):
    with open(r1_path) as f1, open(r2_path) as f2:
        for _ in itertools.islice(f1, start_idx):
            pass
        for _ in itertools.islice(f2, start_idx):
            pass
        chunk_start = start_idx
        while True:
            meta_data1 = defaultdict(list)
            meta_data2 = defaultdict(list)
            n_lines = 0
            for line1, line2 in itertools.islice(zip(f1, f2), n_line_chunk):
                ok1 = _add_sgrep_line(meta_data1, line1)
                ok2 = _add_sgrep_line(meta_data2, line2)
                if ok1 and ok2:
                    n_lines += 1
            if n_lines == 0:
                break
            yield chunk_start, meta_data1, meta_data2
            chunk_start += n_lines


def _create_bins(start_idx, end_idx, bin_size):
    bins_ = []
    for i in range(start_idx, end_idx, bin_size):
        bins_.append((i, min(i + bin_size, end_idx)))
    return bins_


def _limit_umis(meta_data, remaining):
    """
    Keep at most `remaining` barcodes from meta_data (dict insertion order).
    Returns (possibly-truncated meta_data, number of barcodes kept).
    remaining=None means no limit.
    """
    if remaining is None or len(meta_data) <= remaining:
        return meta_data, len(meta_data)
    limited = defaultdict(list)
    for i, bc in enumerate(meta_data.keys()):
        if i >= remaining:
            break
        limited[bc] = meta_data[bc]
    return limited, len(limited)


class _NullLock(object):
    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        return False


def _process_pe_metadata(meta_data1, meta_data2, num_processes):
    if num_processes == 1:
        lock = _NullLock()
        for barcode in meta_data2.keys():
            process_barcode_pe(barcode, meta_data1, meta_data2, lock)
    else:
        import multiprocessing as mp
        with mp.Manager() as manager:
            shared1 = manager.dict(meta_data1)
            shared2 = manager.dict(meta_data2)
            lock = manager.Lock()
            with mp.Pool(num_processes) as pool:
                pool.starmap(process_barcode_pe,
                             [(bc, shared1, shared2, lock) for bc in meta_data2.keys()])
    print("denovo_BC_counts={}".format(len(meta_data2)))
    return sum(len(v) // 2 for v in meta_data2.values())


def _process_se_metadata(meta_data2, num_processes):
    if num_processes == 1:
        lock = _NullLock()
        for barcode in meta_data2.keys():
            process_barcode_se(barcode, meta_data2, lock)
    else:
        import multiprocessing as mp
        with mp.Manager() as manager:
            shared2 = manager.dict(meta_data2)
            lock = manager.Lock()
            with mp.Pool(num_processes) as pool:
                pool.starmap(process_barcode_se,
                             [(bc, shared2, lock) for bc in meta_data2.keys()])
    print("denovo_BC_counts={}".format(len(meta_data2)))
    return sum(len(v) // 2 for v in meta_data2.values())


# ── standalone pipeline CLI ───────────────────────────────────────────────────

def _main_cli():
    import argparse
    import datetime
    import subprocess

    ap = argparse.ArgumentParser(
        description="Standalone per-UMI seed-extension assembler (no megahit, no subprocess fork)")
    ap.add_argument("--sequence_type", choices=["se", "pe"], required=True)
    ap.add_argument("--num_processes", type=int, default=1)
    ap.add_argument("--n_line_chunk", type=int, default=2000000)
    ap.add_argument("--start_idx", type=int, default=0)
    ap.add_argument("--end_idx", type=int, default=None)
    ap.add_argument("--min_ctg_len", type=int, default=400)
    ap.add_argument("--min_overlap", type=int, default=20)
    ap.add_argument("--max_mismatch", type=float, default=0.05)
    ap.add_argument("--nth_of_nodes", type=int, default=0)
    ap.add_argument("--r1", type=str, default="denovo/data_R1_sgrep.tsv")
    ap.add_argument("--r2", type=str, default="denovo/data_R2_sgrep.tsv")
    ap.add_argument("--n", type=int, default=None,
                    help="only assemble the first N UMIs total, across all chunks "
                         "(config: frag_de_novo.assembly_N_umi); default/empty = all UMIs")
    args = ap.parse_args()

    configure(min_ctg_len=args.min_ctg_len, min_overlap=args.min_overlap,
              max_mismatch=args.max_mismatch, out_id=args.nth_of_nodes)

    if not os.path.isdir("denovo"):
        os.makedirs("denovo")

    print("start={}".format(datetime.datetime.now()), flush=True)
    if args.n is not None:
        print("assembly_N_umi={} (denovo limited to first N UMIs)".format(args.n), flush=True)
    processed_umi = 0

    if args.end_idx is None:
        print("end_idx not specified; streaming all reads from start_idx={}".format(args.start_idx),
              flush=True)
        if args.sequence_type == "pe":
            for chunk_start, m1, m2 in _iter_pe_chunks(args.r1, args.r2, args.start_idx, args.n_line_chunk):
                if args.n is not None:
                    remaining = args.n - processed_umi
                    if remaining <= 0:
                        break
                    m2, kept = _limit_umis(m2, remaining)
                    m1 = defaultdict(list, {bc: m1[bc] for bc in m2.keys()})
                else:
                    kept = len(m2)
                print("processing chunk start_idx={} reads={}".format(
                    chunk_start, sum(len(v) // 2 for v in m2.values())), flush=True)
                _process_pe_metadata(m1, m2, args.num_processes)
                processed_umi += kept
                if args.n is not None and processed_umi >= args.n:
                    break
        else:
            for chunk_start, m2 in _iter_se_chunks(args.r2, args.start_idx, args.n_line_chunk):
                if args.n is not None:
                    remaining = args.n - processed_umi
                    if remaining <= 0:
                        break
                    m2, kept = _limit_umis(m2, remaining)
                else:
                    kept = len(m2)
                print("processing chunk start_idx={} reads={}".format(
                    chunk_start, sum(len(v) // 2 for v in m2.values())), flush=True)
                _process_se_metadata(m2, args.num_processes)
                processed_umi += kept
                if args.n is not None and processed_umi >= args.n:
                    break
    else:
        if args.end_idx <= args.start_idx:
            raise SystemExit("Invalid denovo range: start_idx={} end_idx={}".format(
                args.start_idx, args.end_idx))
        bins_ = _create_bins(args.start_idx, args.end_idx, args.n_line_chunk)
        if args.sequence_type == "pe":
            for s, e in bins_:
                if args.n is not None and processed_umi >= args.n:
                    break
                m1, m2 = defaultdict(list), defaultdict(list)
                with open(args.r1) as f1, open(args.r2) as f2:
                    for line1, line2 in itertools.islice(zip(f1, f2), s, e):
                        _add_sgrep_line(m1, line1)
                        _add_sgrep_line(m2, line2)
                if args.n is not None:
                    m2, kept = _limit_umis(m2, args.n - processed_umi)
                    m1 = defaultdict(list, {bc: m1[bc] for bc in m2.keys()})
                else:
                    kept = len(m2)
                _process_pe_metadata(m1, m2, args.num_processes)
                processed_umi += kept
        else:
            for s, e in bins_:
                if args.n is not None and processed_umi >= args.n:
                    break
                m2 = defaultdict(list)
                with open(args.r2) as f2:
                    for line in itertools.islice(f2, s, e):
                        _add_sgrep_line(m2, line)
                if args.n is not None:
                    m2, kept = _limit_umis(m2, args.n - processed_umi)
                else:
                    kept = len(m2)
                _process_se_metadata(m2, args.num_processes)
                processed_umi += kept

    if args.n is not None:
        print("total_umi_assembled={}".format(processed_umi), flush=True)
    print("end={}".format(datetime.datetime.now()), flush=True)
    subprocess.call("touch denovo/frag_denovo_done", shell=True)


# ── CLI self-test ─────────────────────────────────────────────────────────────

def _run_selftest():
    import sys
    import random

    random.seed(42)
    BASES = "ACGT"

    def _rand_seq(n):
        return "".join(random.choice(BASES) for _ in range(n))

    def _make_reads(frag, read_len, step):
        reads = []
        for s in range(0, len(frag) - read_len + 1, step):
            reads.append(frag[s:s + read_len])
        return reads

    FRAG = _rand_seq(600)

    reads_hi = _make_reads(FRAG, 150, 30)
    ctg_hi = assemble_umi(reads_hi, min_ctg=400)
    print("[hi-depth] reads={}  contig_len={}".format(len(reads_hi), len(ctg_hi[0]) if ctg_hi else 0))

    reads_lo = _make_reads(FRAG, 200, 150)
    ctg_lo = assemble_umi(reads_lo, min_ctg=100)
    print("[lo-depth] reads={}  contig_len={}".format(len(reads_lo), len(ctg_lo[0]) if ctg_lo else 0))

    ctg_1 = assemble_umi([FRAG[:250]], min_ctg=100)
    print("[1 read]   reads=1  contig_len={}".format(len(ctg_1[0]) if ctg_1 else 0))

    frag2 = _rand_seq(400)
    r1 = frag2[:200]
    r2 = rc(frag2[150:])
    ctg_rc = assemble_umi([r1, r2], min_ov=30, min_ctg=200)
    print("[RC ext]   reads=2  contig_len={}".format(len(ctg_rc[0]) if ctg_rc else 0))

    ok = True
    if not ctg_hi or (ctg_hi[0] not in FRAG and FRAG not in ctg_hi[0]):
        print("FAIL: hi-depth contig does not match fragment", file=sys.stderr)
        ok = False
    if not ctg_lo:
        print("FAIL: lo-depth assembly produced nothing", file=sys.stderr)
        ok = False
    if not ctg_rc or len(ctg_rc[0]) < 350:
        print("FAIL: RC extension too short ({})".format(len(ctg_rc[0]) if ctg_rc else 0), file=sys.stderr)
        ok = False
    if ok:
        print("correctness: OK")


if __name__ == "__main__":
    import sys
    if len(sys.argv) == 1:
        _run_selftest()
    else:
        _main_cli()