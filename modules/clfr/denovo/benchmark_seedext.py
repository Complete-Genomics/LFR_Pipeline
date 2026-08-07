"""
Benchmark + smoke-test denovo_seed_olc.assemble_umi on real data.

Usage (run from the Snakemake work dir that contains denovo/data_R2_sorted.tsv):

    # timing only, single core
    python3 /path/to/benchmark_seedext.py --n 1000 --r2 denovo/data_R2_sorted.tsv

    # timing with N worker processes -- single-core numbers alone don't
    # tell you real production throughput; the deployed path always runs
    # with multiprocessing (denovo_seed_olc.py --num_processes), and the
    # achievable speedup from more cores has to be measured on the actual
    # target box (hyperthreading / NUMA / other load all affect it --
    # e.g. on one dev machine 12 processes was slower than 8, from
    # oversubscription). Try a few values and look for where it stops
    # helping, don't just take the highest core count.
    python3 /path/to/benchmark_seedext.py --n 1000 --r2 denovo/data_R2_sorted.tsv \\
        --num_processes 30

    # smoke test: write FASTA and inspect
    python3 /path/to/benchmark_seedext.py --n 100 --r2 denovo/data_R2_sorted.tsv \\
        --output smoke_contigs.fa
    grep -c '>' smoke_contigs.fa          # count contigs
    awk '/^>/{next}{print length}' smoke_contigs.fa | sort -n  # contig lengths

Prints: per-UMI latency, throughput, contig yield, 1M/3M extrapolation.
"""

import argparse
import multiprocessing as mp
import time
import sys
import os
from collections import defaultdict
from functools import partial
from pathlib import Path

# ── import assembler ──────────────────────────────────────────────────────────
# module level (not inside main()): must be importable/picklable by name for
# spawned worker processes (macOS/Windows default start method) to find it.
script_dir = Path(__file__).parent
sys.path.insert(0, str(script_dir))
from denovo_seed_olc import assemble_umi

BC_START = 5
BC_LEN = 15


def _assemble_one(seqs, min_ov, min_ctg):
    """Top-level (picklable) per-barcode worker for the --num_processes>1 path."""
    return assemble_umi(seqs, min_ov=min_ov, min_ctg=min_ctg)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--n",  type=int, default=1000, help="number of UMIs to benchmark [1000]")
    parser.add_argument("--r2", type=str, default="denovo/data_R2_sorted.tsv",
                        help="path to sorted TSV (R2) [denovo/data_R2_sorted.tsv]")
    parser.add_argument("--min_ctg", type=int, default=400)
    parser.add_argument("--min_ov",  type=int, default=20)
    parser.add_argument("--num_processes", type=int, default=1,
                        help="parallel worker processes [1]. Only affects the "
                             "assemble_umi loop below -- data loading is always "
                             "single-threaded.")
    parser.add_argument("--output",  type=str, default=None,
                        help="write assembled contigs to this FASTA file (smoke test)")
    args = parser.parse_args()

    # ── load first N barcodes from sgrep TSV ─────────────────────────────────
    meta = defaultdict(list)
    r2_path = args.r2

    if not os.path.exists(r2_path):
        sys.exit(f"ERROR: file not found: {r2_path}\n"
                 f"Run from the work dir containing {r2_path}, or pass --r2 <path>")

    print(f"Reading {r2_path} ...", flush=True)
    with open(r2_path) as f:
        for line in f:
            info = line.rstrip("\n").split("\t")
            if len(info) < 2:
                continue
            bc  = info[0][BC_START : BC_START + BC_LEN]
            seq = info[1]
            meta[bc].append(seq)
            if len(meta) >= args.n and bc in meta:
                # stop once we have collected enough distinct barcodes
                # (keep reading until the current bc is fully loaded)
                pass
            if len(meta) > args.n:
                break

    barcodes = list(meta.keys())[: args.n]
    print(f"Loaded {len(barcodes)} UMIs  "
          f"(reads/UMI: min={min(len(meta[b]) for b in barcodes)}, "
          f"max={max(len(meta[b]) for b in barcodes)}, "
          f"mean={sum(len(meta[b]) for b in barcodes)/len(barcodes):.1f})")

    # ── benchmark ─────────────────────────────────────────────────────────────
    t0 = time.perf_counter()
    if args.num_processes == 1:
        results = [_assemble_one(meta[bc], args.min_ov, args.min_ctg) for bc in barcodes]
    else:
        worker = partial(_assemble_one, min_ov=args.min_ov, min_ctg=args.min_ctg)
        with mp.Pool(args.num_processes) as pool:
            results = pool.map(worker, [meta[bc] for bc in barcodes])
    elapsed = time.perf_counter() - t0

    # ── optional FASTA output + yield stats (single-threaded, after timing) ──
    n_contigs = 0
    n_empty = 0
    contig_lens = []
    out_fh = open(args.output, "w") if args.output else None
    for bc, ctgs in zip(barcodes, results):
        if ctgs:
            n_contigs += len(ctgs)
            contig_lens.extend(len(c) for c in ctgs)
            if out_fh:
                for i, seq in enumerate(ctgs):
                    out_fh.write(">{}>k41_{} len={} reads={}\n{}\n".format(
                        bc, i, len(seq), len(meta[bc]), seq))
        else:
            n_empty += 1
    if out_fh:
        out_fh.close()
        print("FASTA written to: {}  ({} contigs)".format(args.output, n_contigs))

    # ── report ────────────────────────────────────────────────────────────────
    n = len(barcodes)
    print(f"\n{'='*50}")
    print(f"UMIs benchmarked : {n}")
    print(f"Worker processes : {args.num_processes}")
    print(f"Total time       : {elapsed:.3f} s")
    print(f"Per-UMI latency  : {elapsed/n*1000:.3f} ms")
    label = "single core" if args.num_processes == 1 else f"{args.num_processes} processes"
    print(f"Throughput       : {n/elapsed:.0f} UMI/s  ({label})")
    print(f"\nContig yield:")
    print(f"  UMIs with contig : {n - n_empty} / {n}  ({100*(n-n_empty)/n:.1f}%)")
    if contig_lens:
        print(f"  Contig len mean  : {sum(contig_lens)/len(contig_lens):.0f} bp")
        print(f"  Contig len range : {min(contig_lens)} – {max(contig_lens)} bp")
    print(f"{'='*50}")

    # ── extrapolate ───────────────────────────────────────────────────────────
    for total in (1_000_000, 3_000_000):
        t_est = elapsed / n * total
        print(f"Estimated time for {total//1_000_000}M UMIs, "
              f"{label} : {t_est/60:.1f} min")


if __name__ == "__main__":
    main()
