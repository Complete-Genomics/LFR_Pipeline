#!/usr/bin/env python3
"""Pre-assembly read QUALITY filter (reference-free).

Named for what it was measured to do, not what it was designed to do. The
original intent was to separate reads from different source molecules sharing
a barcode; tracing it against the ZymoBIOMICS reference showed it actually
removes indel-rich, error-rich reads, and that this is where its benefit comes
from (denovo.md sec 38).

The evidence, on reads aligned to the TRUE Zymo reference:
  dropped reads   mean identity 92.84, 68% carry indels, 1.16 gaps/read
  kept reads      mean identity 95.81, 52% carry indels, 0.67 gaps/read
and separately, dropped reads align to their own barcode's contig at >=97%
identity in 100% of cases once gaps are allowed -- i.e. they are the same
molecule, just noisier. So this is a quality filter, not a contamination
filter, and the reads it removes are real 16S from the right organism.

Removing them still helps, because indels in reads propagate into the greedy
consensus: on the Zymo control it moves assembly identity 94.13 -> 95.02 and
raises the >=97% fraction from 30.2% to 40.8%.

Mechanics: build a graph over a barcode's reads where an overlapping pair
scoring below `same_molecule_id` is a CONFLICT edge, then greedily remove the
highest-conflict-degree read until no conflicts remain. Reads that do not
overlap at all get no edge -- treating "unconnected" as "different" would
shred clean barcodes, which is how an earlier attempt failed (sec 30).

Two scorers, selected by --scoring:
  ungapped (default) -- column-by-column at a single k-mer offset. An indel
      shifts every downstream base and tanks the score, so indel-bearing reads
      become conflicts. That is exactly the behaviour that makes this a
      quality filter, and it measured best.
  kmer -- indel-tolerant k-mer containment (see overlap_identity_kmer). It is
      the honest way to ask "are these different molecules", and it drops far
      fewer reads (5.7% vs 13.6% on Zymo), but assembly identity comes out
      worse (94.54 vs 95.02) precisely because it spares the noisy reads.
      Use it when the question really is molecule separation, not quality.
"""
import argparse
import csv
import itertools
import os
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + "/test")
from analyze_olc_shortfall import kmer_index, place, rc  # noqa: E402

BC_START = 5
BC_LEN = 15


def _kmers(seq, k):
    return {seq[i:i + k] for i in range(len(seq) - k + 1)}


def overlap_identity_ungapped(read_a, index_b, read_b, min_overlap):
    """Column-by-column identity at a single k-mer offset (no indel tolerance).

    Default scorer. An indel shifts everything downstream, so a read carrying
    one scores far below threshold and is dropped -- which is why this doubles
    as an indel/error-rich read filter and why it outperforms the
    indel-tolerant scorer on assembly accuracy (denovo.md sec 38).
    """
    best = None
    for oriented in (read_a, rc(read_a)):
        hit = place(oriented, read_b, index_b)
        if hit is None:
            continue
        start = hit[0]
        lo, hi = max(0, start), min(len(read_b), start + len(oriented))
        if hi - lo < min_overlap:
            continue
        matches = sum(1 for pos in range(lo, hi)
                      if oriented[pos - start] == read_b[pos])
        ident = matches / (hi - lo)
        if best is None or ident > best:
            best = ident
    return best


def overlap_identity_kmer(read_a, index_b, read_b, min_overlap, k=17):
    """Indel-tolerant estimate of overlap identity, or None when the two reads
    do not overlap by at least min_overlap.

    `place` is still used to FIND the overlapping region -- k-mer anchors
    locate it correctly even with indels. What changed is how the region is
    SCORED. Column-by-column comparison at a single offset was measured to be
    the wrong tool here (denovo.md sec 38): one indel shifts every downstream
    base, identity collapses far below any sane threshold, and a read from the
    same molecule gets called a conflict. On real data that made this filter an
    indel-read remover rather than the contamination remover it claimed to be.

    Instead the overlap is scored by k-mer containment, which an indel only
    damages locally (it destroys the ~k k-mers spanning the indel, not
    everything after it). Containment is converted back onto the identity scale
    with the standard estimator identity ~= containment**(1/k), since a
    per-base identity p leaves roughly p**k of the k-mers intact -- so the
    caller's threshold keeps meaning roughly what it did before.

    Calibrated on Zymo reads labelled by which reference species they align to:
    same-species pairs score median 0.984, different-species pairs 0.910, so a
    0.90 cut mislabels 6.7% of same-species pairs while catching 40% of
    different-species ones. Separation is limited because 16S conserved regions
    keep containment high even between species.
    """
    best = None
    for oriented in (read_a, rc(read_a)):
        hit = place(oriented, read_b, index_b)
        if hit is None:
            continue
        start = hit[0]
        lo, hi = max(0, start), min(len(read_b), start + len(oriented))
        if hi - lo < min_overlap:
            continue
        sub_b = read_b[lo:hi]
        sub_a = oriented[max(0, lo - start):max(0, lo - start) + (hi - lo)]
        if len(sub_a) < k or len(sub_b) < k:
            continue
        ka, kb = _kmers(sub_a, k), _kmers(sub_b, k)
        denom = min(len(ka), len(kb))
        if not denom:
            continue
        containment = len(ka & kb) / denom
        ident = containment ** (1.0 / k) if containment > 0 else 0.0
        if best is None or ident > best:
            best = ident
    return best


def find_contaminants(reads, same_molecule_id, min_overlap, max_reads,
                       scoring="ungapped"):
    """Return the set of indices to drop (indices into `reads`)."""
    score = (overlap_identity_kmer if scoring == "kmer"
             else overlap_identity_ungapped)
    n = min(len(reads), max_reads)
    if n < 3:
        return set()
    indexes = [kmer_index(reads[i], 17) for i in range(n)]
    conflicts = defaultdict(set)
    for i, j in itertools.combinations(range(n), 2):
        if len(reads[i]) < min_overlap or len(reads[j]) < min_overlap:
            continue
        ident = score(reads[i], indexes[j], reads[j], min_overlap)
        if ident is not None and ident < same_molecule_id:
            conflicts[i].add(j)
            conflicts[j].add(i)

    drop = set()
    while True:
        worst, worst_deg = None, 0
        for i in range(n):
            if i in drop:
                continue
            deg = len(conflicts[i] - drop)
            if deg > worst_deg:
                worst, worst_deg = i, deg
        if worst is None:
            break
        drop.add(worst)
    return drop


# denovo_qc_probe.py imports this name; keep it pointing at the default scorer
overlap_identity = overlap_identity_ungapped


_CFG = {}


def _init_worker(cfg):
    _CFG.update(cfg)


def _filter_one(job):
    """(barcode, seqs) -> (barcode, sorted drop indices). Pure function of its
    input so it is safe to farm out to a process pool."""
    barcode, seqs = job
    drop = find_contaminants(seqs, _CFG["same_molecule_id"],
                              _CFG["min_overlap"], _CFG["max_reads"],
                              _CFG.get("scoring", "ungapped"))
    return barcode, sorted(drop)


def iter_barcode_groups(path):
    """Yield (barcode, lines, seqs) groups. Relies on the TSV already being
    barcode-sorted, which denovo_preprocess.smk guarantees."""
    cur_bc, lines, seqs = None, [], []
    with open(path) as fh:
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2:
                continue
            bc = fields[0][BC_START:BC_START + BC_LEN]
            if bc != cur_bc:
                if lines:
                    yield cur_bc, lines, seqs
                cur_bc, lines, seqs = bc, [], []
            lines.append(line)
            seqs.append(fields[1])
    if lines:
        yield cur_bc, lines, seqs


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--r2", required=True, help="sorted sgrep TSV")
    ap.add_argument("--out", required=True, help="filtered sorted sgrep TSV")
    ap.add_argument("--dropped-out", help="sgrep TSV of the reads this pass "
                     "dropped, same format as --out, for salvage/inspection "
                     "(e.g. reassembling just a barcode's dropped reads to see "
                     "if they form a second real molecule rather than noise)")
    ap.add_argument("--report", help="per-barcode drop counts")
    ap.add_argument("--same-molecule-id", type=float, default=0.90,
                     help="overlapping reads below this identity are treated as "
                          "coming from different source molecules. Validated on "
                          "the ZymoBIOMICS control (denovo.md sec 32): 0.90 gives "
                          "+0.89 mean identity, while 0.97 over-fires (30%% of "
                          "reads dropped) and makes the assembly worse")
    ap.add_argument("--scoring", choices=["ungapped", "kmer"], default="ungapped",
                     help="how an overlapping read pair is scored; see module "
                          "docstring. ungapped also removes indel-rich reads and "
                          "measured best for assembly accuracy")
    ap.add_argument("--min-overlap", type=int, default=150)
    ap.add_argument("--max-reads", type=int, default=25,
                     help="cap on reads examined per barcode (pairwise cost is "
                          "quadratic); reads beyond the cap are always kept")
    ap.add_argument("--num_processes", type=int, default=1)
    ap.add_argument("--batch-barcodes", type=int, default=20000,
                     help="barcodes held in memory per parallel batch")
    args = ap.parse_args()

    cfg = {"same_molecule_id": args.same_molecule_id,
           "min_overlap": args.min_overlap,
           "max_reads": args.max_reads,
           "scoring": args.scoring}

    n_bc = n_reads = n_dropped = 0
    report = None
    rep = None
    if args.report:
        report = open(args.report, "w", newline="")
        rep = csv.writer(report, delimiter="\t")
        rep.writerow(["barcode", "n_reads", "n_checked", "n_dropped"])

    pool = None
    if args.num_processes > 1:
        import multiprocessing as mp
        pool = mp.Pool(args.num_processes, initializer=_init_worker, initargs=(cfg,))
    else:
        _init_worker(cfg)

    def handle_batch(batch, out, dropped_out):
        nonlocal n_bc, n_reads, n_dropped
        if not batch:
            return
        jobs = [(bc, seqs) for bc, _lines, seqs in batch]
        if pool is not None:
            # chunksize=1: per-barcode cost is quadratic in read count, so a
            # few deep barcodes dominate -- the same imbalance denovo_seed_olc
            # hit with default chunking.
            results = pool.map(_filter_one, jobs, chunksize=1)
        else:
            results = [_filter_one(j) for j in jobs]
        for (bc, lines, _seqs), (_bc, drop) in zip(batch, results):
            dropset = set(drop)
            n_bc += 1
            n_reads += len(lines)
            n_dropped += len(dropset)
            for idx, line in enumerate(lines):
                if idx in dropset:
                    if dropped_out:
                        dropped_out.write(line)
                else:
                    out.write(line)
            if rep:
                rep.writerow([bc, len(lines), min(len(lines), args.max_reads),
                              len(dropset)])

    dropped_fh = open(args.dropped_out, "w") if args.dropped_out else None
    with open(args.out, "w") as out:
        batch = []
        for group in iter_barcode_groups(args.r2):
            batch.append(group)
            if len(batch) >= args.batch_barcodes:
                handle_batch(batch, out, dropped_fh)
                batch = []
        handle_batch(batch, out, dropped_fh)
    if dropped_fh:
        dropped_fh.close()

    if pool is not None:
        pool.close()
        pool.join()
    if report:
        report.close()
    print(f"barcodes={n_bc}")
    print(f"reads_in={n_reads}")
    print(f"reads_dropped={n_dropped} ({100*n_dropped/n_reads:.2f}%)" if n_reads else "")


if __name__ == "__main__":
    main()