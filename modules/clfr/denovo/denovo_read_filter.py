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

Optional --ml-model/--ml-features: the greedy removal above breaks ties in
conflict degree arbitrarily (lowest read index). Passing a trained
mlpf/model_identity.lgb and its matching denovo_read_features.py output lets
the model's predicted per-read identity break those ties instead -- among
reads tied for worst, drop the one predicted lower-identity. Validated
denovo.md sec 51/59/62 (four samples, +0.03..+0.07 identity against a
greengenes local truth, all significant) before this file had any code path
for it; find_contaminants()/main() now implement that design.

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
from collections import defaultdict

# analyze_olc_shortfall.py is a tracked sibling module (same directory,
# already on sys.path when this file runs) -- no path manipulation needed.
# A second, untracked copy used to live under test/ with sys.path rigged to
# prefer it over this one; the two had drifted identical again by the time
# this was noticed, but a local dev machine silently running a different copy
# than a fresh checkout is exactly the kind of divergence that is easy to not
# notice until results stop matching.
from analyze_olc_shortfall import kmer_index, place, rc

BC_START = 5
BC_LEN = 15
# reformat_fasta2's awk emits "BX:Z:<barcode> @<header>\t<seq>" (barcode at
# BC_START:BC_START+BC_LEN, then one space, then the '@' the header line
# started with) -- the read id that denovo_read_features.py's read_id column
# actually contains is the header WITHOUT that "BX:Z:<barcode> @" prefix.
# Mirrors denovo_seed_olc.py's _add_sgrep_line, which parses the same field
# the same way (info[0][22:]) to build the same id.
ID_START = BC_START + BC_LEN + 2


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
                       scoring="ungapped", ml_scores=None,
                       coverage_guard=False, min_unique_gap=30):
    """Return the set of indices to drop (indices into `reads`).

    ml_scores, if given, is a per-read predicted identity (same order/length
    as `reads`, e.g. from mlpf/model_identity.lgb via denovo_read_features.py's
    feature set) used only to break ties in the greedy worst-degree removal
    below: among reads tied for the highest conflict degree, the one the model
    predicts as lower-identity is dropped. Without it (the default), ties
    resolve to the lowest index, which carries no meaning -- it is just
    iteration order. Validated denovo.md sec 51/59/62: +0.03..+0.07 identity
    against a greengenes local truth across four samples.

    coverage_guard (denovo.md sec 87, "P0"): skip a tie-break candidate if
    dropping it would leave a stretch of its own length unconfirmed by any
    other surviving read -- i.e. it is the pool's only evidence for that
    stretch, so dropping it is a coverage loss the ML score alone can't see.
    See greedy_conflict_drop for the mechanics.

    NOT VALIDATED -- NO-GO on hs8 3000-UMI (denovo.md sec 87): fixed 4-5 of
    the 13 hand-diagnosed severe-loss cases it was built for, but the same
    heuristic misfires often enough elsewhere in the population (16-27 newly
    introduced severe-loss cases, more than it fixes) that mean gain vs plain
    and the severe-loss rate both come out WORSE than the unguarded ML arm.
    Root cause: pool-internal pairwise agreement (<=25 reads, no reference)
    can't reliably tell "this read is the pool's only real coverage of a
    region" apart from "this read just doesn't agree with anything, which is
    what a noisy/indel read this filter is supposed to catch looks like too."
    Left here (default off) as a validated-negative result, not a shipped
    feature -- don't turn this on without new evidence.
    """
    if coverage_guard:
        n, conflicts, agree_windows, read_lens = build_conflict_graph(
            reads, same_molecule_id, min_overlap, max_reads, scoring,
            track_coverage=True)
        return greedy_conflict_drop(conflicts, n, ml_scores,
                                     agree_windows=agree_windows,
                                     read_lens=read_lens,
                                     min_unique_gap=min_unique_gap)
    n, conflicts = build_conflict_graph(
        reads, same_molecule_id, min_overlap, max_reads, scoring)
    return greedy_conflict_drop(conflicts, n, ml_scores)


def overlap_span_identity(read_a, index_b, read_b, min_overlap):
    """Like overlap_identity_ungapped, but also returns the overlap window in
    EACH read's own coordinates for the winning orientation: (ident, a_lo,
    a_hi, b_lo, b_hi). read_b's window falls out of the same place() call for
    free (lo, hi are already in read_b's own coordinates); read_a's window
    needs unflipping back out of rc() when that was the winning orientation.
    Used only by the coverage guard (track_coverage=True below) -- the
    everyday conflict-graph path doesn't need per-read windows, just the
    identity value, so it keeps using the cheaper overlap_identity_ungapped.
    """
    best = None
    for is_rc, oriented in ((False, read_a), (True, rc(read_a))):
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
        a_lo, a_hi = lo - start, hi - start
        if is_rc:
            L = len(read_a)
            a_lo, a_hi = L - a_hi, L - a_lo
        if best is None or ident > best[0]:
            best = (ident, a_lo, a_hi, lo, hi)
    return best


def build_conflict_graph(reads, same_molecule_id, min_overlap, max_reads,
                         scoring="ungapped", track_coverage=False):
    """Return ``(n_checked, adjacency)`` for one barcode's read pool, or with
    track_coverage=True, ``(n_checked, adjacency, agree_windows, read_lens)``
    where agree_windows[i] is a list of (j, lo, hi) overlap windows -- in read
    i's own coordinates -- from every j that AGREES with i (identity >=
    same_molecule_id), for the coverage guard in greedy_conflict_drop.
    track_coverage always uses the ungapped window math (overlap_span_identity)
    regardless of `scoring`, since kmer containment has no single coordinate
    window; combining --coverage-guard with --scoring kmer is rejected in
    main().

    Kept separate from the greedy policy so offline policy-learning benchmarks
    can evaluate alternative choices on the exact production graph instead of
    carrying a subtly divergent copy of graph construction.
    """
    n = min(len(reads), max_reads)
    if n < 3:
        if track_coverage:
            return n, defaultdict(set), defaultdict(list), [len(r) for r in reads[:n]]
        return n, defaultdict(set)
    indexes = [kmer_index(reads[i], 17) for i in range(n)]
    conflicts = defaultdict(set)
    agree_windows = defaultdict(list) if track_coverage else None
    score = (overlap_identity_kmer if scoring == "kmer"
             else overlap_identity_ungapped)
    for i, j in itertools.combinations(range(n), 2):
        if len(reads[i]) < min_overlap or len(reads[j]) < min_overlap:
            continue
        if track_coverage:
            hit = overlap_span_identity(reads[i], indexes[j], reads[j], min_overlap)
            if hit is None:
                continue
            ident, a_lo, a_hi, b_lo, b_hi = hit
            if ident < same_molecule_id:
                conflicts[i].add(j)
                conflicts[j].add(i)
            else:
                agree_windows[i].append((j, a_lo, a_hi))
                agree_windows[j].append((i, b_lo, b_hi))
        else:
            ident = score(reads[i], indexes[j], reads[j], min_overlap)
            if ident is not None and ident < same_molecule_id:
                conflicts[i].add(j)
                conflicts[j].add(i)
    if track_coverage:
        return n, conflicts, agree_windows, [len(r) for r in reads[:n]]
    return n, conflicts


def _uncovered_gap(i, drop, agree_windows, read_lens):
    """Longest stretch of read i's own length not spanned by any surviving
    (not-yet-dropped) agreement neighbor's window.

    Returns 0 -- i.e. no protection -- when i has ZERO surviving agreement
    neighbors, not len(i). An hs8 3000-UMI test of the first version (which
    returned the full length here) made things WORSE overall (mean gain vs
    plain went from +0.066 to -0.064, severe-loss rate 1.28%->1.68%): tracing
    the newly-introduced severe-loss cases showed 6/7 protected reads had
    n_agree_neighbors==0 -- a read that disagrees with EVERY other read in the
    pool is the isolated/indel-rich read this filter exists to remove (module
    docstring), not proof of unique real coverage nothing else happened to
    capture. That signal only means something when it's a PARTIAL gap next to
    real agreement elsewhere on the read; a read with no agreement anywhere
    gets no benefit of the doubt. denovo.md sec 87."""
    L = read_lens[i]
    ivs = sorted((lo, hi) for (j, lo, hi) in agree_windows.get(i, [])
                 if j not in drop)
    if not ivs:
        return 0
    cur = 0
    gap = 0
    for lo, hi in ivs:
        lo, hi = max(0, lo), min(L, hi)
        if lo > cur:
            gap = max(gap, lo - cur)
        cur = max(cur, hi)
    gap = max(gap, L - cur)
    return gap


def greedy_conflict_drop(conflicts, n, ml_scores=None, agree_windows=None,
                          read_lens=None, min_unique_gap=30):
    """Resolve a pre-built conflict graph with the production greedy policy:
    repeatedly drop the highest-conflict-degree read (ties broken by
    ml_scores, lowest predicted identity first, or by index if ml_scores is
    None) until no conflicts remain.

    agree_windows/read_lens (denovo.md sec 87, "P0" coverage guard): when
    given, a tied candidate is skipped -- falling through to the next-worst
    tied candidate -- if dropping it would leave >= min_unique_gap bp of its
    own length unconfirmed by any other surviving agreement neighbor, i.e. it
    is the pool's only evidence for that stretch. Motivated by trace_ties.py
    evidence (sec 86 subsection 4, failure mode 2): the ML tie-break has no
    positional information and can drop a read that happens to be the pool's
    only coverage of a region, even though its predicted identity looks fine
    in isolation. If EVERY tied candidate would violate the guard, drop
    whichever leaves the smallest gap instead of stalling, so the loop still
    terminates -- resolving the conflict always takes priority over the
    guard, the guard only picks WHICH tied read pays for it.
    """
    drop = set()
    guarded = agree_windows is not None and read_lens is not None
    while True:
        candidates = []
        for i in range(n):
            if i in drop:
                continue
            deg = len(conflicts[i] - drop)
            if deg == 0:
                continue
            i_score = ml_scores[i] if ml_scores is not None else 0.0
            candidates.append((deg, i_score, i))
        if not candidates:
            break
        max_deg = max(c[0] for c in candidates)
        tied = [c for c in candidates if c[0] == max_deg]
        tied.sort(key=lambda c: c[1])  # worst (lowest) predicted identity first
        chosen = None
        if guarded:
            for _, _, i in tied:
                if _uncovered_gap(i, drop, agree_windows, read_lens) < min_unique_gap:
                    chosen = i
                    break
            if chosen is None:
                chosen = min(
                    tied, key=lambda c: _uncovered_gap(c[2], drop, agree_windows, read_lens)
                )[2]
        else:
            chosen = tied[0][2]
        drop.add(chosen)
    return drop


# denovo_qc_probe.py imports this name; keep it pointing at the default scorer
overlap_identity = overlap_identity_ungapped


_CFG = {}


def _init_worker(cfg):
    _CFG.update(cfg)


def _filter_one(job):
    """(barcode, seqs, ml_scores) -> (barcode, sorted drop indices). Pure
    function of its input so it is safe to farm out to a process pool.
    ml_scores is None unless --ml-model/--ml-features were given; workers
    never load the model themselves -- scoring happens once in the main
    process (see handle_batch in main()) and only the small per-barcode score
    list travels through the job tuple."""
    barcode, seqs, ml_scores = job
    drop = find_contaminants(seqs, _CFG["same_molecule_id"],
                              _CFG["min_overlap"], _CFG["max_reads"],
                              _CFG.get("scoring", "ungapped"),
                              ml_scores=ml_scores,
                              coverage_guard=_CFG.get("coverage_guard", False),
                              min_unique_gap=_CFG.get("min_unique_gap", 30))
    return barcode, sorted(drop)


def iter_feature_groups(path):
    """Yield (barcode, [(read_id, [13 feature floats]), ...]) from a
    denovo_read_features.py output TSV. That file is barcode-sorted by
    construction (it streams a barcode-sorted FASTQ), the same invariant
    iter_barcode_groups relies on below, so the two can be merge-joined
    without loading either into memory."""
    import denovo_read_features as drf
    n_cols = len(drf.FIELDS)
    cur_bc, rows = None, []
    with open(path) as fh:
        fh.readline()  # header
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < n_cols:
                continue
            rid, bc = f[0], f[1]
            vals = [float(x) for x in f[2:n_cols]]
            if bc != cur_bc:
                if rows:
                    yield cur_bc, rows
                cur_bc, rows = bc, []
            rows.append((rid, vals))
    if rows:
        yield cur_bc, rows


def iter_joined_groups(r2_path, features_path):
    """Zip iter_barcode_groups(r2_path) with iter_feature_groups(features_path)
    on barcode via an ordinary sorted merge (both streams share that order).
    Yields (barcode, lines, seqs, feat_rows); feat_rows is None when a barcode
    has no match on the features side, so its group falls back to the
    default (unscored) tie-break rather than guessing."""
    feat_iter = iter_feature_groups(features_path)
    fbc, frows = next(feat_iter, (None, None))
    for bc, lines, seqs in iter_barcode_groups(r2_path):
        while fbc is not None and fbc < bc:
            fbc, frows = next(feat_iter, (None, None))
        yield (bc, lines, seqs, frows) if fbc == bc else (bc, lines, seqs, None)


def score_batch(booster, batch):
    """[(bc, lines, seqs, feat_rows), ...] -> ml_scores list, one entry per
    group (None where that group can't be scored: no --ml-model, or no/partial
    match on the features side). One booster.predict() call for the whole
    batch instead of one per barcode, matching the batching this file already
    does for I/O -- millions of individual predict() calls would each carry
    fixed native-call overhead on top of a tiny (~pool-size x 13) matrix,
    adding up at production scale."""
    scores = [None] * len(batch)
    if booster is None:
        return scores
    X, spans = [], []
    for gi, (_bc, _lines, _seqs, frows) in enumerate(batch):
        if not frows:
            continue
        start = len(X)
        X.extend(vals for _rid, vals in frows)
        spans.append((gi, start, len(X)))
    if not X:
        return scores
    preds = booster.predict(X)
    for gi, start, end in spans:
        _bc, lines, _seqs, frows = batch[gi]
        ids = [ln.split("\t", 1)[0][ID_START:] for ln in lines]
        rid_score = {rid: preds[start + k] for k, (rid, _vals) in enumerate(frows)}
        group_scores = [rid_score.get(rid) for rid in ids]
        if all(s is not None for s in group_scores):
            scores[gi] = group_scores
    return scores


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
                     help="barcodes held in memory per parallel batch, and (with "
                          "--ml-model) reads scored per model.predict() call")
    ap.add_argument("--ml-model", help="lightgbm model file (e.g. "
                     "mlpf/model_identity.lgb) predicting per-read identity; "
                     "when set, its prediction breaks conflict-graph ties "
                     "instead of the arbitrary default. Requires --ml-features. "
                     "denovo.md sec 51/59/62")
    ap.add_argument("--ml-features", help="denovo_read_features.py output TSV "
                     "for the same reads as --r2 (barcode-sorted, from the "
                     "same source FASTQ); required with --ml-model")
    ap.add_argument("--coverage-guard", action="store_true",
                     help="skip a tie-break candidate (falling through to the "
                          "next-worst tied one) if dropping it would leave a "
                          "stretch of its own length with no other confirming "
                          "read in the pool; see greedy_conflict_drop. "
                          "denovo.md sec 87 (P0). Requires --scoring ungapped")
    ap.add_argument("--min-unique-gap", type=int, default=30,
                     help="bp threshold for --coverage-guard: a candidate is "
                          "guarded if dropping it leaves this many bp of its "
                          "own length unconfirmed by any surviving neighbor")
    args = ap.parse_args()
    if bool(args.ml_model) != bool(args.ml_features):
        ap.error("--ml-model and --ml-features must be given together")
    if args.coverage_guard and args.scoring != "ungapped":
        ap.error("--coverage-guard requires --scoring ungapped (kmer "
                  "containment has no single coordinate window)")

    cfg = {"same_molecule_id": args.same_molecule_id,
           "min_overlap": args.min_overlap,
           "max_reads": args.max_reads,
           "scoring": args.scoring,
           "coverage_guard": args.coverage_guard,
           "min_unique_gap": args.min_unique_gap}

    booster = None
    if args.ml_model:
        import lightgbm as lgb
        booster = lgb.Booster(model_file=args.ml_model)

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
        ml_scores = score_batch(booster, batch)
        jobs = [(bc, seqs, ml_scores[gi])
                for gi, (bc, _lines, seqs, _frows) in enumerate(batch)]
        if pool is not None:
            # chunksize=1: per-barcode cost is quadratic in read count, so a
            # few deep barcodes dominate -- the same imbalance denovo_seed_olc
            # hit with default chunking.
            results = pool.map(_filter_one, jobs, chunksize=1)
        else:
            results = [_filter_one(j) for j in jobs]
        for (bc, lines, _seqs, _frows), (_bc, drop) in zip(batch, results):
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
        groups = (iter_joined_groups(args.r2, args.ml_features) if booster is not None
                  else ((bc, lines, seqs, None)
                        for bc, lines, seqs in iter_barcode_groups(args.r2)))
        batch = []
        for group in groups:
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
