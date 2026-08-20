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
    python3 denovo_seed_olc.py

Standalone CLI (drop-in replacement for denovo_clfr_ram.py --module denovo_parallel)
-------------------------------------------------------------------------------------
Run directly from a Snakemake work dir containing denovo/data_R1_sorted.tsv
and denovo/data_R2_sorted.tsv (same layout denovo_clfr_ram.py expects).
No megahit binary, no tmp_dir, no subprocess fork per UMI.

    python3 denovo_seed_olc.py \\
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
        --r2 denovo/data_R2_sorted.tsv \\
        --min_ctg 400

Output printed:
    per-UMI latency, throughput (UMI/s), contig yield, 1M/3M extrapolation

Polish step (post-assembly majority-vote consensus correction)
------------------------------------------------------------------
Borrows the "majority vote per position" idea from
cgi/pipeline/tests/dev/LFR_pipeline/src/sparse_denovo_trasm's CallArray
pileup, without its k-mer-index / numpy / homopolymer machinery.

After assemble_umi() (or the mappy path) produces a draft contig,
polish_contig() re-aligns every UMI read against it via a cheap k-mer
offset vote (no full DP alignment), tallies per-position base votes,
and flips a position only when an alternative base out-votes the
original (which starts with 1 implicit vote) by >= vote_concordance
with >= min_coverage total votes. Corrects substitution errors from
the single-read greedy merge; does NOT handle indels (e.g. homopolymer
length errors) -- a read with an indel relative to the contig simply
fails the concordance check and is excluded from voting rather than
corrupting the consensus. On by default; disable with
configure(polish=False) or CLI --no_polish.

Internal-anchor extension (fallback for boundary-overlap failures)
------------------------------------------------------------------
suffix_prefix_overlap only ever compares a candidate read's literal
first/last N bases against the contig's boundary. Real data (verified
on a real 16S rRNA barcode, BLAST-confirmed truth sequence) can have
reads whose genuinely-matching content does not start at the read's
own edge -- e.g. a short PCR-chimera artifact or quality-degraded
stretch fused onto an otherwise-accurate long middle section. Such
reads pass right by suffix_prefix_overlap since it never looks inside
a read for a usable anchor.

internal_anchor_extend_indexed() is a fallback tried only after
ordinary boundary extension is fully exhausted for a contig: using a
k-mer index built once per UMI (_build_pool_kmer_index), it looks up
matches to the contig's boundary region anywhere in the remaining
candidates (not just their start), verifies a real, honest overlap from
there (>= internal_min_verify, deliberately longer than the default
boundary min_ov to guard against short universally-conserved-region
false positives -- e.g. bacterial 16S primer sites shared across
unrelated organisms), and on success extends the contig with only the
NEW sequence past the verified overlap -- discarding whatever noise
came before the anchor in the candidate.

Verified on real data: recovered a barcode's longest contig from 575bp
to 1115bp (97.2% identity, correct monotonic alignment, to the
BLAST-confirmed 1368bp megahit truth), with no measured increase in
chimeric-merge rate on synthetic shared-conserved-motif stress tests.

The forward pair above only searches for the anchor inside the
candidate, assuming the contig's own boundary is reliable -- true except
for the very first seed's own raw, unverified edge (rare, but the one
case where this assumption breaks: sorting by post-trim length picks the
longest read as seed, which is at best a weak proxy for "this read's own
edges are clean" -- trimming mostly reflects adapter/insert-size
geometry, not base-level error rate at the very ends). A reverse pair
(_internal_anchor_extend_3prime_reverse_indexed /
_internal_anchor_extend_5prime_reverse_indexed) trusts a candidate's own
content instead and searches for where it anchors inside the contig,
truncating and replacing the contig's own noisy edge rather than
carrying it forward forever. Tried only as a last resort, after both
forward directions fail, so existing behavior is unchanged whenever a
forward match exists. On by default; disable with
configure(use_internal_anchor=False) or CLI --no_internal_anchor.
"""

import gzip
import itertools
import os
from collections import defaultdict, Counter
from functools import lru_cache

# ── module-level config (set via configure() before multiprocessing) ──────────

_CFG = {
    "min_ctg":   400,
    "min_ov":    20,
    "max_mm":    0.05,
    "out_id":    0,
    "out_file":  "denovo/final_contigs_{id}.fa",
    "seed_k":    10,    # k-mer size for overlap pre-filter
    "adaptive_seed_k": False,  # retry a strict seed_k with fallback_seed_k only if no valid contig
    "fallback_seed_k": 10,
    "use_mappy": False,  # OFF by default: _assemble_umi_mappy has none of assemble_umi's
                         # correctness hardening (no mismatch-rate check, no chimeric-merge
                         # safety, always collapses a barcode to a single contig with no
                         # multi-fragment/collision handling, no low-depth single-read
                         # support) and every fix validated across this project's history
                         # only applies to assemble_umi. A real 1.5M production run was
                         # confirmed to silently take this path whenever the environment
                         # happened to have `mappy` importable (auto-detect via None),
                         # making every one of those fixes a no-op in production without
                         # anyone intending it. None = auto-detect (legacy; do not use for
                         # production), True/False = force.
    "polish":              True,  # majority-vote consensus correction after assembly
    "polish_min_coverage": 3,     # min total votes (incl. 1 implicit vote for original base) to consider flipping
    "polish_vote_concordance": 0.6,  # winning base must hold >= this fraction of votes to flip
    "polish_kmer_step":    5,     # stride for sampling k-mers when re-aligning reads to the contig
    "use_internal_anchor":  True,  # fallback: k-mer anchor anywhere in a read, not just its literal ends
    "internal_min_verify":  60,    # min confirmed overlap length for the internal-anchor fallback
    "max_contigs":          8,     # maximum final contigs emitted per UMI
    "use_minimizer_dedup":  False,  # collapse near-identical reads before assembly -- OFF: confirmed conserved-motif collision bug, see _canonical_minimizer
    "use_cross_attempt_evidence": True,  # collective rescue may vote using reads already consumed by earlier raw-building attempts, not just the current attempt's leftovers -- see assemble_umi docstring for the isolation A/B this default is pending on
    "min_join_support": 1,  # verified raw reads required across a newly created join; 1 preserves legacy assembly
    "join_trace": None,  # path to write a per-join TSV (see _write_join_trace) -- None (default) disables ALL join-level bookkeeping/bridge-support computation, so a barcode run with this off pays zero cost and produces byte-identical assembly output to before this instrumentation existed. Explicit opt-IN only (CLI --join-trace PATH), matching this file's existing convention (see _configure_from_args's docstring for why opt-OUT flags have twice caused a silent-default production incident here) -- pure observation, never feeds back into any join/extension decision.
}


def configure(min_ctg_len=400, min_overlap=20, max_mismatch=0.05,
              out_id=0, out_file="denovo/final_contigs_{id}.fa", use_mappy=False,
              polish=True, polish_min_coverage=3, polish_vote_concordance=0.6,
              polish_kmer_step=5, use_internal_anchor=True, internal_min_verify=60,
              max_contigs=8, use_minimizer_dedup=False,
              use_cross_attempt_evidence=True, adaptive_seed_k=False,
              fallback_seed_k=10, seed_k=10, min_join_support=1,
              join_trace=None):
    """Call once in the parent process before spawning Pool workers."""
    if min_join_support < 1:
        raise ValueError("min_join_support must be >= 1")
    _CFG["min_ctg"]   = min_ctg_len
    _CFG["min_ov"]    = min_overlap
    _CFG["max_mm"]    = max_mismatch
    _CFG["seed_k"]    = seed_k
    _CFG["out_id"]    = out_id
    _CFG["out_file"]  = out_file
    _CFG["use_mappy"] = use_mappy
    _CFG["polish"]                   = polish
    _CFG["polish_min_coverage"]      = polish_min_coverage
    _CFG["polish_vote_concordance"]  = polish_vote_concordance
    _CFG["polish_kmer_step"]         = polish_kmer_step
    _CFG["use_internal_anchor"]      = use_internal_anchor
    _CFG["internal_min_verify"]      = internal_min_verify
    _CFG["max_contigs"]              = max_contigs
    _CFG["use_minimizer_dedup"]      = use_minimizer_dedup
    _CFG["use_cross_attempt_evidence"] = use_cross_attempt_evidence
    _CFG["adaptive_seed_k"] = adaptive_seed_k
    _CFG["fallback_seed_k"] = fallback_seed_k
    _CFG["min_join_support"] = min_join_support
    _CFG["join_trace"] = join_trace


# ── sequence utilities ────────────────────────────────────────────────────────

_RC = str.maketrans("ACGT", "TGCA")

_DEDUP_MINIMIZER_K = 21
_DEDUP_KEEP_PER_CLUSTER = 2


def rc(seq):
    return seq.translate(_RC)[::-1]


def _canonical_minimizer(seq, k=_DEDUP_MINIMIZER_K):
    """
    Cheapest possible near-duplicate signature: the lexicographically
    smallest k-mer across both orientations of the read. O(len(seq)), no
    pairwise comparison -- same principle as BBMap Clumpify's reference-free
    duplicate detection (canonical k-mer binning), not the O(n^2)
    difflib.SequenceMatcher pairwise comparison this replaces (measured at
    5.1s / 30.1s of clustering cost alone for 27 / 64 real reads -- a net
    loss bigger than any downstream saving; the minimizer version costs
    1.9ms / 4.6ms for the same UMIs).

    CONFIRMED FAILURE MODE, why _minimizer_dedupe defaults to off: a SINGLE
    global minimizer per read is a much weaker signal than real minimizer-
    based dedup tools use (they window it / require several agreeing
    minimizers). On 16S data specifically, an AT-rich conserved motif can
    easily BE the lexicographically smallest k-mer in many otherwise
    completely different reads that tile unrelated true positions and share
    nothing else -- confirmed directly: 5 synthetic reads with distinct
    500bp unique payloads but one common 21-A conserved motif all produced
    the identical minimizer and collapsed to 2 "kept" reads, silently
    discarding 3 with real, non-redundant information. A full-scale rerun
    with this dedup enabled by default showed exactly this signature at
    scale (vs the prior baseline): total contigs 124,240 -> 109,039, >=1kb
    contigs 21,139 -> 17,595, net length change ~-1.06 Mb. Do not re-enable
    by default without fixing the single-minimizer collision risk (e.g.
    windowed/multiple minimizers, or verifying actual shared sequence length
    before collapsing) and re-validating at full scale, not just the
    hand-picked high-PCR-redundancy cases this was originally tuned on.
    """
    best = None
    for variant in (seq, rc(seq)):
        for i in range(len(variant) - k + 1):
            kmer = variant[i:i + k]
            if best is None or kmer < best:
                best = kmer
    return best


def _minimizer_dedupe(seqs, k=_DEDUP_MINIMIZER_K, keep_n=_DEDUP_KEEP_PER_CLUSTER):
    """
    Collapse near-identical reads sharing a canonical minimizer down to at
    most `keep_n` representatives (longest first) per bucket, before the
    greedy assembler ever sees them. Reads too short for a full k-mer have
    no minimizer to bucket by and are always kept as-is.

    This is a distinct problem from the exact-duplicate dedup already done
    by assemble_umi's own `sorted(set(seqs), ...)`: real production UMIs
    routinely carry many NEAR-identical reads (same narrow region, distinct
    only by independent sequencing errors -- e.g. one real 64-read UMI had
    a cluster of 16 such reads) that exact-string dedup can't touch. Highly
    redundant reads don't just cost time -- they can actively hurt the
    greedy assembler: each one that doesn't merge into an existing attempt
    can spawn ANOTHER redundant raw contig over the same already-covered
    region instead of the pool's attention going toward genuinely
    unexplored territory. Verified on real UMIs: the 64-read case above
    went from 555ms/1121bp-best to 149ms/1516bp-best after this dedup --
    faster AND a better assembly, not a trade-off.

    Does NOT by itself bound worst-case cost: UMIs with genuine (not
    redundant) high molecular diversity keep most of their reads through
    this filter -- see _COLLECTIVE_MAX_CROSS_ATTEMPT_READS and
    _MAX_RAW_CONTIGS_FOR_WIDE_MERGE for the fragmentation-driven pathological
    cases this alone does not fix.
    """
    seqs_u = list(dict.fromkeys(seqs))
    buckets = defaultdict(list)
    kept = []
    for s in seqs_u:
        if len(s) < k:
            kept.append(s)
            continue
        buckets[_canonical_minimizer(s, k)].append(s)
    for members in buckets.values():
        members.sort(key=len, reverse=True)
        kept.extend(members[:keep_n])
    return kept


def _kmer_positions(seq, k, start=0, end=None):
    """Map kmer -> list of its ABSOLUTE positions within seq[start:end]."""
    s = seq[start:end]
    if len(s) < k:
        return {}
    positions = {}
    for i in range(len(s) - k + 1):
        positions.setdefault(s[i:i + k], []).append(start + i)
    return positions


# ── core overlap ──────────────────────────────────────────────────────────────

def _within_mismatch_budget(a_tail, b_head, ov, max_mm):
    """
    True iff mismatches between a_tail and b_head satisfy mm/ov <= max_mm.

    Exits as soon as the running mismatch count makes that impossible --
    mismatches only accumulate as the loop progresses, so once mm/ov
    exceeds max_mm it can never recover to pass. Profiling showed that,
    for real 16S data, the overwhelming majority of candidate ov values
    tested here are coincidental k-mer hits with no real underlying
    overlap, which mismatch heavily almost immediately. Rejecting those
    early turns what used to be a full O(ov) scan (199M+ character
    comparisons across 227K calls, ~10s of tottime alone) into a
    handful of comparisons per rejected candidate.
    """
    mm = 0
    for x, y in zip(a_tail, b_head):
        if x != y:
            mm += 1
            if mm > ov * max_mm:
                return False
    return True


def suffix_prefix_overlap(a, b, min_ov, max_mm, seed_k=10):
    """
    Return the length of b's prefix that overlaps a's suffix, 0 if none.

    Checks decreasing overlap lengths so returns the longest valid overlap.

    Two-stage: (1) an exact k-mer match at a-position p_a / b-position p_b
    can only be part of a valid suffix/prefix alignment at exactly
    ov = len(a) - p_a + p_b -- profiling on real 16S data showed the
    previous blind "try every ov from limit down to min_ov" scan was 83%
    of total assemble_umi() runtime (200M+ character comparisons for
    227K calls), almost all wasted on reads that share a coincidental
    k-mer near the boundary but have no real overlap. Deriving candidate
    ov values directly from where the shared k-mer actually sits collapses
    that scan to just the handful of lengths real matches imply.
    (2) if none of those candidates pass max_mm, falls back to the
    original exhaustive scan over every remaining ov -- this is only a
    safety net (dense/evenly-spaced mismatches could in principle break
    up every exact k-mer window within the true best ov while an
    unrelated coincidental k-mer elsewhere still passes the existence
    pre-filter below), so stage (1) can only make this function faster,
    never change what it returns relative to before.
    """
    limit = min(len(a), len(b))
    if limit < min_ov:
        return 0

    # existence pre-filter, same semantics as before: no shared k-mer
    # anywhere in the boundary window -> definitely no valid overlap
    check_len = min(limit, max(min_ov * 3, seed_k * 4))
    a_kmers = _kmer_positions(a, seed_k, start=len(a) - check_len)
    b_kmers = _kmer_positions(b, seed_k, end=check_len)
    shared = a_kmers.keys() & b_kmers.keys()
    if not shared:
        return 0

    # stage 1: only test ov values an actual k-mer match implies
    candidates = set()
    for kmer in shared:
        for p_a in a_kmers[kmer]:
            for p_b in b_kmers[kmer]:
                ov = len(a) - p_a + p_b
                if min_ov <= ov <= limit:
                    candidates.add(ov)

    for ov in sorted(candidates, reverse=True):
        if _within_mismatch_budget(a[-ov:], b[:ov], ov, max_mm):
            return ov

    # stage 2 (rare safety net): fall back to the exhaustive scan
    for ov in range(limit, min_ov - 1, -1):
        if ov in candidates:
            continue  # already tested in stage 1
        if _within_mismatch_budget(a[-ov:], b[:ov], ov, max_mm):
            return ov
    return 0


# ── assembler ─────────────────────────────────────────────────────────────────

def _build_pool_kmer_index(pool, seed_k):
    """
    One-time k-mer index over an entire UMI's read pool (both forward and
    reverse-complement orientation of every read), built ONCE per
    assemble_umi() call and reused by the internal-anchor fallback across
    every contig-building attempt for that UMI.

    Without this, the fallback would rescan every remaining candidate's
    full length from scratch on every single stall -- expensive
    (profiled at ~72% of total runtime on real 16S data) precisely
    because it's real data with lots of boundary-extension failures,
    i.e. exactly the reads this fallback exists to rescue. Building one
    index up front turns each fallback lookup into O(check_len) dict
    lookups instead of O(pool_size * read_length).

    Maps kmer -> list of (pool_index, position_in_variant, is_rc).
    `position_in_variant` indexes into pool[pool_index] if is_rc is
    False, or rc(pool[pool_index]) if is_rc is True.
    """
    index = defaultdict(list)
    for idx, seq in enumerate(pool):
        for variant, is_rc in ((seq, False), (rc(seq), True)):
            for p in range(len(variant) - seed_k + 1):
                index[variant[p:p + seed_k]].append((idx, p, is_rc))
    return index


def _build_boundary_kmer_index(pool, min_ov, seed_k):
    """
    Build a one-time inverted index for ordinary suffix/prefix overlap.

    ``suffix_prefix_overlap`` first requires one shared ``seed_k``-mer
    between the two boundary windows. Indexing those same windows lets a
    contig query only reads that can pass that mandatory pre-filter instead
    of scanning every remaining read on every raw-contig attempt.

    Maps boundary k-mer -> pool indices. Both orientations are indexed, but
    orientation and the final overlap are still checked by
    ``suffix_prefix_overlap`` so this is candidate generation only.
    """
    window = max(min_ov * 3, seed_k * 4)
    prefix = defaultdict(set)
    suffix = defaultdict(set)
    for idx, seq in enumerate(pool):
        for variant in (seq, rc(seq)):
            check_len = min(len(variant), window)
            if check_len < seed_k:
                continue
            for kmer in _kmer_positions(variant, seed_k, end=check_len):
                prefix[kmer].add(idx)
            for kmer in _kmer_positions(
                    variant, seed_k, start=len(variant) - check_len):
                suffix[kmer].add(idx)
    return prefix, suffix


def _boundary_overlap_candidates(contig, unused_set, boundary_index,
                                 min_ov, seed_k):
    """
    Return unused reads that can pass suffix_prefix_overlap's mandatory
    shared-k-mer pre-filter at either end of ``contig``.

    This deliberately returns a superset: final mismatch/overlap validation
    remains unchanged. A candidate shorter than the standard boundary window
    may make the queried contig window wider than strictly necessary, which
    can only add false-positive work, never hide a valid overlap.
    """
    if boundary_index is None:
        return unused_set

    window = min(len(contig), max(min_ov * 3, seed_k * 4))
    if window < seed_k:
        return set()

    prefix, suffix = boundary_index
    candidates = set()
    for kmer in _kmer_positions(
            contig, seed_k, start=len(contig) - window):
        candidates.update(prefix.get(kmer, ()))
    for kmer in _kmer_positions(contig, seed_k, end=window):
        candidates.update(suffix.get(kmer, ()))
    return candidates & unused_set


# A stalled OLC pass needs a stronger signal than the 10-mer pre-filter used
# for fast boundary-overlap rejection.  These are deliberately internal
# constants, not another user-facing tuning surface: the rescue is only for
# dense, error-bearing UMI read pools, where several longer exact anchors are
# expected even though the full raw-read overlap exceeds max_mm.
_COLLECTIVE_ANCHOR_K = 17
_COLLECTIVE_MIN_ANCHORS = 3
_COLLECTIVE_MIN_ANCHOR_SPAN = 40
_COLLECTIVE_MIN_EXTENSION_SUPPORT = 2
_COLLECTIVE_MAX_KMER_HITS = 12
_MAX_RAW_CONTIGS_FOR_WIDE_MERGE = 50
_MAX_REVERSE_ANCHOR_CANDIDATES = 150
# Cross-attempt evidence is valuable for ordinary UMI depths but repeatedly
# rescanning every read becomes pathological for rare, extremely deep barcode
# groups. Candidate-only collective rescue remains enabled above this bound.
# The real 20k slice has median=37, P99=130 reads/UMI; 120 covers 98.6% of UMIs
# and 448/449 of the diagnosed bridged_one_end cohort.
_COLLECTIVE_MAX_CROSS_ATTEMPT_READS = 120


def _collective_anchor_extend(contig, pool, anchor_index, candidate_set,
                              evidence_set=None):
    """
    Rescue a stalled greedy extension using collective, not pairwise, evidence.

    Each remaining read/orientation is placed by the offset supported by its
    longest co-linear set of exact 17-mer anchors in ``contig``.  A raw read is
    never accepted merely because its overall mismatch rate is low: independent
    sequencing errors make that test reject genuine overlaps at realistic error
    rates.  Instead, all strongly placed reads vote in a temporary pileup and
    only bases beyond a contig end with support from at least two reads extend
    it.  This retains the existing conservative behaviour for singleton tips
    and barcode collisions while allowing redundant noisy reads to establish a
    shared extension.

    ``candidate_set`` contains reads still available to the current greedy
    attempt; only those reads may be consumed. ``evidence_set`` may be wider:
    assemble_umi passes every read in the UMI so reads consumed by an earlier
    attempt can still vote for a bridge. This shares evidence across attempts
    without putting raw reads into the post-assembly contig merge pool.

    ``anchor_index`` maps k-mer to (pool_idx, read_pos, is_rc), as built by
    _build_pool_kmer_index(). Repetitive/low-complexity anchors are ignored so
    a conserved short motif cannot by itself join two molecules.
    """
    k = _COLLECTIVE_ANCHOR_K
    n = len(contig)
    if n < k:
        return None, set()
    if evidence_set is None:
        evidence_set = candidate_set

    # (read index, orientation, implied offset) -> contig anchor positions.
    # Positions rather than a scalar count let us reject several copies of one
    # repeated motif: trustworthy anchors must also span sequence.
    offset_hits = defaultdict(list)
    for contig_pos in range(n - k + 1):
        kmer = contig[contig_pos:contig_pos + k]
        if len(set(kmer)) < 3:
            continue
        hits = anchor_index.get(kmer, ())
        if len(hits) > _COLLECTIVE_MAX_KMER_HITS:
            continue
        for idx, read_pos, is_rc in hits:
            if idx in evidence_set:
                offset_hits[(idx, is_rc, contig_pos - read_pos)].append(contig_pos)

    placements = []
    best_by_read = {}
    for (idx, is_rc, offset), positions in offset_hits.items():
        unique_positions = sorted(set(positions))
        if len(unique_positions) < _COLLECTIVE_MIN_ANCHORS:
            continue
        span = unique_positions[-1] - unique_positions[0] + k
        if span < _COLLECTIVE_MIN_ANCHOR_SPAN:
            continue

        read = rc(pool[idx]) if is_rc else pool[idx]
        # A fully-contained read may still be useful for polishing, but it
        # must not consume an unused read in an extension rescue.
        if offset >= 0 and offset + len(read) <= n:
            continue
        key = (len(unique_positions), span, -abs(offset), not is_rc)
        previous = best_by_read.get(idx)
        if previous is None or key > previous[0]:
            best_by_read[idx] = (key, offset, read)

    for idx in sorted(best_by_read):
        _, offset, read = best_by_read[idx]
        placements.append((idx, offset, read))
    if len(placements) < _COLLECTIVE_MIN_EXTENSION_SUPPORT:
        return None, set()

    votes = defaultdict(Counter)
    for pos, base in enumerate(contig):
        votes[pos][base] = 1  # preserve the draft on ties inside the contig
    for _, offset, read in placements:
        for read_pos, base in enumerate(read):
            votes[offset + read_pos][base] += 1

    left = 0
    for pos in range(-1, min(votes) - 1, -1):
        if max(votes[pos].values()) < _COLLECTIVE_MIN_EXTENSION_SUPPORT:
            break
        left = pos
    right = n
    for pos in range(n, max(votes) + 1):
        if max(votes[pos].values()) < _COLLECTIVE_MIN_EXTENSION_SUPPORT:
            break
        right = pos + 1
    if left == 0 and right == n:
        return None, set()

    assembled = []
    for pos in range(left, 0):
        counter = votes[pos]
        best_count = max(counter.values())
        assembled.append(min(base for base, count in counter.items()
                             if count == best_count))
    # Collective rescue is an extension mechanism, not a second polisher.
    # Keeping the existing draft byte-identical avoids changing downstream
    # overlap decisions; substitution correction remains polish_contig's job.
    assembled.extend(contig)
    for pos in range(n, right):
        counter = votes[pos]
        best_count = max(counter.values())
        assembled.append(min(base for base, count in counter.items()
                             if count == best_count))
    # Do not consume merely-contained placements.  They were useful evidence
    # while evaluating this rescue, but they add no accepted new boundary and
    # may still seed/extend a separate component if this tentative layout was
    # not the best use of that read.  Only reads that vote on the retained
    # left/right extension become part of this greedy component.
    used = set()
    for idx, offset, read in placements:
        read_end = offset + len(read)
        supports_left = offset < 0 and read_end > left
        supports_right = offset < right and read_end > n
        if idx in candidate_set and (supports_left or supports_right):
            used.add(idx)
    return "".join(assembled), used


# ── join-level spanning-guard instrumentation (P1: observation only) ──────────
#
# Everything below exists to answer one question -- is denovo_junction_qc.py's
# validated post-assembly chimera signal (AUC 0.827: at a real chimera
# junction no genuine read spans it with real margin AND real two-sided
# identity, because no physical molecule contains both sides) ALSO available
# at join time, so a future guard could refuse a chimeric join instead of
# discarding the whole contig afterwards. P1 only measures it -- every
# function here is reachable only when join_map is not None (CLI
# --join-trace / configure(join_trace=...)), which is None/off by default,
# so none of this changes a single assembly decision.
#
# CRITICAL prior failure this must not repeat: an earlier spanning check
# required only that a read PLACE across a position, with no identity check,
# and scored AUC ~= 0.5 -- 16S conserved regions let a read from species A
# anchor onto species B's sequence and appear to span. The two-sided local
# identity test in _bridge_support_at is what makes the signal real; do not
# simplify it away.

_BRIDGE_MARGIN = 30
_BRIDGE_MIN_SIDE_IDENTITY = 0.95
_BRIDGE_WINDOW = 100
# Placement k/thresholds mirror denovo_junction_qc.py / analyze_olc_shortfall.py's
# place()/place_either() exactly (k=17, >=3 anchors spanning >=40bp, hit fan-out
# capped at 12, ambiguous top-2 placements rejected) -- reimplemented locally
# rather than imported (denovo_junction_qc.py reaches its placement helpers via
# a sys.path insert into this package's own test/ directory, which is fine for
# an offline QC script but not a dependency this production assembler should
# carry) rather than copied verbatim.
_BRIDGE_ANCHOR_K = 17
_BRIDGE_MIN_ANCHORS = 3
_BRIDGE_MIN_SPAN = 40
_BRIDGE_MAX_KMER_HITS = 12


def _jt_kmer_index(seq, k):
    index = defaultdict(list)
    for pos in range(len(seq) - k + 1):
        index[seq[pos:pos + k]].append(pos)
    return index


def _jt_place(query, target, target_index, k=_BRIDGE_ANCHOR_K,
              min_anchors=_BRIDGE_MIN_ANCHORS, min_span=_BRIDGE_MIN_SPAN):
    """Mirrors analyze_olc_shortfall.place(): ungapped placement of `query`
    on `target` by its longest co-linear set of exact k-mer anchors, None if
    no single placement dominates (ties on (anchor_count, span) are treated
    as ambiguous, not guessed at)."""
    offsets = defaultdict(list)
    for qpos in range(len(query) - k + 1):
        hits = target_index.get(query[qpos:qpos + k], ())
        if len(hits) > _BRIDGE_MAX_KMER_HITS:
            continue
        for tpos in hits:
            offsets[tpos - qpos].append(qpos)
    candidates = []
    for offset, positions in offsets.items():
        unique = sorted(set(positions))
        if len(unique) >= min_anchors and unique[-1] - unique[0] >= min_span:
            candidates.append((len(unique), unique[-1] - unique[0], offset))
    if not candidates:
        return None
    candidates.sort(reverse=True)
    best = candidates[0]
    if len(candidates) > 1 and candidates[1][:2] == best[:2]:
        return None
    _anchors, _span, start = best
    return start, start + len(query)


def _jt_place_oriented(read, contig, index):
    """Mirrors denovo_junction_qc.place_oriented(): try both orientations,
    require exactly one to place unambiguously."""
    hits = []
    for oriented in (read, rc(read)):
        result = _jt_place(oriented, contig, index)
        if result is not None:
            hits.append((result, oriented))
    if len(hits) != 1:
        return None
    (start, end), oriented = hits[0]
    return start, end, oriented


def _bridge_support_at(contig, junction_pos, evidence_reads,
                        margin=_BRIDGE_MARGIN,
                        min_side_identity=_BRIDGE_MIN_SIDE_IDENTITY,
                        window=_BRIDGE_WINDOW):
    """Count of evidence_reads that verifiably bridge contig[junction_pos]:
    placed (both orientations tried, exactly one unambiguous placement),
    extending >= margin bp past junction_pos on BOTH sides, AND scoring
    >= min_side_identity on a local `window`-bp match profile on each side
    (a windowed test, not whole-side, since placement is ungapped -- a single
    indel anywhere would frameshift and wreck whole-side identity even for an
    otherwise-good read; see denovo_junction_qc.spanning_profile).

    evidence_reads should be the WHOLE UMI read pool (reads already consumed
    by earlier greedy attempts AND still-unused ones) -- a read that could
    bridge a junction has often already been consumed by an earlier boundary
    merge, or rejected there by the pairwise mismatch budget the collective
    rescue exists to work around; restricting evidence to only-unused reads
    would systematically undercount true bridging support.

    This is denovo_junction_qc.spanning_profile specialized to ONE position
    instead of a full per-position profile (this is a per-join query, not a
    whole-contig scan), reimplemented locally -- see the module-level note
    above this function for why not imported.
    """
    n = len(contig)
    if n == 0 or not evidence_reads:
        return 0
    index = _jt_kmer_index(contig, _BRIDGE_ANCHOR_K)
    count = 0
    for read in set(evidence_reads):
        hit = _jt_place_oriented(read, contig, index)
        if hit is None:
            continue
        start, end, oriented = hit
        lo, hi = max(0, start), min(n, end)
        if lo >= hi:
            continue
        if not (lo + margin <= junction_pos < hi - margin):
            continue
        matches = [1 if oriented[pos - start] == contig[pos] else 0
                   for pos in range(lo, hi)]
        k = junction_pos - lo
        m_len = len(matches)
        l_start = max(0, k - window)
        r_end = min(m_len, k + window)
        l_len, r_len = k - l_start, r_end - k
        if l_len < margin or r_len < margin:
            continue
        left_id = sum(matches[l_start:k]) / l_len
        right_id = sum(matches[k:r_end]) / r_len
        if left_id >= min_side_identity and right_id >= min_side_identity:
            count += 1
    return count


_JUNCTION_CONTEXT_HALF = 20


def _junction_context(contig, pos, half=_JUNCTION_CONTEXT_HALF):
    """(lo, context): 40bp sequence context centered on `pos`, captured once
    at join time so _write_join_trace can independently verify the position
    bookkeeping below by re-locating this exact substring in the final
    emitted contig. `lo` (the context's own start offset) is returned
    alongside the string and must be carried -- and shifted in lockstep with
    `pos` -- through every later join this junction survives: near a
    contig's own edge the window is clamped asymmetrically (e.g. pos=4 with
    half=20 clamps lo to 0, not pos-20=-16), and re-deriving lo from a LATER,
    already-shifted pos via the same max(0, pos-half) formula silently gives
    the wrong answer whenever that original clamp fired -- confirmed on real
    16S data: a junction captured near a short draft's edge, then shifted by
    a later collective-rescue prepend, validated against the wrong position
    when lo was recomputed instead of shifted alongside pos."""
    lo = max(0, pos - half)
    hi = min(len(contig), pos + half)
    return lo, contig[lo:hi]


def _flip_junctions(junctions, length):
    """Remap junction records recorded against a pool entry's forward
    orientation onto its reverse-complement of the given `length`: a cut
    just before forward position p sits just after reverse position
    length - p. The frozen context snippet is describing physical bases,
    not a coordinate, so it must itself be reverse-complemented -- and its
    own `lo` recomputed from where THAT flipped snippet starts -- to still
    be a valid re-locatable substring of the RC-oriented sequence."""
    out = []
    for j in junctions:
        ctx = j.get("context", "")
        if ctx:
            new_ctx = rc(ctx)
            new_lo = length - (j["lo"] + len(ctx))
        else:
            new_ctx = ctx
            new_lo = length - j["pos"]
        out.append(dict(j, pos=length - j["pos"], lo=new_lo, context=new_ctx))
    return out


def _context_survives(j, keep_from, keep_to):
    """True iff a junction's ALREADY-CAPTURED context window (tracked via
    its own `lo` field, not re-derived from `pos` -- see _junction_context)
    lies entirely within [keep_from, keep_to) of a join that keeps only that
    sub-range of the source it was captured against. Checking just the
    junction's single position against the keep range is not enough: a
    later truncating join can keep a junction's position while overwriting
    content next to it that its frozen context snippet still depends on,
    silently invalidating that snippet without moving the position at all."""
    return j["lo"] >= keep_from and j["lo"] + len(j["context"]) <= keep_to


def _diff_join_geometry(old_contig, new_contig):
    """
    Reverse-engineer a join's position-bookkeeping shape purely from the
    contig strings before/after the join, so junction tracing never needs
    internal offsets threaded back out of suffix_prefix_overlap's or
    internal_anchor_extend_indexed's existing return contracts (both
    2-tuples, both exercised by test_denovo_seed_olc.py -- changing their
    arity is exactly the kind of "harmless" change this file's history
    warns against). Every join site handled this way has one of two shapes:

    - the old contig survives completely, as a prefix (boundary 3' /
      internal-anchor forward-3') or a suffix (boundary 5' / forward-5') of
      the new one -- "which end grew" is directly len(new) - len(old).
    - only a prefix or suffix of the old contig survives, the rest replaced
      by a candidate's own content (internal-anchor's reverse-3'/reverse-5'
      variants, which trust a candidate over a noisy contig edge -- see
      internal_anchor_extend_indexed's docstring). Detected via longest
      common prefix/suffix between old and new.

    Collective rescue is handled separately at its call site (it can extend
    both ends in one call, which this two-shape model does not cover).

    Returns dict(kind, new_junction_pos, keep_from, keep_to, shift):
    an existing junction at old-contig position p survives iff
    keep_from <= p <= keep_to, and maps to p + shift in the new contig.
    """
    if new_contig.startswith(old_contig):
        return {"kind": "append", "new_junction_pos": len(old_contig),
                "keep_from": 0, "keep_to": len(old_contig), "shift": 0}
    if new_contig.endswith(old_contig):
        shift = len(new_contig) - len(old_contig)
        return {"kind": "prepend", "new_junction_pos": shift,
                "keep_from": 0, "keep_to": len(old_contig), "shift": shift}

    max_cp = min(len(old_contig), len(new_contig))
    cp = 0
    while cp < max_cp and old_contig[cp] == new_contig[cp]:
        cp += 1
    cs = 0
    while cs < max_cp - cp and old_contig[-1 - cs] == new_contig[-1 - cs]:
        cs += 1
    if cp >= cs:
        return {"kind": "truncate_tail", "new_junction_pos": cp,
                "keep_from": 0, "keep_to": cp, "shift": 0}
    drop_before = len(old_contig) - cs
    shift = len(new_contig) - cs - drop_before
    return {"kind": "truncate_head", "new_junction_pos": drop_before + shift,
            "keep_from": drop_before, "keep_to": len(old_contig), "shift": shift}


def _apply_join_shift(junctions, geometry):
    out = []
    for j in junctions:
        if _context_survives(j, geometry["keep_from"], geometry["keep_to"]):
            shift = geometry["shift"]
            out.append(dict(j, pos=j["pos"] + shift, lo=j["lo"] + shift))
    return out


def _candidate_junction_carry(new_contig, geometry, cand_fwd, cand_rc, join_map):
    """Best-effort: detect which orientation of a pool candidate contributed
    a join's new territory and remap that candidate's OWN inherited
    junctions (site-4 contig-merge only -- join_map has no entries for raw
    reads) into new_contig's frame. Detection is by exact substring match
    (the new territory must literally be a suffix/prefix of the candidate in
    one orientation); if neither orientation matches -- should not happen
    given how these joins are constructed, but this is trace bookkeeping,
    not an assembly decision -- carry-over is silently skipped rather than
    guessed at. See _context_survives for why a candidate junction's whole
    context window, not just its position, must fall inside the retained
    sub-range."""
    pos = geometry["new_junction_pos"]
    if geometry["kind"] in ("append", "truncate_tail"):
        territory = new_contig[pos:]
        for cand, is_rc in ((cand_fwd, False), (cand_rc, True)):
            if territory and cand.endswith(territory):
                inherited = join_map.get(cand_fwd, [])
                if is_rc:
                    inherited = _flip_junctions(inherited, len(cand_fwd))
                local_start = len(cand) - len(territory)
                shift = pos - local_start
                return [dict(j, pos=j["pos"] + shift, lo=j["lo"] + shift)
                        for j in inherited
                        if _context_survives(j, local_start, len(cand))]
        return []
    else:
        territory = new_contig[:pos]
        for cand, is_rc in ((cand_fwd, False), (cand_rc, True)):
            if territory and cand.startswith(territory):
                inherited = join_map.get(cand_fwd, [])
                if is_rc:
                    inherited = _flip_junctions(inherited, len(cand_fwd))
                return [dict(j) for j in inherited
                        if _context_survives(j, 0, len(territory))]
        return []


def _make_junction(contig, pos, join_type, ov_len, evidence_reads):
    """Build one junction record, computing its verified bridge-support
    count AT JOIN TIME against `contig` -- the draft as it exists right
    after this specific join, not any later, further-extended state -- plus
    a 40bp context snippet (and its own start offset `lo`) for the
    emission-time position-relocation check (see _junction_context).

    Every call site is required to only reach here with a strictly interior
    position: a junction at or past a contig's own edge can never be
    genuinely spanned by any read (bridge_reads would be forced to 0
    regardless of biology, not because no read exists), and one earlier,
    real bug here (a no-op internal-anchor/boundary "success" -- see
    _record_join's docstring) produced exactly that. This assertion is the
    regression guard for that class of bug: it must never fire once every
    caller correctly skips recording a non-interior position, so it fires
    loudly (this is diagnostic tracing, --join-trace opt-in only -- never
    reachable when join_map is None, i.e. never in the default zero-cost
    path) rather than silently emitting a row that can't mean what its
    columns claim.
    """
    assert 0 < pos < len(contig), (
        "junction position must be strictly interior: pos=%d contig_len=%d "
        "join_type=%s" % (pos, len(contig), join_type))
    bridge = _bridge_support_at(contig, pos, evidence_reads) if evidence_reads else 0
    lo, context = _junction_context(contig, pos)
    return {"pos": pos, "lo": lo, "join_type": join_type, "ov_len": ov_len,
            "bridge_reads": bridge, "context": context}


def _has_min_join_support(contig, junction_pos, evidence_reads, min_join_support):
    """Whether a genuinely new interior join has enough independent support.

    ``_bridge_support_at`` requires a read to place unambiguously and match
    both sides of the join with real margin; it also de-duplicates identical
    sequence strings.  Thus a threshold of two means two distinct raw-read
    sequences, not merely the candidate read plus a duplicated copy.  A
    contained/redundant candidate creates no new interior join and is safe to
    consume without this test.
    """
    if min_join_support <= 1 or not (0 < junction_pos < len(contig)):
        return True
    if evidence_reads is None:
        raise ValueError("join evidence is required when min_join_support > 1")
    return (_bridge_support_at(contig, junction_pos, evidence_reads)
            >= min_join_support)


def _record_join(junctions, join_map, old_contig, new_contig, new_pos, carried,
                  join_type, ov_len, evidence_reads):
    """Shared tail of every join-recording call site in _extend_one_contig:
    build the new junction record, append it (plus any carried-over inherited
    junctions from a candidate that was itself an already-assembled contig)
    and register the result into join_map under new_contig's own string
    value so a later call/the caller can look it up by the final returned
    contig.

    `junctions` must already reflect this join's effect on any PRIOR
    junctions (shifted/dropped per _diff_join_geometry or the boundary-5p
    call site's own equivalent maths) -- this function only adds the new
    junction plus `carried`, it never re-shifts `junctions` itself.

    A join is only recorded as a NEW junction when new_pos is strictly
    interior (0 < new_pos < len(new_contig)): a boundary/internal-anchor
    "success" can still contribute zero new territory when the consumed
    candidate is fully redundant with content already there (verified only
    via the shared max_mm mismatch BUDGET, not exact-match -- a duplicate
    read can pass that budget while replacing a stretch with byte-identical
    or equal-length content). That is correct, harmless assembler behavior
    -- a redundant read consumed, nothing else changes -- but it is not a
    join this trace cares about, and recording it anyway would place
    junction_pos at or past the contig's own edge, where by construction no
    read can ever bridge it (bridge_reads forced to 0 regardless of
    biology) and _write_join_trace's own interior-position assertion would
    reject it. `carried` and any pre-shifted `junctions` are still folded
    in and registered either way, since a candidate's own content can be
    genuinely new even when this specific join creates no fresh boundary.
    """
    result = list(junctions) + list(carried)
    if 0 < new_pos < len(new_contig):
        result.append(_make_junction(new_contig, new_pos, join_type, ov_len, evidence_reads))
    join_map[new_contig] = result
    return result


def _internal_anchor_extend_3prime_indexed(contig, pool, kmer_index, unused_set,
                                            min_ov, max_mm, seed_k, min_verify,
                                            check_len_override=None):
    """
    Index-accelerated 3' (right-end) internal-anchor extension: looks up
    contig's tail k-mers directly in the pre-built pool-wide index
    instead of rescanning every remaining candidate's full length.
    Only considers candidates whose pool index is still in unused_set.

    check_len_override: search window size in place of the default
    max(min_ov*3, seed_k*4) -- see internal_anchor_extend_indexed's
    docstring for why contig-merging needs this to cover the whole
    contig instead of a read-scale edge window.

    Returns (new_contig, used_pool_index) or (None, None).
    """
    n = len(contig)
    requested = check_len_override if check_len_override is not None else max(min_ov * 3, seed_k * 4)
    check_len = min(n, requested)
    if check_len < seed_k:
        return None, None

    tail = contig[-check_len:]
    best = None  # (sort_key, new_contig, pool_idx)

    for p in range(len(tail) - seed_k + 1):
        for idx, j, is_rc in kmer_index.get(tail[p:p + seed_k], ()):
            if idx not in unused_set:
                continue
            candidate = rc(pool[idx]) if is_rc else pool[idx]
            L = len(candidate)
            contig_start = n - check_len + p
            overlap_len = min(check_len - p, L - j)
            if overlap_len < min_verify:
                continue
            if j + overlap_len >= L:
                # candidate contributes NO new trailing content -- its match
                # reaches all the way to its own end, so accepting this would
                # consume the candidate's pool slot for zero gain while
                # silently discarding whatever unique content it has BEFORE
                # the anchor (never used by this direction). Only relevant
                # with a wide check_len_override (contig-merging): read-scale
                # windows rarely let a whole short read get swallowed this
                # way, but a shorter sibling contig can look "fully contained"
                # in exactly this spot and disappear without contributing its
                # own unique edge. Skip so the candidate stays available for
                # a direction that might actually use that edge.
                continue
            contig_region = contig[contig_start:contig_start + overlap_len]
            cand_region = candidate[j:j + overlap_len]
            if _within_mismatch_budget(contig_region, cand_region, overlap_len, max_mm):
                # Selection rule deliberately matches the pre-index
                # implementation's behavior (verified on real data: 97.2%
                # identity to a BLAST-confirmed truth over a 1115bp rescue)
                # rather than a locally "smarter" one: prefer the SMALLEST
                # pool index with any valid anchor (not the biggest single
                # new_seq_len gain) -- greedy assembly is order-sensitive,
                # and taking the single best-looking step here can strand
                # a read that would have combined better with others in a
                # LATER contig-building attempt within the same UMI. Within
                # one candidate, prefer its longest confirmed overlap.
                sort_key = (idx, -overlap_len)
                if best is None or sort_key < best[0]:
                    best = (sort_key, contig + candidate[j + overlap_len:], idx)

    if best is not None:
        return best[1], best[2]
    return None, None


def _internal_anchor_extend_5prime_indexed(contig, pool, kmer_index, unused_set,
                                            min_ov, max_mm, seed_k, min_verify,
                                            check_len_override=None):
    """
    Index-accelerated 5' (left-end) internal-anchor extension: symmetric
    counterpart of the 3' version, looking up contig's head k-mers and
    checking whether a candidate's content BEFORE the anchor can be
    prepended.

    check_len_override: see _internal_anchor_extend_3prime_indexed.

    Returns (new_contig, used_pool_index) or (None, None).
    """
    n = len(contig)
    requested = check_len_override if check_len_override is not None else max(min_ov * 3, seed_k * 4)
    check_len = min(n, requested)
    if check_len < seed_k:
        return None, None

    head = contig[:check_len]
    best = None  # (sort_key, new_contig, pool_idx)

    for p in range(len(head) - seed_k + 1):
        for idx, j, is_rc in kmer_index.get(head[p:p + seed_k], ()):
            if idx not in unused_set:
                continue
            candidate = rc(pool[idx]) if is_rc else pool[idx]
            L = len(candidate)
            cand_prefix_start = j - p
            if cand_prefix_start <= 0:
                # <0: candidate doesn't reach back to contig's own start.
                # ==0: candidate contributes NO new leading content -- see
                # _internal_anchor_extend_3prime_indexed's symmetric check
                # for why accepting this would silently discard whatever
                # unique content the candidate has AFTER the anchor.
                continue
            # unlike the 3' case (where the verified window naturally runs
            # from the anchor p out to contig's end), here it runs from
            # contig's own start (position 0) out to the anchor -- NOT
            # "check_len - p" (that formula belongs to the 3' direction;
            # using it here silently truncated the verified region and
            # caused real anchors to be missed, regressing a real-data
            # rescue from 1115bp down to ~580bp before this fix).
            overlap_len = min(check_len, L - cand_prefix_start)
            if overlap_len < min_verify:
                continue
            cand_region = candidate[cand_prefix_start:cand_prefix_start + overlap_len]
            if len(cand_region) < overlap_len:
                continue
            contig_region = contig[:overlap_len]
            if _within_mismatch_budget(contig_region, cand_region, overlap_len, max_mm):
                # see 3' version: prefer smallest pool index (matches
                # pre-index behavior), not the single biggest gain
                sort_key = (idx, -overlap_len)
                if best is None or sort_key < best[0]:
                    best = (sort_key, candidate[:cand_prefix_start] + contig, idx)

    if best is not None:
        return best[1], best[2]
    return None, None


def _internal_anchor_extend_3prime_reverse_indexed(contig, pool, kmer_index, unused_set,
                                                     min_ov, max_mm, seed_k, min_verify,
                                                     check_len_override=None):
    """
    Reverse-direction mirror of _internal_anchor_extend_3prime_indexed:
    that function trusts the contig's own tail and searches for an anchor
    inside a candidate; this one instead trusts a candidate's own content
    and searches for where it anchors somewhere INSIDE the contig's tail
    region -- so a noisy contig tail (in practice: the very first seed's
    own raw, unverified 3' edge -- see internal_anchor_extend_indexed's
    former KNOWN LIMITATION note) gets truncated and replaced instead of
    being carried forward unverified forever, the way the forward-only
    version does (it only ever appends new candidate content past an
    anchor -- it never touches contig content already committed before
    that anchor).

    Only accepts a candidate that reaches forward far enough to cover the
    contig all the way from the anchor through its current end -- if the
    candidate is shorter than that, there's no way to tell whether the
    contig's own un-covered remainder beyond the candidate is noise or
    genuine, so it's left alone rather than guessed at.

    Uses the shared pool-wide k-mer index to look up candidates from the
    contig's tail. The former implementation inverted this lookup: it built
    a contig index, then rescanned every position of up to 150 candidates on
    every stalled raw-contig attempt. That retained an
    O(attempts * remaining-pool) pathological path after ordinary boundary
    overlap had already been indexed.

    check_len_override: see _internal_anchor_extend_3prime_indexed.

    Returns (new_contig, used_pool_index) or (None, None).
    """
    n = len(contig)
    requested = check_len_override if check_len_override is not None else max(min_ov * 3, seed_k * 4)
    check_len = min(n, requested)
    if check_len < seed_k:
        return None, None

    best = None  # (sort_key, new_contig, pool_idx)

    candidates = (set(sorted(unused_set)[:_MAX_REVERSE_ANCHOR_CANDIDATES])
                  if len(unused_set) > _MAX_REVERSE_ANCHOR_CANDIDATES else unused_set)
    for anchor_pos in range(n - check_len, n - seed_k + 1):
        kmer = contig[anchor_pos:anchor_pos + seed_k]
        for idx, p, is_rc in kmer_index.get(kmer, ()):
            if idx not in candidates:
                continue
            cand = rc(pool[idx]) if is_rc else pool[idx]
            L = len(cand)
            # The anchor can sit anywhere in the candidate; the shared
            # full-pool index retains every position, not just an end window.
            if L - p < n - anchor_pos:
                continue  # candidate doesn't reach forward to contig's own end
            overlap_len = n - anchor_pos
            if overlap_len < min_verify:
                continue
            contig_region = contig[anchor_pos:n]
            cand_region = cand[p:p + overlap_len]
            if _within_mismatch_budget(contig_region, cand_region, overlap_len, max_mm):
                new_contig = contig[:anchor_pos] + cand[p:]
                sort_key = (idx, -overlap_len)
                if best is None or sort_key < best[0]:
                    best = (sort_key, new_contig, idx)

    if best is not None:
        return best[1], best[2]
    return None, None


def _internal_anchor_extend_5prime_reverse_indexed(contig, pool, kmer_index, unused_set,
                                                     min_ov, max_mm, seed_k, min_verify,
                                                     check_len_override=None):
    """
    Reverse-direction mirror of _internal_anchor_extend_5prime_indexed
    (see _internal_anchor_extend_3prime_reverse_indexed for the general
    idea): trusts a candidate's own content and searches for where it
    anchors somewhere INSIDE the contig's HEAD region, so a noisy contig
    head -- in practice: the very first seed's own raw, unverified 5' edge
    -- gets truncated and replaced instead of carried forward unverified.

    Only accepts a candidate that reaches back far enough to cover the
    contig all the way from position 0 through the anchor -- otherwise
    there's no way to tell whether the contig's own un-covered head before
    that point is noise or genuine.

    check_len_override: see _internal_anchor_extend_3prime_indexed.

    Returns (new_contig, used_pool_index) or (None, None).
    """
    n = len(contig)
    requested = check_len_override if check_len_override is not None else max(min_ov * 3, seed_k * 4)
    check_len = min(n, requested)
    if check_len < seed_k:
        return None, None

    best = None  # (sort_key, new_contig, pool_idx)

    # see _internal_anchor_extend_3prime_reverse_indexed for why this is capped.
    candidates = (set(sorted(unused_set)[:_MAX_REVERSE_ANCHOR_CANDIDATES])
                  if len(unused_set) > _MAX_REVERSE_ANCHOR_CANDIDATES else unused_set)
    for anchor_pos in range(check_len - seed_k + 1):
        kmer = contig[anchor_pos:anchor_pos + seed_k]
        for idx, p, is_rc in kmer_index.get(kmer, ()):
            if idx not in candidates:
                continue
            cand = rc(pool[idx]) if is_rc else pool[idx]
            L = len(cand)
            # The anchor can sit anywhere in the candidate; the shared
            # full-pool index retains every position, not just an end window.
            cand_prefix_start = p - anchor_pos
            if cand_prefix_start < 0:
                continue  # candidate doesn't reach back to contig's own start
            overlap_len = anchor_pos + seed_k
            if overlap_len < min_verify:
                continue
            contig_region = contig[:overlap_len]
            cand_region = cand[cand_prefix_start:cand_prefix_start + overlap_len]
            if _within_mismatch_budget(contig_region, cand_region, overlap_len, max_mm):
                new_contig = cand[:p + seed_k] + contig[overlap_len:]
                sort_key = (idx, -overlap_len)
                if best is None or sort_key < best[0]:
                    best = (sort_key, new_contig, idx)

    if best is not None:
        return best[1], best[2]
    return None, None


def internal_anchor_extend_indexed(contig, pool, kmer_index, unused_set,
                                    min_ov, max_mm, seed_k=10, min_verify=60,
                                    check_len_override=None):
    """
    Fallback extension tried only once ordinary boundary suffix/prefix
    extension is fully exhausted for a contig. suffix_prefix_overlap only
    ever compares a candidate's literal first/last N bases against the
    contig's boundary -- real reads can have their genuinely matching
    region start partway in (e.g. a short PCR-chimera artifact or
    quality-degraded stretch fused onto an otherwise-accurate read),
    which suffix_prefix_overlap structurally cannot see.

    min_verify: minimum confirmed overlap length to accept -- deliberately
    longer than the default boundary min_ov, since allowing the anchor to
    sit anywhere in a read is more permissive about WHERE a match can
    start; a longer required confirmed stretch guards against spurious
    short matches (e.g. a universally-conserved primer region shared by
    unrelated templates, not a genuine single-molecule overlap).

    check_len_override: search-window size in place of the default
    max(min_ov*3, seed_k*4). Left at the default (None) for read-scale
    calls during initial assembly, where a noisy edge is expected to sit
    only a short distance in from the literal boundary. Passed as a large
    value (effectively "the whole contig") by _dedupe_and_merge_contigs,
    which reuses this same fallback to merge already-assembled contigs --
    real production data showed two contigs' true connection point is
    routinely hundreds of bp deep into one or both of them (not close to
    either edge at all, unlike the read-noise case this was originally
    built for), which the small default window structurally cannot reach.
    Searching the full contig is affordable here because there are only
    a handful of contigs per UMI (<=4), unlike the full read pool.

    Tries, in order: forward 3' (trust contig's tail, search inside
    candidates), forward 5' (trust contig's head, search inside
    candidates), then -- only if both of those fail -- reverse 3' and
    reverse 5' (trust a CANDIDATE's own content instead, search for where
    it anchors inside the contig itself). The reverse pair fixes what
    used to be a known gap: the forward-only pair assumes the contig's
    own boundary is reliable and never corrects it, which mostly doesn't
    matter (a contig's boundary is either the original longest raw read
    or the product of a prior verified merge) except for the very first
    seed's own raw, unverified edge. Ordering forward before reverse
    keeps existing behavior byte-identical whenever a forward match
    exists; reverse only ever fires as a true last resort.

    Returns (new_contig, used_pool_index) or (None, None).
    """
    result = _internal_anchor_extend_3prime_indexed(
        contig, pool, kmer_index, unused_set, min_ov, max_mm, seed_k, min_verify,
        check_len_override)
    if result[0] is not None:
        return result
    result = _internal_anchor_extend_5prime_indexed(
        contig, pool, kmer_index, unused_set, min_ov, max_mm, seed_k, min_verify,
        check_len_override)
    if result[0] is not None:
        return result
    result = _internal_anchor_extend_3prime_reverse_indexed(
        contig, pool, kmer_index, unused_set, min_ov, max_mm, seed_k, min_verify,
        check_len_override)
    if result[0] is not None:
        return result
    return _internal_anchor_extend_5prime_reverse_indexed(
        contig, pool, kmer_index, unused_set, min_ov, max_mm, seed_k, min_verify,
        check_len_override)


def _extend_one_contig(pool, available, min_ov, max_mm, seed_k, use_internal_anchor=True,
                       internal_min_verify=60, kmer_index_holder=None,
                       internal_check_len=None, collective_index_holder=None,
                       collective_evidence_set=None, boundary_index=None,
                       max_contig_len=None, join_map=None, join_evidence_reads=None,
                       is_merge=False, min_join_support=1):
    """
    Build a single greedy-extended contig from the longest still-available
    read in `pool` (a FIXED list of sequences that is never reindexed --
    see assemble_umi for why). `available` is the set of indices into
    `pool` this attempt may draw from; since pool is globally sorted
    longest-first, min(available) is always the longest remaining read.

    Returns (contig, used_indices): used_indices is the subset of
    `available` this attempt consumed (always includes at least the
    seed's own index, even on a failed/orphan attempt, so the caller can
    still make progress).

    Boundary suffix/prefix extension (suffix_prefix_overlap) is always
    tried first. Internal-anchor extension (internal_anchor_extend_indexed)
    is a fallback tried only once a full sweep finds no more boundary
    extensions -- and after every internal-anchor success, boundary
    extension is retried first again before falling back further, since
    the newly-extended boundary may unlock ordinary merges.

    kmer_index_holder: a 1-element list acting as a lazy, shared cache for
    the pool-wide k-mer index (built via _build_pool_kmer_index). Building
    it costs O(pool_size * read_length) -- worth avoiding entirely for
    barcodes where boundary extension already resolves everything and the
    fallback never triggers. The index is built on the first actual
    fallback attempt (across possibly several _extend_one_contig calls
    for the same UMI) and reused after that. Required (non-None) when
    use_internal_anchor is True.

    internal_check_len: forwarded to internal_anchor_extend_indexed's
    check_len_override -- see that function's docstring.

    collective_index_holder: lazy shared index for the multi-read rescue.
    None disables that rescue (used for the small post-assembly contig merge).

    collective_evidence_set: reads allowed to vote in collective rescue. The
    top-level assembler supplies the full UMI pool so evidence is shared across
    greedy attempts; only ``unused`` reads can still be consumed.

    boundary_index: one-time index from boundary k-mers to pool indices. It
    removes reads that cannot pass suffix_prefix_overlap's existing mandatory
    pre-filter; overlap validation and deterministic candidate order are
    unchanged. None retains the legacy full-pool scan.

    max_contig_len: stop extending once the contig reaches this length,
    returning whatever was built so far, instead of continuing until
    `unused` is exhausted. None (default) is unbounded -- the legacy
    behavior. Backstops against pathological input: a near-monobase read
    pool (e.g. a corrupted barcode's reads, or any other source of
    homopolymer-heavy sequence) gives boundary/internal-anchor overlap
    detection spurious matches at nearly every offset, so the loop below
    can keep finding a "valid" next read indefinitely. A production
    barcode this shape produced a 127kb single contig for what should be
    a <=1.6kb amplicon (denovo.md sec 60) -- filtering likely-corrupted
    barcodes upstream catches most of this class, but not all of it (that
    same incident had two contigs from barcodes with no detectable string
    corruption), so this is the backstop that doesn't depend on
    recognizing the cause in advance.

    join_map: None (default, zero cost) disables ALL join-level bookkeeping
    below -- observation only, see the "join-level spanning-guard
    instrumentation" section above internal_anchor_extend_indexed. When
    provided, a dict from contig string -> list of junction records valid in
    THAT contig's own coordinate frame; this call reads pool[seed_idx]'s
    entry as the seed's inherited junctions (non-empty only when pool holds
    already-assembled contigs, i.e. the site-4 contig-contig merge in
    _dedupe_and_merge_contigs) and registers the working contig's own
    junction list back into join_map under its current string value after
    every successful join, so later calls / the caller can look it up by the
    final returned contig string.

    join_evidence_reads: the whole UMI's raw read pool, used to compute each
    new junction's verified bridge-read count for join tracing and, when
    min_join_support > 1, for acceptance itself.

    min_join_support: require this many verified, distinct raw reads across a
    newly created interior join. The default of 1 preserves legacy output;
    values >1 are an explicit experimental chimera guard.

    is_merge: True labels every join recorded during this call "contig_merge"
    regardless of which sub-mechanism (boundary/internal-anchor) found it --
    matches the trace schema's join_type, which treats "two contigs merged"
    as one category independent of how the connection was found. Only
    meaningful when join_map is not None.
    """
    seed_idx = min(available)
    contig = pool[seed_idx]
    used = {seed_idx}
    unused = set(available) - {seed_idx}
    junctions = list(join_map.get(contig, [])) if join_map is not None else None

    changed = True
    while changed:
        if max_contig_len is not None and len(contig) > max_contig_len:
            break
        changed = False
        # iterate in ascending index order (== pool's longest-first sort
        # order): sets don't preserve insertion order, and trying
        # candidates in a different order than before changes which
        # merge happens first in this greedy algorithm -- which can change
        # the final assembled contig even when every individual merge is
        # independently valid. sorted() restores the original, deterministic
        # "prefer the longest remaining read" trial order.
        boundary_candidates = _boundary_overlap_candidates(
            contig, unused, boundary_index, min_ov, seed_k)
        for i in sorted(boundary_candidates):
            seq = pool[i]
            extended = False

            for cand, cand_is_rc in ((seq, False), (rc(seq), True)):
                # try extending contig at 3' end
                ov = suffix_prefix_overlap(contig, cand, min_ov, max_mm, seed_k)
                if ov:
                    proposed = contig + cand[ov:]
                    if not _has_min_join_support(
                            proposed, len(contig), join_evidence_reads,
                            min_join_support):
                        continue
                    old_contig = contig
                    contig = proposed
                    unused.remove(i)
                    used.add(i)
                    changed = True
                    extended = True
                    # ov can legitimately equal len(cand) (candidate fully
                    # contained in contig's own tail, e.g. a short redundant
                    # read) -- cand[ov:] is then empty and this "extension"
                    # adds no new territory (new_pos would equal len(contig),
                    # past the contig's own last interior offset). Real,
                    # harmless assembler behavior (a redundant read
                    # consumed) -- _record_join's own interior-position
                    # check skips recording a bogus junction for it.
                    if join_map is not None:
                        junctions = _record_join(
                            junctions, join_map, old_contig, contig,
                            new_pos=len(old_contig),
                            carried=_candidate_junction_carry(
                                contig,
                                {"kind": "append", "new_junction_pos": len(old_contig)},
                                seq, rc(seq), join_map),
                            join_type="contig_merge" if is_merge else "boundary_3p",
                            ov_len=ov, evidence_reads=join_evidence_reads)
                    break

                # try extending contig at 5' end
                ov = suffix_prefix_overlap(cand, contig, min_ov, max_mm, seed_k)
                if ov:
                    shift = len(cand) - ov
                    proposed = cand + contig[ov:]
                    if not _has_min_join_support(
                            proposed, shift, join_evidence_reads,
                            min_join_support):
                        continue
                    old_contig = contig
                    contig = proposed
                    unused.remove(i)
                    used.add(i)
                    changed = True
                    extended = True
                    if join_map is not None:
                        # existing junctions: old_contig[ov:] survives (see
                        # boundary 5' case in _diff_join_geometry's docstring)
                        # -- keep only those whose full 40bp context also
                        # survives, not just their bare position (see
                        # _context_survives).
                        kept = [dict(j, pos=j["pos"] + shift, lo=j["lo"] + shift)
                                for j in junctions
                                if _context_survives(j, ov, len(old_contig))]
                        inherited = join_map.get(seq, [])
                        if cand_is_rc:
                            inherited = _flip_junctions(inherited, len(seq))
                        junctions = kept + [dict(j) for j in inherited]
                        # shift == 0 (ov == len(cand)) means cand is fully
                        # contained in contig's own head -- a real, harmless
                        # consumed-redundant-read extension with no new
                        # territory (new_pos would be 0). _record_join's own
                        # interior-position check (0 < new_pos < len(new))
                        # skips recording a bogus junction for that case
                        # while still folding in and registering the
                        # existing/carried junctions above, which legitimately
                        # need updating regardless -- cand may differ from
                        # old_contig[:ov] within the mismatch budget, so its
                        # own content is still genuinely now part of the draft.
                        junctions = _record_join(
                            junctions, join_map, old_contig, contig,
                            new_pos=shift, carried=[],
                            join_type="contig_merge" if is_merge else "boundary_5p",
                            ov_len=ov, evidence_reads=join_evidence_reads)
                    break

            if extended:
                # restart scan so new contig ends are retried against all unused
                break

        if changed or not unused:
            continue

        # boundary extension is fully exhausted -- try the internal-anchor
        # fallback once before giving up on this contig. build the shared
        # index lazily, only now that it's actually needed.
        if use_internal_anchor and unused:
            if kmer_index_holder[0] is None:
                kmer_index_holder[0] = _build_pool_kmer_index(pool, seed_k)
            new_contig, used_idx = internal_anchor_extend_indexed(
                contig, pool, kmer_index_holder[0], unused, min_ov, max_mm, seed_k,
                min_verify=internal_min_verify, check_len_override=internal_check_len)
            if new_contig is not None:
                old_contig = contig
                geometry = _diff_join_geometry(old_contig, new_contig)
                if not _has_min_join_support(
                        new_contig, geometry["new_junction_pos"],
                        join_evidence_reads, min_join_support):
                    new_contig = None
                if new_contig is None:
                    pass
                else:
                    contig = new_contig
                    unused.remove(used_idx)
                    used.add(used_idx)
                    changed = True
                # The reverse-3'/reverse-5' variants (see their docstrings)
                # replace a stretch of contig with a candidate's own content
                # verified only via the shared max_mm mismatch BUDGET, not
                # zero-difference equality -- a candidate that happens to
                # match that stretch exactly (a redundant read genuinely
                # duplicating already-correct content, common in a
                # redundant UMI pool) still "succeeds" and consumes a pool
                # index, but contributes NO new sequence at all. Confirmed
                # on real 16S data (barcode AAAAAATGTATATGT): three separate
                # redundant reads each producing a byte-identical contig
                # this way in succession. That is correct, harmless
                # assembler behavior (a redundant read consumed, nothing
                # else changes) -- but it is not a join in the sense this
                # trace cares about (no new boundary, no chimeric-merge risk
                # was created), and _diff_join_geometry's startswith check
                # trivially classifies old==new as "append" with
                # new_junction_pos == len(contig) -- a position past the
                # contig's own last valid interior offset. Recording it
                # produced exactly that (junction_pos == contig_len,
                # bridge_reads forced to 0 since nothing can span past a
                # contig's own end) and, since the outer while loop retries
                # internal-anchor again next iteration, could repeat for
                # every further redundant candidate -- the duplicate-row
                # source this was traced back to. Guard on real content
                # change instead of trusting internal_anchor_extend_indexed's
                # success return alone.
                if join_map is not None and contig != old_contig:
                    cand_fwd = pool[used_idx]
                    carried = _candidate_junction_carry(
                        contig, geometry, cand_fwd, rc(cand_fwd), join_map)
                    junctions = _record_join(
                        _apply_join_shift(junctions, geometry), join_map,
                        old_contig, contig, geometry["new_junction_pos"], carried,
                        join_type="contig_merge" if is_merge else "internal_anchor",
                        ov_len=-1, evidence_reads=join_evidence_reads)
                if changed:
                    continue

        # Pairwise raw-read mismatch is the wrong evidence model once a
        # redundant UMI contains independent sequencing errors: it combines
        # both reads' errors and can reject every genuine bridge.  The
        # collective rescue is deliberately *after* the existing single-read
        # fallback, preserving a previously accepted conservative extension
        # whenever it is available.  It is lazy and only runs after both
        # normal paths have stalled.
        if collective_index_holder is not None:
            if collective_index_holder[0] is None:
                collective_index_holder[0] = _build_pool_kmer_index(
                    pool, _COLLECTIVE_ANCHOR_K)
            # Preserve the previously accepted candidate-only path exactly.
            # Cross-attempt evidence is a fallback, not a competing vote pool:
            # mixing both up front can change a rescue that already worked and
            # steer the later greedy partition down a different path.
            new_contig, rescue_used = _collective_anchor_extend(
                contig, pool, collective_index_holder[0], unused)
            if (new_contig is None and collective_evidence_set is not None
                    and collective_evidence_set != unused):
                new_contig, rescue_used = _collective_anchor_extend(
                    contig, pool, collective_index_holder[0], unused,
                    evidence_set=collective_evidence_set)
            # Every accepted greedy step must consume at least one currently
            # unused read. Cross-attempt evidence may vote, but an
            # evidence-only extension has no new candidate to claim and can
            # be replayed at another repetitive offset forever. A real
            # 28-read UMI reproduced exactly that cycle: rescue_used stayed
            # empty while the contig grew by alternating 130/352 bp steps.
            # Requiring progress both prevents the loop and preserves the
            # candidate/evidence separation contract.
            if new_contig is not None and rescue_used:
                old_contig = contig
                contig = new_contig
                unused -= rescue_used
                used |= rescue_used
                changed = True
                if join_map is not None:
                    # _collective_anchor_extend always retains the prior
                    # draft byte-identical in the middle (assembled.extend
                    # (contig) in its own source -- it is an extension
                    # mechanism, not a second polisher) and may grow either
                    # or both ends in this one call, unlike every other join
                    # site here -- so it needs its own bookkeeping rather
                    # than _diff_join_geometry's single-new-junction model.
                    prepend_len = contig.find(old_contig)
                    if prepend_len < 0:
                        # Invariant above broke; every existing junction's
                        # pos/lo is only valid in old_contig's frame, so
                        # registering them unshifted under the NEW contig's
                        # key would be actively wrong, not just incomplete --
                        # leave join_map untouched (a later lookup on this
                        # contig value correctly returns no junctions) rather
                        # than guess at a position. This is observation, not
                        # an assembly decision.
                        pass
                    else:
                        append_len = len(contig) - prepend_len - len(old_contig)
                        jtype = "contig_merge" if is_merge else "collective"
                        shifted = [dict(j, pos=j["pos"] + prepend_len, lo=j["lo"] + prepend_len)
                                   for j in junctions]
                        new_js = []
                        if prepend_len > 0:
                            new_js.append(_make_junction(
                                contig, prepend_len, jtype, -1, join_evidence_reads))
                        if append_len > 0:
                            new_js.append(_make_junction(
                                contig, prepend_len + len(old_contig), jtype, -1,
                                join_evidence_reads))
                        junctions = shifted + new_js
                        join_map[contig] = junctions

    return contig, used


def _assemble_umi_once(seqs, min_ov=20, max_mm=0.05, min_ctg=400, seed_k=10,
                       use_internal_anchor=True, internal_min_verify=60,
                       use_collective_rescue=True, use_minimizer_dedup=False,
                       use_cross_attempt_evidence=True, join_map=None,
                       min_join_support=1):
    """
    Greedy seed-extension assembly for one UMI's reads.

    A barcode can genuinely carry reads from more than one physical DNA
    fragment (barcode reuse/collision is a known stLFR/cLFR reality) --
    megahit's de Bruijn graph naturally splits into multiple connected
    components in that case, one contig per fragment. This loops the
    same greedy seed-extension over whatever reads are left after each
    contig, so it does the same instead of silently discarding every
    read that didn't merge into the first (longest) seed's chain.

    seqs        : list of DNA sequences (forward strand, no quality)
    min_ov      : minimum overlap length to merge two reads [20]
    max_mm      : maximum mismatch rate in overlap region [0.05]
    min_ctg     : discard contigs shorter than this [400]
    seed_k      : k-mer length for overlap pre-filter [10]
    use_internal_anchor : also try internal_anchor_extend_indexed() as a fallback
                  when boundary suffix/prefix extension stalls [True]
    internal_min_verify : min confirmed overlap length for the internal
                  anchor fallback to accept a match [60]
    use_collective_rescue : on a stalled OLC pass, jointly place multiple
                  reads by long k-mer anchors and extend only their supported
                  pileup consensus [True]
    use_minimizer_dedup : collapse near-identical reads (same canonical
                  minimizer) down to _DEDUP_KEEP_PER_CLUSTER before assembly
                  -- OFF by default: confirmed to collapse genuinely-distinct
                  reads that share a conserved motif on real 16S data
                  (see _canonical_minimizer's docstring) [False]
    use_cross_attempt_evidence : when a stalled contig's collective rescue
                  (use_collective_rescue) finds no support among the current
                  attempt's own leftover reads, also let it vote using reads
                  already consumed by EARLIER raw-building attempts for this
                  same UMI (only up to _COLLECTIVE_MAX_CROSS_ATTEMPT_READS
                  total pool size). True reproduces the wider-evidence
                  behaviour that rescued two known hard cases (580->1433bp,
                  800->1121bp); False restricts every rescue to
                  candidate-only evidence (the original, narrower semantics).
                  A full 20k run with this True (mixed in with other changes)
                  showed both lengthened and shortened main contigs relative
                  to the old baseline (see lfr.md 21) -- not a clean isolation.
                  This toggle exists so a dedicated candidate-only 20k A/B can
                  attribute that effect to cross-attempt evidence specifically
                  instead of conflating it with the unrelated endpoint-index
                  speed fix [True]
    join_map    : None (default, zero cost -- see the "join-level
                  spanning-guard instrumentation" section above
                  internal_anchor_extend_indexed) or a dict from contig
                  string -> junction-record list, threaded into every
                  _extend_one_contig call (raw-building AND the post-merge
                  call inside _dedupe_and_merge_contigs) so joins are
                  traceable end to end. Bridge-support evidence is always
                  the UMI's WHOLE original read pool (this function's own
                  `seqs` argument, captured before any minimizer dedup --
                  a read that could bridge a junction is often one already
                  consumed by an earlier attempt, or rejected there by the
                  pairwise mismatch budget, so restricting evidence to
                  fewer/deduped reads would undercount).
    min_join_support : verified distinct raw reads required at each newly
                  created interior join; 1 preserves legacy behavior.

    Returns list of contig sequences (0 or more per UMI).
    """
    if not seqs:
        return []

    join_evidence_reads = (seqs if join_map is not None or min_join_support > 1
                           else None)

    if use_minimizer_dedup:
        seqs = _minimizer_dedupe(seqs)

    # deduplicate, sort longest-first -> longest read is seed.
    # Secondary sort key (the string itself) makes tie-breaking deterministic:
    # set() iteration order depends on Python's per-process string hash
    # randomization, so without this, inputs with many same-length reads
    # (e.g. fixed-length test data) would pick a different, effectively
    # random seed read -- and thus a different assembly -- on every rerun.
    # `pool` is FIXED (never reindexed) so kmer_index's pool-indices stay
    # valid across every contig-building attempt below; which reads are
    # still up for grabs is tracked separately via the shrinking
    # `available` index set instead of physically re-slicing pool.
    pool = sorted(set(seqs), key=lambda s: (-len(s), s))

    # lazy, shared across every contig-building attempt below: built on
    # the first actual fallback trigger, not unconditionally up front --
    # see _extend_one_contig's kmer_index_holder docstring for why (skips
    # the O(pool_size * read_length) index-build cost entirely for
    # barcodes where boundary extension already resolves everything).
    kmer_index_holder = [None]
    collective_index_holder = [None] if use_collective_rescue else None
    collective_evidence_set = (
        set(range(len(pool)))
        if use_cross_attempt_evidence and len(pool) <= _COLLECTIVE_MAX_CROSS_ATTEMPT_READS
        else None
    )
    boundary_index = _build_boundary_kmer_index(pool, min_ov, seed_k)

    available = set(range(len(pool)))
    raw_contigs = []
    # Do not cap this pre-merge phase's LENGTH TARGET (min_ctg) here.  A greedy
    # attempt may strand a valid component until the fifth or later seed;
    # applying the user-visible final output floor here would prevent
    # _dedupe_and_merge_contigs from seeing it. max_contig_len below is a
    # different thing -- a runaway-growth backstop, not a target -- so it
    # applies at this raw stage where the pathological growth actually
    # happens (denovo.md sec 60).
    max_contig_len = max(min_ctg * 5, 3000)
    while available:
        contig, used = _extend_one_contig(pool, available, min_ov, max_mm, seed_k,
                                          use_internal_anchor, internal_min_verify,
                                          kmer_index_holder=kmer_index_holder,
                                          collective_index_holder=collective_index_holder,
                                          collective_evidence_set=collective_evidence_set,
                                          boundary_index=boundary_index,
                                          max_contig_len=max_contig_len,
                                          join_map=join_map,
                                          join_evidence_reads=join_evidence_reads,
                                          min_join_support=min_join_support)
        # collect every attempt regardless of length -- filtering by min_ctg
        # here (before _dedupe_and_merge_contigs runs) would silently
        # discard pieces that individually fall short but would merge with
        # another attempt into something that clears the floor. Verified
        # on real production data: contigs that look "separate" routinely
        # share hundreds of bp with each other (see _dedupe_and_merge_contigs),
        # so min_ctg must apply to the POST-merge result, not each raw attempt.
        raw_contigs.append(contig)
        # always drop every read the attempt consumed (even just the seed
        # itself, on a failed/orphan attempt) so available strictly shrinks
        # and a genuinely separate fragment among the rest still gets a shot
        available -= used

    merged = _dedupe_and_merge_contigs(raw_contigs, min_ov, max_mm, seed_k,
                                       use_internal_anchor=use_internal_anchor,
                                       internal_min_verify=internal_min_verify,
                                       max_contig_len=max_contig_len,
                                       join_map=join_map,
                                       join_evidence_reads=join_evidence_reads,
                                       min_join_support=min_join_support)
    return [c for c in merged if len(c) >= min_ctg]


def assemble_umi(seqs, min_ov=20, max_mm=0.05, min_ctg=400, seed_k=10,
                 use_internal_anchor=True, internal_min_verify=60,
                 use_collective_rescue=True, use_minimizer_dedup=False,
                 use_cross_attempt_evidence=True, adaptive_seed_k=False,
                 fallback_seed_k=10, join_map=None, min_join_support=1):
    """Assemble one UMI, optionally retrying a strict k-mer run conservatively.

    With ``adaptive_seed_k=True`` and a primary ``seed_k`` larger than
    ``fallback_seed_k``, the UMI is assembled with the strict k-mer first.  A
    complete fallback run at the smaller k is attempted only when the strict
    run produces no contig meeting ``min_ctg``.  The fallback replaces the
    strict result only when it itself produces a valid contig; this keeps the
    rule from merging two independently chosen candidate sets or selecting a
    longer but unsupported fragment.

    join_map: forwarded to _assemble_umi_once -- None (default) disables all
    join-level tracing/bridge-support computation. When a fallback run
    replaces the primary result, join_map may still carry leftover entries
    for the discarded primary run's contigs; harmless, since only the
    actually-returned contigs' strings are ever looked up.

    min_join_support: minimum number of distinct raw reads that must bridge a
    newly created interior join. Defaults to 1 to preserve legacy behavior.
    """
    result = _assemble_umi_once(
        seqs, min_ov, max_mm, min_ctg, seed_k, use_internal_anchor,
        internal_min_verify, use_collective_rescue, use_minimizer_dedup,
        use_cross_attempt_evidence, join_map=join_map,
        min_join_support=min_join_support)
    if (not adaptive_seed_k or seed_k <= fallback_seed_k
            or any(len(c) >= min_ctg for c in result)):
        return result
    fallback = _assemble_umi_once(
        seqs, min_ov, max_mm, min_ctg, fallback_seed_k, use_internal_anchor,
        internal_min_verify, use_collective_rescue, use_minimizer_dedup,
        use_cross_attempt_evidence, join_map=join_map,
        min_join_support=min_join_support)
    return fallback or result


def _dedupe_and_merge_contigs(contigs, min_ov, max_mm, seed_k,
                              use_internal_anchor=True, internal_min_verify=60,
                              max_contig_len=None, join_map=None,
                              join_evidence_reads=None, min_join_support=1):
    """
    Post-process the contig list from one UMI: the outer loop in
    assemble_umi builds each contig from a single greedy pass, so a read
    that "missed" merging on one pass (e.g. its bridging partner was
    already claimed) can end up starting a second, spurious contig that
    actually belongs to the same fragment as an earlier one. Two cases:

    1. Genuine boundary overlap (one contig's end matches another's
       start) -- these get merged into one longer contig via the same
       suffix/prefix extension logic used for raw reads (a contig is
       just a longer sequence).
    2. A pair's TRUE connection point sits deep inside one or both of
       them, not at either edge -- ordinary suffix_prefix_overlap only
       checks literal boundaries so it can't see this. Verified on real
       production data (a 1.5kb-library run): most multi-contig UMIs had
       200-750bp of genuine shared sequence between their separate
       contigs, sitting hundreds of bp from either edge -- these are
       independently-grown pieces of the same underlying molecule that
       ordinary boundary merging can't reconnect. Uses the same
       internal-anchor fallback built for noisy read edges, but with
       check_len_override effectively covering each WHOLE contig (a
       handful of contigs per UMI, so full-length search is cheap) --
       the default read-scale window is far too narrow to reach a
       connection point hundreds of bp deep into a several-hundred-bp
       contig.

       "A handful of contigs" is the load-bearing assumption: cProfile on
       a real pathological UMI (1065 raw reads, genuinely diverse rather
       than redundant) showed the raw-building phase alone produced 568
       raw contigs, and merging them at full-contig check_len drove
       _internal_anchor_extend_*_reverse_indexed (which rebuilds a k-mer
       index from scratch on every single call, with no caching across
       candidates) to 60%+ of a 140-second single-UMI runtime. Capped at
       _MAX_RAW_CONTIGS_FOR_WIDE_MERGE: past that many raw contigs, the
       UMI is already so fragmented that reconnecting distant pairs via a
       full-contig scan is no longer "a handful of cheap lookups", so this
       falls back to the read-scale default window instead of guessing at
       a merge that costs more than the whole rest of the pipeline.
    3. Pure internal containment (one contig sits entirely inside
       another, not at either edge -- suffix_prefix_overlap only checks
       boundaries so it won't catch this) -- the contained one adds no
       new sequence, so it's simply dropped.

    max_contig_len: forwarded to _extend_one_contig (see its docstring).
    Required here too, not just on the raw-building phase in
    _assemble_umi_once: a pathological pool that hits the cap on every raw
    attempt produces several capped-length raw contigs, and merging those
    back together uncapped defeats the raw-phase cap entirely (a real
    incident barcode did exactly this -- 8 raw attempts at ~3kb each
    re-merged into one 24kb contig before this parameter existed here).

    join_map / join_evidence_reads: forwarded to _extend_one_contig with
    is_merge=True (see its docstring) -- None (default) disables all
    join-level tracing/bridge-support computation for this merge phase.
    Note `contigs` here are themselves already-assembled contigs, not raw
    reads: this is the ONE call site where pool entries can carry their own
    inherited junctions (from the raw-building phase in _assemble_umi_once),
    which is why _extend_one_contig looks candidates' prior junctions up in
    join_map by their own string value rather than needing a second parallel
    list threaded through this dedupe/sort.
    """
    if len(contigs) <= 1:
        return contigs

    # 1+2. merge any pair with a real boundary overlap, or a genuine
    # internal connection point found via the internal-anchor fallback
    # (see case 2 above for why this needs use_internal_anchor=True with
    # a full-contig check_len, unlike the read-level default -- and why
    # that check_len is capped rather than unconditional).
    pool = sorted(set(contigs), key=lambda s: (-len(s), s))
    wide_check_len = (max(len(s) for s in pool)
                      if len(pool) <= _MAX_RAW_CONTIGS_FOR_WIDE_MERGE else None)
    available = set(range(len(pool)))
    kmer_index_holder = [None]
    boundary_index = _build_boundary_kmer_index(pool, min_ov, seed_k)
    merged = []
    while available:
        contig, used = _extend_one_contig(pool, available, min_ov, max_mm, seed_k,
                                          use_internal_anchor=use_internal_anchor,
                                          internal_min_verify=internal_min_verify,
                                          kmer_index_holder=kmer_index_holder,
                                          internal_check_len=wide_check_len,
                                          boundary_index=boundary_index,
                                          max_contig_len=max_contig_len,
                                          join_map=join_map,
                                          join_evidence_reads=join_evidence_reads,
                                          is_merge=True,
                                          min_join_support=min_join_support)
        merged.append(contig)
        available -= used

    # 3. drop pure containment (substring anywhere, not just at a boundary)
    merged.sort(key=len, reverse=True)
    final = []
    for c in merged:
        if not any(c != kept and c in kept for kept in final):
            final.append(c)
    return final


# ── post-assembly polish (majority-vote consensus correction) ─────────────────

def polish_contig(contig, seqs, min_ov=20, max_mm=0.05, seed_k=10,
                   min_coverage=3, vote_concordance=0.6, kmer_step=5):
    """
    Correct substitution errors in a draft contig via per-position majority
    vote across the UMI's reads, borrowing the idea from sparse_denovo_trasm's
    CallArray pileup (see module docstring) without its k-mer-index / numpy /
    homopolymer machinery.

    Re-aligns each read against the contig with a cheap k-mer offset vote
    (no full DP): sample k-mers from the read, look up matching positions in
    a one-time contig k-mer index, and take the most common implied offset.
    Reads/orientations that don't clear max_mm over the resulting overlap are
    skipped. Surviving reads cast one vote per covered position; the contig's
    own base gets 1 implicit vote so a single dissenting read can't flip
    anything. A position flips only when the alternative base reaches
    >= vote_concordance of >= min_coverage total votes.

    Substitution-only, like assemble_umi's overlap check: an indel-bearing
    read (e.g. homopolymer-length error) fails the concordance check for
    everything downstream of the indel and is simply excluded from voting,
    rather than corrupting the consensus.

    Returns the polished contig (same length as input; only base
    substitutions, no indels are ever introduced by this step).
    """
    n = len(contig)
    if n == 0 or not seqs:
        return contig

    k = seed_k
    contig_kmers = defaultdict(list)
    for i in range(0, n - k + 1):
        contig_kmers[contig[i:i + k]].append(i)

    votes = [None] * n  # lazy per-position Counter

    for raw in seqs:
        for cand in (raw, rc(raw)):
            L = len(cand)
            if L < min_ov:
                continue

            offset_counts = {}
            for j in range(0, L - k + 1, kmer_step):
                for pos in contig_kmers.get(cand[j:j + k], ()):
                    off = pos - j
                    offset_counts[off] = offset_counts.get(off, 0) + 1
            if not offset_counts:
                continue

            best_offset = max(offset_counts, key=offset_counts.get)

            read_start = max(0, -best_offset)
            contig_start = max(0, best_offset)
            overlap_len = min(L - read_start, n - contig_start)
            if overlap_len < min_ov:
                continue

            read_region = cand[read_start:read_start + overlap_len]
            contig_region = contig[contig_start:contig_start + overlap_len]
            if not _within_mismatch_budget(contig_region, read_region, overlap_len, max_mm):
                continue

            for idx in range(overlap_len):
                pos = contig_start + idx
                if votes[pos] is None:
                    votes[pos] = Counter()
                votes[pos][read_region[idx]] += 1
            break  # this orientation aligned; don't also try the RC of the same read

    polished = list(contig)
    for i, counter in enumerate(votes):
        if counter is None:
            continue
        counter[contig[i]] += 1  # original base's implicit vote
        total = sum(counter.values())
        if total < min_coverage:
            continue
        base, cnt = counter.most_common(1)[0]
        if base != contig[i] and cnt / total >= vote_concordance:
            polished[i] = base

    return "".join(polished)


# ── optional mappy fast path ──────────────────────────────────────────────────

def _assemble_umi_mappy(seqs, min_ctg):
    """
    Overlap detection via mappy (minimap2 Python bindings).
    Returns contig list, or None if mappy unavailable / fails.

    CONFIRMED PRODUCTION INCIDENT: this path has none of assemble_umi's
    correctness hardening and was NOT the intended default, but silently
    became the actual production behavior anyway because _CFG["use_mappy"]
    used to default to None ("auto-detect and prefer if importable") with
    no CLI flag to override it -- and the production server happened to
    have mappy installed. Every fix validated on this project across many
    sessions (mismatch-rate control, internal-anchor/collective-rescue
    fallbacks, chimeric-merge safety, minimizer-dedup bug fix, endpoint/
    reverse k-mer indices, the infinite-extension progress guard) only
    applies to assemble_umi and was a complete no-op for any run that took
    this path -- confirmed via a real 1.5M-UMI run where numerous barcodes
    matched a known PRE-fix baseline almost exactly.

    Known gaps relative to assemble_umi, not yet closed:
    - No mismatch-rate check on accepted overlaps -- relies entirely on
      minimap2's own ava-sr scoring/heuristics.
    - Always collapses to exactly ONE contig per barcode: the greedy merge
      loop below has no equivalent of assemble_umi's multiple raw-attempt
      loop, so a barcode carrying reads from more than one physical
      fragment (a known stLFR/cLFR reality -- barcode reuse/collision) has
      no chance of being split into separate, correct contigs -- it is
      prone to producing a single chimeric merge instead.
    - No single-read UMI support (returns None below len(seqs) < 2),
      the exact low-depth case assemble_umi was originally built for.
    - No internal-anchor fallback, no collective/cross-attempt rescue for
      redundant low-overlap noisy reads, no post-assembly polish step.

    use_mappy therefore defaults to False everywhere (configure()/_CFG/CLI
    --mappy is opt-in). Do not re-enable by default without giving this
    path the same validation assemble_umi has been through.
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

def _write_contigs(barcode, contigs, out_file, lock, max_contigs):
    """
    Append contigs to final_contigs_{id}.fa.
    Header: >{barcode}>k41_{i}  (first 15 chars = barcode,
    matching denovo_supp.py record.id[:CBC_LEN]; second '>' marks the
    barcode/UMI boundary, matching megahit-path convention).
    """
    if not contigs:
        return
    lines = []
    for i, seq in enumerate(contigs[:max_contigs]):
        lines.append(">{barcode}>k41_{i}\n{seq}\n".format(barcode=barcode, i=i, seq=seq))
    block = "".join(lines)
    with lock:
        with open(out_file, "a") as fh:
            fh.write(block)


_JOIN_TRACE_HEADER = ("barcode\tcontig_idx\tjoin_type\tov_len\tjunction_pos\t"
                      "contig_len\tbridge_reads\tn_evidence_reads\t"
                      "ctx_pos_agrees\tctx_pos_agrees_fuzzy\n")

_FUZZY_MAX_MISMATCHES = 3
_FUZZY_SEARCH_RADIUS = 10


def _fuzzy_context_agrees(final, context, lo, max_mismatches=_FUZZY_MAX_MISMATCHES,
                          radius=_FUZZY_SEARCH_RADIUS):
    """Best-scoring alignment of the recorded context over `final` within
    +/- radius of the tracked position `lo`, tolerating up to
    max_mismatches substitutions -- never indels, since polish_contig (the
    only source of drift between the context's own draft and `final`) is
    documented substitution-only, same length, no indels ever. Returns 1
    iff some offset in that window scores <= max_mismatches.

    This is the discriminator between "bookkeeping put the position in the
    wrong place" (fails at every nearby offset too) and "the position is
    right but polish edited a few bases inside the window" (succeeds at or
    very near the recorded lo, within the mismatch budget) -- deliberately
    independent of ctx_pos_agrees's exact match and of lo itself: it
    re-searches nearby rather than trusting the recorded position outright.
    """
    if not context:
        return 0
    n, L = len(final), len(context)
    if L > n:
        return 0
    best = None
    for cand_lo in range(max(0, lo - radius), min(n - L, lo + radius) + 1):
        window = final[cand_lo:cand_lo + L]
        hd = sum(1 for a, b in zip(window, context) if a != b)
        if best is None or hd < best:
            best = hd
            if best == 0:
                break
    return 1 if best is not None and best <= max_mismatches else 0


def _write_join_trace(barcode, draft_contigs, final_contigs, join_map,
                      n_evidence_reads, out_path, lock, max_contigs):
    """
    Append one TSV row per join that survives into an EMITTED contig (a join
    belonging to a raw/merged contig dropped by containment filtering or the
    max_contigs truncation never appears here) to --join-trace's output.
    Multiprocessing-safe via the same lock pattern as _write_contigs.

    draft_contigs: assemble_umi's own return value, BEFORE _polish_all.
    join_map is keyed by the exact draft contig string each join produced,
    and polish_contig's substitution-only correction (see its docstring --
    same length, no indels, ever) means the post-polish string is almost
    never byte-identical to any join_map key even though junction POSITIONS
    are unaffected by it.

    final_contigs: the POST-polish contigs actually handed to _write_contigs
    -- contig_idx below matches _write_contigs' own [:max_contigs]
    truncation and the emitted header's k41_{i} exactly, since ground truth
    for any downstream consumer of this trace is built on that same index.
    Context relocation (ctx_pos_agrees) is deliberately checked against this
    final, actually-emitted string, not the pre-polish draft: that is the
    real end-to-end validation, including whatever a handful of polish
    substitutions inside a 40bp context window may have done to it.
    """
    if not draft_contigs:
        return
    rows = []
    for i, (draft, final) in enumerate(zip(draft_contigs[:max_contigs],
                                           final_contigs[:max_contigs])):
        for j in join_map.get(draft, []):
            pos = j["pos"]
            # _make_junction already asserts every junction is interior at
            # CREATION time; re-asserting here against the actual emitted
            # contig_len (post-polish, same length as draft -- see
            # polish_contig's docstring) is the end-to-end regression guard:
            # a bug in a LATER shift/carry step could in principle move a
            # once-valid position out of range even though it started fine.
            # This is diagnostic tracing (--join-trace opt-in only), so it
            # fails loudly rather than silently writing a row whose
            # bridge_reads is structurally forced to 0 by construction, not
            # by biology.
            assert 0 < pos < len(final), (
                "emitted junction position must be strictly interior: "
                "barcode=%s contig_idx=%d join_type=%s pos=%d contig_len=%d"
                % (barcode, i, j["join_type"], pos, len(final)))
            # j["lo"] is the context's own tracked start offset -- shifted in
            # lockstep with pos through every join this junction survived
            # (see _junction_context) -- NOT re-derived from pos here, since
            # a junction captured near an earlier, shorter draft's edge has
            # an asymmetrically-clamped window that max(0, pos-half) cannot
            # reconstruct after later shifts.
            found = final.find(j["context"]) if j["context"] else -1
            agrees = 1 if found == j["lo"] else 0
            agrees_fuzzy = _fuzzy_context_agrees(final, j["context"], j["lo"])
            rows.append("\t".join(str(x) for x in (
                barcode, i, j["join_type"], j["ov_len"], pos, len(final),
                j["bridge_reads"], n_evidence_reads, agrees, agrees_fuzzy)) + "\n")
    if not rows:
        return
    block = "".join(rows)
    with lock:
        write_header = (not os.path.exists(out_path)
                        or os.path.getsize(out_path) == 0)
        with open(out_path, "a") as fh:
            if write_header:
                fh.write(_JOIN_TRACE_HEADER)
            fh.write(block)


# ── per-barcode workers (drop-in replacements) ────────────────────────────────

def _seqs_from_meta(meta, barcode):
    """Extract sequence strings from meta_data dict (strips header lines)."""
    entries = meta.get(barcode)
    if not entries:
        return []
    # entries = ['>id0', 'seq0', '>id1', 'seq1', ...]
    return entries[1::2]


def _polish_all(contigs, seqs, min_ov, max_mm, seed_k):
    if not contigs or not _CFG["polish"]:
        return contigs
    return [
        polish_contig(
            c, seqs, min_ov=min_ov, max_mm=max_mm, seed_k=seed_k,
            min_coverage=_CFG["polish_min_coverage"],
            vote_concordance=_CFG["polish_vote_concordance"],
            kmer_step=_CFG["polish_kmer_step"],
        )
        for c in contigs
    ]


def process_barcode_se(barcode, shared_meta_data2, lock):
    """SE drop-in for denovo_clfr_ram.process_barcode_se."""
    min_ctg  = _CFG["min_ctg"]
    min_ov   = _CFG["min_ov"]
    max_mm   = _CFG["max_mm"]
    seed_k   = _CFG["seed_k"]
    out_file = _CFG["out_file"].format(id=_CFG["out_id"])
    use_mp   = _CFG["use_mappy"]
    use_anchor  = _CFG["use_internal_anchor"]
    anchor_verify = _CFG["internal_min_verify"]
    max_contigs = _CFG["max_contigs"]
    use_dedup = _CFG["use_minimizer_dedup"]
    use_cross_evidence = _CFG["use_cross_attempt_evidence"]
    adaptive_seed_k = _CFG["adaptive_seed_k"]
    fallback_seed_k = _CFG["fallback_seed_k"]
    min_join_support = _CFG["min_join_support"]
    join_trace_path = _CFG["join_trace"]

    seqs = _seqs_from_meta(shared_meta_data2, barcode)
    if not seqs:
        return

    # None (the default) keeps assemble_umi on its exact pre-instrumentation
    # path -- every join_map-gated branch it touches is skipped outright, so
    # this flag OFF is byte-identical output, not just "traces nothing".
    join_map = {} if join_trace_path else None

    contigs = None
    if use_mp is not False:
        contigs = _assemble_umi_mappy(seqs, min_ctg)
    if contigs is None:
        contigs = assemble_umi(seqs, min_ov, max_mm, min_ctg, seed_k,
                               use_internal_anchor=use_anchor,
                               internal_min_verify=anchor_verify,
                               use_minimizer_dedup=use_dedup,
                               use_cross_attempt_evidence=use_cross_evidence,
                               adaptive_seed_k=adaptive_seed_k,
                               fallback_seed_k=fallback_seed_k,
                               join_map=join_map,
                               min_join_support=min_join_support)

    draft_contigs = contigs
    contigs = _polish_all(contigs, seqs, min_ov, max_mm, seed_k)
    _write_contigs(barcode, contigs, out_file, lock, max_contigs)
    if join_map is not None:
        _write_join_trace(barcode, draft_contigs, contigs, join_map, len(seqs),
                          join_trace_path, lock, max_contigs)


def process_barcode_pe(barcode, shared_meta_data1, shared_meta_data2, lock):
    """PE drop-in for denovo_clfr_ram.process_barcode_pe."""
    min_ctg  = _CFG["min_ctg"]
    min_ov   = _CFG["min_ov"]
    max_mm   = _CFG["max_mm"]
    seed_k   = _CFG["seed_k"]
    out_file = _CFG["out_file"].format(id=_CFG["out_id"])
    use_mp   = _CFG["use_mappy"]
    use_anchor  = _CFG["use_internal_anchor"]
    anchor_verify = _CFG["internal_min_verify"]
    max_contigs = _CFG["max_contigs"]
    use_dedup = _CFG["use_minimizer_dedup"]
    use_cross_evidence = _CFG["use_cross_attempt_evidence"]
    adaptive_seed_k = _CFG["adaptive_seed_k"]
    fallback_seed_k = _CFG["fallback_seed_k"]
    min_join_support = _CFG["min_join_support"]
    join_trace_path = _CFG["join_trace"]

    r1 = _seqs_from_meta(shared_meta_data1, barcode)
    r2 = _seqs_from_meta(shared_meta_data2, barcode)
    seqs = r1 + r2
    if not seqs:
        return

    join_map = {} if join_trace_path else None

    contigs = None
    if use_mp is not False:
        contigs = _assemble_umi_mappy(seqs, min_ctg)
    if contigs is None:
        contigs = assemble_umi(seqs, min_ov, max_mm, min_ctg, seed_k,
                               use_internal_anchor=use_anchor,
                               internal_min_verify=anchor_verify,
                               use_minimizer_dedup=use_dedup,
                               use_cross_attempt_evidence=use_cross_evidence,
                               adaptive_seed_k=adaptive_seed_k,
                               fallback_seed_k=fallback_seed_k,
                               join_map=join_map,
                               min_join_support=min_join_support)

    draft_contigs = contigs
    contigs = _polish_all(contigs, seqs, min_ov, max_mm, seed_k)
    _write_contigs(barcode, contigs, out_file, lock, max_contigs)
    if join_map is not None:
        _write_join_trace(barcode, draft_contigs, contigs, join_map, len(seqs),
                          join_trace_path, lock, max_contigs)


# ── sgrep TSV parsing (same format as denovo_clfr_ram.add_sgrep_line) ──────────

@lru_cache(maxsize=4096)
def _is_low_complexity_barcode(barcode):
    """Reject barcodes containing an ambiguous base call ('N'), which
    should be unreachable for a barcode that survived correction against a
    fixed whitelist and signals a base-calling artifact, not a real UMI
    (denovo.md sec 29 first described these as "AAAA...-style artifacts";
    sec 60 traced a concrete production failure to them).

    NOT a general low-complexity/homopolymer-fraction filter: an earlier
    version rejected any barcode >=80% one base, copying
    denovo_qc_probe.py's is_low_complexity() threshold verbatim. That
    threshold is fine for is_low_complexity()'s own job (nudge a diagnostic
    probe's random sample away from artifacts -- over-excluding there just
    means a slightly different sample, never data loss), but reused as a
    hard rejection gate on production assembly it was a false-positive
    generator: real barcodes in this UMI design are legitimately
    single-base-dominated (e.g. a plain "AAAAAAAAAAAAAAA" with normal
    reads and a normal ~1.3kb contig is common), and thread-through-the-
    fraction-threshold checking on a real 3000-barcode sample (denovo.md
    sec 60) found 339 of them (11.3%) -- all producing ordinary contigs,
    none resembling the actual incident -- would have been silently
    dropped. Worse, fraction alone can't even cleanly separate good from
    bad: one of the actual incident barcodes and the ordinary
    "AAAAAAAAAAAAAAA" case sit at the same >=0.93 single-base fraction, so
    no threshold choice avoids both false positives and false negatives.
    N-in-barcode has no such ambiguity and produced zero false positives
    on the same sample. The two non-N incident barcodes this therefore
    misses are still bounded by _extend_one_contig's max_contig_len --
    that backstop doesn't need to guess the cause in advance.
    """
    return "N" in barcode


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
    if _is_low_complexity_barcode(bc):
        return False
    rid = ">" + info[0][22:]
    seq = info[1]
    meta_data[bc].append(rid)
    meta_data[bc].append(seq)
    return True


def _add_fastq_record(meta_data, header, seq):
    """FASTQ counterpart of _add_sgrep_line, for --r2_format fastq.

    The sgrep TSV's readname column is the FASTQ header with the leading '@'
    replaced by nothing and the quality lines dropped, so the barcode sits at
    the same offsets once the '@' is accounted for. Parsed from the header
    text rather than by offset, though, since a FASTQ fed here may not have
    gone through reformat_fasta2's column layout at all.
    """
    rid = header.rstrip("\n")[1:].split("\t")[0]
    try:
        bc = rid.split("#")[1].split("/")[0]
    except IndexError:
        return False
    if len(bc) != 15 or _is_low_complexity_barcode(bc):
        return False
    meta_data[bc].append(">" + rid)
    meta_data[bc].append(seq.rstrip("\n"))
    return True


def _iter_se_chunks_fastq(r2_path, start_idx, n_line_chunk):
    """Same contract as _iter_se_chunks but reading a BARCODE-SORTED FASTQ.

    Exists so the learned-score path does not have to materialize a second
    copy of the reads: denovo_read_features.py already needs a barcode-sorted
    FASTQ (it needs the quality strings), and without this the same data would
    then be rewritten as a TSV purely to be assembled. Measured tradeoff in
    denovo.md sec 65.

    start_idx/n_line_chunk stay in units of READS, not lines, so --start_idx
    and --n_line_chunk keep meaning the same thing across both formats.
    """
    opener = gzip.open if r2_path.endswith(".gz") else open
    with opener(r2_path, "rt") as f:
        for _ in range(start_idx):
            if not f.readline():
                break
            f.readline(); f.readline(); f.readline()
        chunk_start = start_idx
        while True:
            meta_data2 = defaultdict(list)
            n_reads = 0
            for _ in range(n_line_chunk):
                header = f.readline()
                if not header:
                    break
                seq = f.readline()
                f.readline()
                f.readline()
                if _add_fastq_record(meta_data2, header, seq):
                    n_reads += 1
            if n_reads == 0:
                break
            yield chunk_start, meta_data2
            chunk_start += n_reads


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


_worker_meta_data1 = None
_worker_meta_data2 = None
_worker_lock = None


def _init_pool_worker(meta_data1, meta_data2, lock, cfg):
    """
    Pool(initializer=...) target: runs once per worker process at pool
    startup, not once per barcode. meta_data1/meta_data2/lock are sent
    through the officially-supported process-bootstrap pickling channel
    (the same one Process(args=...) uses) -- this is the one place a
    Lock object is actually allowed to be pickled/shared at all.

    Replaces mp.Manager().dict(): a Manager dict instead proxies every
    single .get(barcode) call through a separate IPC server process --
    one pickled round trip per barcode, real cost at millions-of-barcodes
    scale (same fix already applied to denovo_clfr_ram.py's megahit
    path; this file's own multiprocessing wiring had not been updated to
    match). NOTE: do NOT instead pass the dicts as per-task starmap/map
    arguments -- Pool distributes each task's args through its own
    internal queue and would re-pickle the whole dict on every single
    task, which is worse than Manager, not better.

    cfg must also be forwarded explicitly (not just meta_data/lock):
    under the 'fork' start method (Linux default) a worker inherits the
    parent's already-configure()'d _CFG for free via copy-on-write
    memory, so this looked unnecessary in local testing -- but under
    'spawn' (macOS/Windows default, verified empirically here) each
    worker re-imports this module fresh in a brand new interpreter,
    resetting _CFG to its hardcoded defaults and silently discarding
    every configure() call the parent made (min_ctg, out_file, polish
    settings, everything). Confirmed by reproduction: with configure()
    called in the parent only, spawned workers tried to write
    'denovo/final_contigs_0.fa' (the hardcoded default) instead of the
    configured path.
    """
    global _worker_meta_data1, _worker_meta_data2, _worker_lock
    _worker_meta_data1 = meta_data1
    _worker_meta_data2 = meta_data2
    _worker_lock = lock
    _CFG.update(cfg)


def _pool_process_barcode_pe(barcode):
    process_barcode_pe(barcode, _worker_meta_data1, _worker_meta_data2, _worker_lock)


def _pool_process_barcode_se(barcode):
    process_barcode_se(barcode, _worker_meta_data2, _worker_lock)


def _process_pe_metadata(meta_data1, meta_data2, num_processes):
    """
    chunksize=1: per-barcode cost varies enormously (a few ms for a typical
    UMI vs multiple seconds for a high-depth/pathological one). Pool.map's
    default chunksize hands each worker a large contiguous slice of the
    barcode list up front; if even one slow barcode lands in a worker's
    slice, that whole slice blocks while every other worker sits idle
    waiting for the next chunk boundary -- confirmed in production via
    `ps` showing the active-worker count visibly drop mid-chunk instead of
    staying near num_processes throughout. chunksize=1 makes every worker
    pull one barcode at a time so slow barcodes only ever block themselves.
    """
    if num_processes == 1:
        lock = _NullLock()
        for barcode in meta_data2.keys():
            process_barcode_pe(barcode, meta_data1, meta_data2, lock)
    else:
        import multiprocessing as mp
        lock = mp.Lock()
        with mp.Pool(num_processes, initializer=_init_pool_worker,
                     initargs=(meta_data1, meta_data2, lock, dict(_CFG))) as pool:
            pool.map(_pool_process_barcode_pe, meta_data2.keys(), chunksize=1)
    print("denovo_BC_counts={}".format(len(meta_data2)))
    return sum(len(v) // 2 for v in meta_data2.values())


def _process_se_metadata(meta_data2, num_processes):
    if num_processes == 1:
        lock = _NullLock()
        for barcode in meta_data2.keys():
            process_barcode_se(barcode, meta_data2, lock)
    else:
        import multiprocessing as mp
        lock = mp.Lock()
        with mp.Pool(num_processes, initializer=_init_pool_worker,
                     initargs=(None, meta_data2, lock, dict(_CFG))) as pool:
            pool.map(_pool_process_barcode_se, meta_data2.keys(), chunksize=1)
    print("denovo_BC_counts={}".format(len(meta_data2)))
    return sum(len(v) // 2 for v in meta_data2.values())


# ── standalone pipeline CLI ───────────────────────────────────────────────────

def _build_arg_parser():
    import argparse

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
    ap.add_argument("--seed_k", type=int, default=10,
                    help="k-mer length for overlap pre-filter [10]")
    ap.add_argument("--adaptive_seed_k", action="store_true",
                    help="retry with --fallback_seed_k only if primary k-mer run has no valid contig")
    ap.add_argument("--fallback_seed_k", type=int, default=10,
                    help="fallback k-mer length for adaptive_seed_k [10]")
    ap.add_argument("--nth_of_nodes", type=int, default=0)
    ap.add_argument("--r1", type=str, default="denovo/data_R1_sorted.tsv")
    ap.add_argument("--r2", type=str, default="denovo/data_R2_sorted.tsv")
    ap.add_argument("--r2_format", choices=["tsv", "fastq"], default="tsv",
                    help="format of --r2. 'fastq' reads a BARCODE-SORTED FASTQ "
                         "(.gz ok) directly, skipping the TSV conversion -- only "
                         "useful when a sorted FASTQ already exists for the "
                         "learned-score path (denovo_read_features.py needs the "
                         "quality strings). SE only; PE still requires TSVs.")
    ap.add_argument("--n", type=int, default=None,
                    help="only assemble the first N UMIs total, across all chunks "
                         "(config: frag_de_novo.assembly_N_umi); default/empty = all UMIs")
    ap.add_argument("--no_polish", action="store_true",
                    help="skip post-assembly majority-vote consensus correction (on by default)")
    ap.add_argument("--no_internal_anchor", action="store_true",
                    help="skip internal k-mer anchor fallback extension (on by default)")
    ap.add_argument("--internal_min_verify", type=int, default=60,
                    help="min confirmed overlap length for the internal-anchor fallback [60]")
    ap.add_argument("--max_contigs", type=int, default=8,
                    help="maximum final contigs emitted per UMI after merge/dedupe [8]")
    ap.add_argument("--minimizer_dedup", action="store_true",
                    help="collapse near-identical reads (same canonical minimizer) before "
                         "assembly -- OFF by default: confirmed to collapse genuinely-distinct "
                         "reads sharing a conserved motif on real 16S data, see "
                         "_canonical_minimizer's docstring. Opt-in only, for experimentation")
    ap.add_argument("--no_cross_attempt_evidence", action="store_true",
                    help="restrict collective rescue to the current attempt's own "
                         "leftover reads only, disabling cross-attempt evidence "
                         "sharing (on by default; see assemble_umi docstring)")
    ap.add_argument("--mappy", action="store_true",
                    help="try mappy/minimap2 (ava-sr) overlap detection before "
                         "assemble_umi -- OFF by default: a real 1.5M production run "
                         "was confirmed to silently take this path whenever `mappy` "
                         "happened to be importable, bypassing every correctness fix "
                         "this project has made (no mismatch-rate check, no chimeric-"
                         "merge safety, always one contig per barcode, no single-read "
                         "UMI support). Opt-in only, for experimentation; see "
                         "_assemble_umi_mappy's docstring")
    ap.add_argument("--min-join-support", type=int, default=1,
                    help="verified distinct raw reads required across a newly created "
                         "interior join [1]. Values above 1 are an experimental "
                         "chimera guard and apply only to the seed-extension path.")
    ap.add_argument("--join-trace", dest="join_trace", type=str, default=None,
                    help="write one TSV row per join that survives into an emitted "
                         "contig (barcode, contig_idx, join_type, ov_len, "
                         "junction_pos, contig_len, bridge_reads, n_evidence_reads, "
                         "ctx_pos_agrees) to PATH -- OFF by default (None): explicit "
                         "opt-IN only, matching this file's own convention (see "
                         "_configure_from_args's docstring for why an opt-OUT default "
                         "has twice caused a silent-default production incident here). "
                         "Pure observation -- computes verified bridge-read support at "
                         "each join for a future spanning guard but never changes an "
                         "assembly decision itself; leaving this unset costs nothing")
    return ap


def _configure_from_args(args):
    """
    Maps parsed CLI args to configure() kwargs -- pulled out of _main_cli so
    tests can verify the actual end-to-end default (e.g. "no flags passed")
    without running the whole read-streaming pipeline. This mapping is the
    exact spot two real bugs lived in before:
    - --no_minimizer_dedup (an opt-OUT flag) meant "flag absent" silently
      forced use_minimizer_dedup=True via `not args.no_minimizer_dedup`,
      re-enabling a feature confirmed to corrupt real 16S assemblies, even
      though assemble_umi/configure/_CFG's own defaults all say False.
      Fixed by making it an explicit opt-IN flag (--minimizer_dedup,
      default False) instead.
    - use_mappy had no CLI flag at all and defaulted to None ("auto-detect
      and prefer if importable"). A real 1.5M production run silently took
      the mappy path the whole time because the server happened to have it
      installed, making every assemble_umi fix validated on this project a
      no-op in production. Fixed the same way: explicit opt-IN flag
      (--mappy, default False).
    """
    configure(min_ctg_len=args.min_ctg_len, min_overlap=args.min_overlap,
              max_mismatch=args.max_mismatch, out_id=args.nth_of_nodes,
              seed_k=args.seed_k, adaptive_seed_k=args.adaptive_seed_k,
              fallback_seed_k=args.fallback_seed_k,
              polish=not args.no_polish,
              use_internal_anchor=not args.no_internal_anchor,
              internal_min_verify=args.internal_min_verify,
              max_contigs=args.max_contigs,
              use_minimizer_dedup=args.minimizer_dedup,
              use_cross_attempt_evidence=not args.no_cross_attempt_evidence,
              use_mappy=args.mappy,
              min_join_support=args.min_join_support,
              join_trace=args.join_trace)


def _git_version_string():
    """
    Identify exactly which commit/working-tree state of this script produced
    a given run's output. Motivated by a real production incident: a full
    1.5M-UMI run reported yield stats identical to a known pre-fix baseline,
    and resolving whether it actually used the latest code took several
    rounds of manually checking `git log`/timestamps after the fact. Runs
    now self-report this instead of relying on out-of-band reconstruction.
    Never raises -- a missing git binary or a non-repo checkout must not
    block an actual production run over a diagnostics nicety.
    """
    import subprocess
    script_dir = os.path.dirname(os.path.abspath(__file__))
    try:
        commit = subprocess.check_output(
            ["git", "-C", script_dir, "rev-parse", "HEAD"],
            stderr=subprocess.DEVNULL).decode().strip()
        branch = subprocess.check_output(
            ["git", "-C", script_dir, "rev-parse", "--abbrev-ref", "HEAD"],
            stderr=subprocess.DEVNULL).decode().strip()
        dirty = subprocess.call(
            ["git", "-C", script_dir, "diff", "--quiet", "HEAD", "--",
             os.path.abspath(__file__)],
            stderr=subprocess.DEVNULL) != 0
        return "commit={} branch={} local_edits_to_this_file={}".format(commit, branch, dirty)
    except Exception as exc:
        return "unknown (git unavailable or not a checkout: {})".format(exc)


def _main_cli():
    import datetime
    import subprocess

    args = _build_arg_parser().parse_args()
    _configure_from_args(args)

    if not os.path.isdir("denovo"):
        os.makedirs("denovo")

    if args.join_trace:
        trace_dir = os.path.dirname(os.path.abspath(args.join_trace))
        if trace_dir and not os.path.isdir(trace_dir):
            os.makedirs(trace_dir)
        print("join_trace={}".format(args.join_trace), flush=True)

    version_line = "run_start={} denovo_seed_olc.py {}".format(
        datetime.datetime.now(), _git_version_string())
    print(version_line, flush=True)
    with open(os.path.join("denovo", "version.txt"), "a") as vf:
        vf.write(version_line + "\n")

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
            se_reader = (_iter_se_chunks_fastq if args.r2_format == "fastq"
                         else _iter_se_chunks)
            for chunk_start, m2 in se_reader(args.r2, args.start_idx, args.n_line_chunk):
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

    # polish: corrupt a correctly-assembled contig at one position, verify
    # majority vote from the (uncorrupted) reads flips it back
    reads_dense = _make_reads(FRAG, 150, 15)  # dense overlap -> high per-position coverage
    ctg_dense = assemble_umi(reads_dense, min_ctg=400)
    polish_ok = False
    if ctg_dense and ctg_dense[0] in FRAG:
        good_contig = ctg_dense[0]
        err_pos = len(good_contig) // 2
        wrong_base = "ACGT"[("ACGT".index(good_contig[err_pos]) + 1) % 4]
        corrupted = good_contig[:err_pos] + wrong_base + good_contig[err_pos + 1:]
        polished = polish_contig(corrupted, reads_dense, min_ov=20, max_mm=0.05, seed_k=10,
                                  min_coverage=3, vote_concordance=0.6, kmer_step=5)
        polish_ok = (polished == good_contig) or (polished[err_pos] == good_contig[err_pos])
    print("[polish]   corrected={}".format(polish_ok))

    # internal-anchor fallback: simulate a real-data pattern found in 16S
    # rRNA barcodes (noisy/chimeric read edge fused onto an otherwise
    # accurate long middle section) -- a boundary-only suffix/prefix
    # overlap can never see the genuine match since it's not at the
    # read's literal start; the internal-anchor fallback should find it.
    # NOTE: assemble_umi always picks the LONGEST read as the initial
    # seed/contig, so the clean read must be the longer of the two here --
    # the current fallback only searches for an anchor inside the
    # candidate against the (assumed-reliable) contig boundary, not the
    # reverse.
    frag3 = _rand_seq(500)
    seed3 = frag3[:300]                          # longer -> becomes the seed/contig
    junk_prefix = _rand_seq(40)                   # unrelated garbage stitched onto a real overlap
    bridging_read = junk_prefix + frag3[200:400]  # shorter -> tested as a candidate; real match starts at position 40, not 0
    ctg_no_anchor = assemble_umi([seed3, bridging_read], min_ov=30, min_ctg=350,
                                  use_internal_anchor=False)
    ctg_with_anchor = assemble_umi([seed3, bridging_read], min_ov=30, min_ctg=350,
                                    use_internal_anchor=True, internal_min_verify=60)
    anchor_ok = bool(
        (not ctg_no_anchor or len(ctg_no_anchor[0]) < 400)
        and ctg_with_anchor and len(ctg_with_anchor[0]) >= 400
        and ctg_with_anchor[0] in frag3
    )
    print("[internal_anchor] no_anchor_len={}  with_anchor_len={}  rescued={}".format(
        len(ctg_no_anchor[0]) if ctg_no_anchor else 0,
        len(ctg_with_anchor[0]) if ctg_with_anchor else 0,
        anchor_ok))

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
    if not polish_ok:
        print("FAIL: polish did not correct the injected substitution error", file=sys.stderr)
        ok = False
    if not anchor_ok:
        print("FAIL: internal-anchor fallback did not rescue the noisy-edge read", file=sys.stderr)
        ok = False
    if ok:
        print("correctness: OK")


if __name__ == "__main__":
    import sys
    if len(sys.argv) == 1:
        _run_selftest()
    else:
        _main_cli()
