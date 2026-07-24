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
Run directly from a Snakemake work dir containing denovo/data_R1_sgrep.tsv
and denovo/data_R2_sgrep.tsv (same layout denovo_clfr_ram.py expects).
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
        --r2 denovo/data_R2_sgrep.tsv \\
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
Known gap: only searches for the anchor inside the candidate, not
inside the contig itself (see _internal_anchor_extend_3prime
docstring). On by default; disable with configure(use_internal_anchor=False)
or CLI --no_internal_anchor.
"""

import itertools
import os
from collections import defaultdict, Counter

# ── module-level config (set via configure() before multiprocessing) ──────────

_CFG = {
    "min_ctg":   400,
    "min_ov":    20,
    "max_mm":    0.05,
    "out_id":    0,
    "out_file":  "denovo/final_contigs_{id}.fa",
    "seed_k":    10,    # k-mer size for overlap pre-filter
    "use_mappy": None,  # None = auto-detect, True/False = force
    "polish":              True,  # majority-vote consensus correction after assembly
    "polish_min_coverage": 3,     # min total votes (incl. 1 implicit vote for original base) to consider flipping
    "polish_vote_concordance": 0.6,  # winning base must hold >= this fraction of votes to flip
    "polish_kmer_step":    5,     # stride for sampling k-mers when re-aligning reads to the contig
    "use_internal_anchor":  True,  # fallback: k-mer anchor anywhere in a read, not just its literal ends
    "internal_min_verify":  60,    # min confirmed overlap length for the internal-anchor fallback
}


def configure(min_ctg_len=400, min_overlap=20, max_mismatch=0.05,
              out_id=0, out_file="denovo/final_contigs_{id}.fa", use_mappy=None,
              polish=True, polish_min_coverage=3, polish_vote_concordance=0.6,
              polish_kmer_step=5, use_internal_anchor=True, internal_min_verify=60):
    """Call once in the parent process before spawning Pool workers."""
    _CFG["min_ctg"]   = min_ctg_len
    _CFG["min_ov"]    = min_overlap
    _CFG["max_mm"]    = max_mismatch
    _CFG["out_id"]    = out_id
    _CFG["out_file"]  = out_file
    _CFG["use_mappy"] = use_mappy
    _CFG["polish"]                   = polish
    _CFG["polish_min_coverage"]      = polish_min_coverage
    _CFG["polish_vote_concordance"]  = polish_vote_concordance
    _CFG["polish_kmer_step"]         = polish_kmer_step
    _CFG["use_internal_anchor"]      = use_internal_anchor
    _CFG["internal_min_verify"]      = internal_min_verify


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


def _internal_anchor_extend_3prime_indexed(contig, pool, kmer_index, unused_set,
                                            min_ov, max_mm, seed_k, min_verify):
    """
    Index-accelerated 3' (right-end) internal-anchor extension: looks up
    contig's tail k-mers directly in the pre-built pool-wide index
    instead of rescanning every remaining candidate's full length.
    Only considers candidates whose pool index is still in unused_set.

    Returns (new_contig, used_pool_index) or (None, None).
    """
    n = len(contig)
    check_len = min(n, max(min_ov * 3, seed_k * 4))
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
            contig_region = contig[contig_start:contig_start + overlap_len]
            cand_region = candidate[j:j + overlap_len]
            mm = sum(x != y for x, y in zip(contig_region, cand_region))
            if mm / overlap_len <= max_mm:
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
                                            min_ov, max_mm, seed_k, min_verify):
    """
    Index-accelerated 5' (left-end) internal-anchor extension: symmetric
    counterpart of the 3' version, looking up contig's head k-mers and
    checking whether a candidate's content BEFORE the anchor can be
    prepended.

    Returns (new_contig, used_pool_index) or (None, None).
    """
    n = len(contig)
    check_len = min(n, max(min_ov * 3, seed_k * 4))
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
            if cand_prefix_start < 0:
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
            mm = sum(x != y for x, y in zip(contig_region, cand_region))
            if mm / overlap_len <= max_mm:
                # see 3' version: prefer smallest pool index (matches
                # pre-index behavior), not the single biggest gain
                sort_key = (idx, -overlap_len)
                if best is None or sort_key < best[0]:
                    best = (sort_key, candidate[:cand_prefix_start] + contig, idx)

    if best is not None:
        return best[1], best[2]
    return None, None


def internal_anchor_extend_indexed(contig, pool, kmer_index, unused_set,
                                    min_ov, max_mm, seed_k=10, min_verify=60):
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

    KNOWN LIMITATION: only searches for the anchor inside candidate
    reads, assuming the contig's own boundary is reliable -- doesn't
    handle the reverse (contig's edge is the noisy one). In practice
    rare, since a contig's boundary was either the original longest raw
    read or the product of a prior verified merge. The gap is only the
    very first seed itself having a noisy edge.

    Returns (new_contig, used_pool_index) or (None, None).
    """
    result = _internal_anchor_extend_3prime_indexed(
        contig, pool, kmer_index, unused_set, min_ov, max_mm, seed_k, min_verify)
    if result[0] is not None:
        return result
    return _internal_anchor_extend_5prime_indexed(
        contig, pool, kmer_index, unused_set, min_ov, max_mm, seed_k, min_verify)


def _extend_one_contig(pool, available, min_ov, max_mm, seed_k, use_internal_anchor=True,
                       internal_min_verify=60, kmer_index_holder=None):
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
    """
    seed_idx = min(available)
    contig = pool[seed_idx]
    used = {seed_idx}
    unused = set(available) - {seed_idx}

    changed = True
    while changed and unused:
        changed = False
        # iterate in ascending index order (== pool's longest-first sort
        # order): sets don't preserve insertion order, and trying
        # candidates in a different order than before changes which
        # merge happens first in this greedy algorithm -- which can change
        # the final assembled contig even when every individual merge is
        # independently valid. sorted() restores the original, deterministic
        # "prefer the longest remaining read" trial order.
        for i in sorted(unused):
            seq = pool[i]
            extended = False

            for cand in (seq, rc(seq)):
                # try extending contig at 3' end
                ov = suffix_prefix_overlap(contig, cand, min_ov, max_mm, seed_k)
                if ov:
                    contig += cand[ov:]
                    unused.remove(i)
                    used.add(i)
                    changed = True
                    extended = True
                    break

                # try extending contig at 5' end
                ov = suffix_prefix_overlap(cand, contig, min_ov, max_mm, seed_k)
                if ov:
                    contig = cand + contig[ov:]
                    unused.remove(i)
                    used.add(i)
                    changed = True
                    extended = True
                    break

            if extended:
                # restart scan so new contig ends are retried against all unused
                break

        if changed or not unused or not use_internal_anchor:
            continue

        # boundary extension is fully exhausted -- try the internal-anchor
        # fallback once before giving up on this contig. build the shared
        # index lazily, only now that it's actually needed.
        if kmer_index_holder[0] is None:
            kmer_index_holder[0] = _build_pool_kmer_index(pool, seed_k)
        new_contig, used_idx = internal_anchor_extend_indexed(
            contig, pool, kmer_index_holder[0], unused, min_ov, max_mm, seed_k,
            min_verify=internal_min_verify)
        if new_contig is not None:
            contig = new_contig
            unused.remove(used_idx)
            used.add(used_idx)
            changed = True

    return contig, used


def assemble_umi(seqs, min_ov=20, max_mm=0.05, min_ctg=400, seed_k=10, max_contigs=4,
                 use_internal_anchor=True, internal_min_verify=60):
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
    max_contigs : stop after this many accepted contigs [4] (matches
                  _write_contigs' existing per-barcode cap)
    use_internal_anchor : also try internal_anchor_extend_indexed() as a fallback
                  when boundary suffix/prefix extension stalls [True]
    internal_min_verify : min confirmed overlap length for the internal
                  anchor fallback to accept a match [60]

    Returns list of contig sequences (0 or more per UMI).
    """
    if not seqs:
        return []

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

    available = set(range(len(pool)))
    contigs = []
    while available and len(contigs) < max_contigs:
        contig, used = _extend_one_contig(pool, available, min_ov, max_mm, seed_k,
                                          use_internal_anchor, internal_min_verify,
                                          kmer_index_holder=kmer_index_holder)
        if len(contig) >= min_ctg:
            contigs.append(contig)
        # always drop every read the attempt consumed (even just the seed
        # itself, on a failed/orphan attempt) so available strictly shrinks
        # and a genuinely separate fragment among the rest still gets a shot
        available -= used

    return _dedupe_and_merge_contigs(contigs, min_ov, max_mm, seed_k)


def _dedupe_and_merge_contigs(contigs, min_ov, max_mm, seed_k):
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
    2. Pure internal containment (one contig sits entirely inside
       another, not at either edge -- suffix_prefix_overlap only checks
       boundaries so it won't catch this) -- the contained one adds no
       new sequence, so it's simply dropped.
    """
    if len(contigs) <= 1:
        return contigs

    # 1. merge any pair with a real boundary overlap. use_internal_anchor=False:
    # already-assembled contigs have clean, verified boundaries (either from
    # normal extension or a verified anchor merge), so the noisy-read-edge
    # rescue this fallback exists for doesn't apply here.
    pool = sorted(set(contigs), key=lambda s: (-len(s), s))
    available = set(range(len(pool)))
    merged = []
    while available:
        contig, used = _extend_one_contig(pool, available, min_ov, max_mm, seed_k,
                                          use_internal_anchor=False)
        merged.append(contig)
        available -= used

    # 2. drop pure containment (substring anywhere, not just at a boundary)
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
            mismatches = sum(x != y for x, y in zip(contig_region, read_region))
            if mismatches / overlap_len > max_mm:
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

    seqs = _seqs_from_meta(shared_meta_data2, barcode)
    if not seqs:
        return

    contigs = None
    if use_mp is not False:
        contigs = _assemble_umi_mappy(seqs, min_ctg)
    if contigs is None:
        contigs = assemble_umi(seqs, min_ov, max_mm, min_ctg, seed_k,
                               use_internal_anchor=use_anchor,
                               internal_min_verify=anchor_verify)

    contigs = _polish_all(contigs, seqs, min_ov, max_mm, seed_k)
    _write_contigs(barcode, contigs, out_file, lock)


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

    r1 = _seqs_from_meta(shared_meta_data1, barcode)
    r2 = _seqs_from_meta(shared_meta_data2, barcode)
    seqs = r1 + r2
    if not seqs:
        return

    contigs = None
    if use_mp is not False:
        contigs = _assemble_umi_mappy(seqs, min_ctg)
    if contigs is None:
        contigs = assemble_umi(seqs, min_ov, max_mm, min_ctg, seed_k,
                               use_internal_anchor=use_anchor,
                               internal_min_verify=anchor_verify)

    contigs = _polish_all(contigs, seqs, min_ov, max_mm, seed_k)
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
    ap.add_argument("--no_polish", action="store_true",
                    help="skip post-assembly majority-vote consensus correction (on by default)")
    ap.add_argument("--no_internal_anchor", action="store_true",
                    help="skip internal k-mer anchor fallback extension (on by default)")
    ap.add_argument("--internal_min_verify", type=int, default=60,
                    help="min confirmed overlap length for the internal-anchor fallback [60]")
    args = ap.parse_args()

    configure(min_ctg_len=args.min_ctg_len, min_overlap=args.min_overlap,
              max_mismatch=args.max_mismatch, out_id=args.nth_of_nodes,
              polish=not args.no_polish,
              use_internal_anchor=not args.no_internal_anchor,
              internal_min_verify=args.internal_min_verify)

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