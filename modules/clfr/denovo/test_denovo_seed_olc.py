"""
Unit tests for denovo_seed_olc.py.

Run:
    python3 -m unittest test_denovo_seed_olc -v
or, from this directory:
    python3 test_denovo_seed_olc.py

Pure stdlib (unittest + random + multiprocessing), matching the module's
own zero-dependency design.
"""
import random
import sys
import os
import tempfile
import unittest
from unittest import mock

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import denovo_seed_olc as m


def rand_seq(rng, n):
    return "".join(rng.choice("ACGT") for _ in range(n))


class TestSuffixPrefixOverlap(unittest.TestCase):
    def test_exact_overlap_found(self):
        a = "AAAA" + "ACGTACGTAC"
        b = "ACGTACGTAC" + "TTTT"
        self.assertEqual(m.suffix_prefix_overlap(a, b, min_ov=10, max_mm=0.05), 10)

    def test_no_overlap_below_min_ov(self):
        a = "AAAAAAAAAA"
        b = "TTTTTTTTTT"
        self.assertEqual(m.suffix_prefix_overlap(a, b, min_ov=5, max_mm=0.05), 0)

    def test_tolerates_mismatch_within_budget(self):
        rng = random.Random(1)
        core = rand_seq(rng, 40)
        a = rand_seq(rng, 20) + core
        mutated = list(core)
        mutated[0] = "A" if core[0] != "A" else "C"  # exactly 1 mismatch / 40 = 2.5% <= 5%
        b = "".join(mutated) + rand_seq(rng, 20)
        self.assertEqual(m.suffix_prefix_overlap(a, b, min_ov=20, max_mm=0.05), 40)

    def test_rejects_mismatch_over_budget(self):
        rng = random.Random(2)
        core = rand_seq(rng, 20)
        a = rand_seq(rng, 20) + core
        mutated = list(core)
        for i in range(4):  # 4/20 = 20% > 5%
            mutated[i] = "A" if core[i] != "A" else "C"
        b = "".join(mutated) + rand_seq(rng, 20)
        self.assertEqual(m.suffix_prefix_overlap(a, b, min_ov=20, max_mm=0.05), 0)

    def test_boundary_index_never_hides_a_valid_overlap(self):
        pool = [
            "TTTTTTTTTTACGTACGTACGTACGTAC",
            "ACGTACGTACGTACGTACGGGGGGGGGG",
            "CCCCCCCCCCGTACGTACGTACGTACGT",
            "GATCGATCGATCGATCGATCGATCGATC",
        ]
        index = m._build_boundary_kmer_index(pool, min_ov=20, seed_k=10)
        available = set(range(len(pool)))

        for contig in pool:
            indexed = m._boundary_overlap_candidates(
                contig, available, index, min_ov=20, seed_k=10)
            for idx, seq in enumerate(pool):
                valid = any(
                    m.suffix_prefix_overlap(contig, cand, 20, 0.05, 10)
                    or m.suffix_prefix_overlap(cand, contig, 20, 0.05, 10)
                    for cand in (seq, m.rc(seq))
                )
                if valid:
                    self.assertIn(idx, indexed)


class TestWithinMismatchBudget(unittest.TestCase):
    """
    _within_mismatch_budget replaced a plain sum(x!=y for x,y in zip(...))
    with an early-exit version for speed (see denovo_seed_olc.py commit
    history: suffix_prefix_overlap was 83-85% of total assemble_umi()
    runtime on real 16S data). This must return exactly the same verdict
    as the naive full count for every case, or the speedup would be a
    silent correctness regression.
    """
    def _naive(self, a, b, ov, max_mm):
        mm = sum(x != y for x, y in zip(a, b))
        return mm / ov <= max_mm

    def test_matches_naive_full_scan_on_random_cases(self):
        rng = random.Random(3)
        for _ in range(500):
            n = rng.randint(1, 200)
            a = rand_seq(rng, n)
            # b = a with a random number of substitutions
            b = list(a)
            n_mut = rng.randint(0, n)
            for i in rng.sample(range(n), n_mut):
                b[i] = rng.choice([c for c in "ACGT" if c != a[i]])
            b = "".join(b)
            max_mm = rng.choice([0.0, 0.01, 0.05, 0.1, 0.5])
            expected = self._naive(a, b, n, max_mm)
            actual = m._within_mismatch_budget(a, b, n, max_mm)
            self.assertEqual(actual, expected, (a, b, n, max_mm))


class TestAssembleUmi(unittest.TestCase):
    def test_all_raw_components_reach_post_assembly_merge(self):
        rng = random.Random(47)
        reads = [rand_seq(rng, 100) for _ in range(8)]
        # The deterministic random reads have no 20 bp suffix/prefix overlap,
        # so each must reach post-assembly merge.  The final-output cap belongs
        # in _write_contigs(), not this raw-component phase.
        self.assertEqual(len(m.assemble_umi(reads, min_ctg=100)), 8)

    def test_writer_keeps_all_eight_configured_components(self):
        tmp = tempfile.NamedTemporaryFile(delete=False)
        tmp.close()
        try:
            m._write_contigs("BC", ["A" * 100] * 8, tmp.name, m._NullLock(), 8)
            with open(tmp.name) as handle:
                self.assertEqual(sum(line.startswith(">") for line in handle), 8)
        finally:
            os.unlink(tmp.name)

    def test_hi_depth_full_reconstruction(self):
        rng = random.Random(4)
        truth = rand_seq(rng, 600)
        reads = [truth[i:i + 100] for i in range(0, 501, 20)]
        contigs = m.assemble_umi(reads, min_ov=20, max_mm=0.05, min_ctg=400)
        self.assertEqual(len(contigs), 1)
        self.assertEqual(len(contigs[0]), 600)

    def test_single_read_seed_only(self):
        rng = random.Random(5)
        read = rand_seq(rng, 250)
        contigs = m.assemble_umi([read], min_ctg=200)
        self.assertEqual(contigs, [read])

    def test_below_min_ctg_is_discarded(self):
        rng = random.Random(6)
        read = rand_seq(rng, 100)
        contigs = m.assemble_umi([read], min_ctg=400)
        self.assertEqual(contigs, [])


    def test_reverse_complement_read_merges(self):
        rng = random.Random(7)
        truth = rand_seq(rng, 400)
        reads = [truth[:220], m.rc(truth[180:])]
        contigs = m.assemble_umi(reads, min_ov=20, max_mm=0.05, min_ctg=300)
        self.assertEqual(len(contigs), 1)
        self.assertEqual(len(contigs[0]), 400)

    def test_empty_input(self):
        self.assertEqual(m.assemble_umi([]), [])

    def test_min_ctg_applies_after_merge_not_before(self):
        """
        Regression test for a real bug: assemble_umi used to filter each
        raw contig-building attempt by min_ctg BEFORE ever calling
        _dedupe_and_merge_contigs, so two individually-short attempts
        that would genuinely merge into something over the floor were
        each discarded first and never got the chance. On a real
        1.5kb-library production run this was a meaningful contributor
        to low yield: many barcodes' raw attempts individually missed
        the --min_ctg floor while their merged result would have cleared
        it easily. Fixed by collecting every raw attempt regardless of
        length, merging first, and filtering by min_ctg only on the
        final, post-merge contigs.
        """
        rng = random.Random(41)
        shared = rand_seq(rng, 100)
        longer = rand_seq(rng, 150) + shared                      # 250bp, shared = its own tail
        shorter = rand_seq(rng, 30) + shared + rand_seq(rng, 80)   # 210bp, shared in its middle
        # both individual raw attempts (250, 210) are below min_ctg;
        # only their merge (330bp) clears it.
        contigs = m.assemble_umi([longer, shorter], min_ov=20, max_mm=0.05,
                                 min_ctg=260, seed_k=10)
        self.assertEqual(len(contigs), 1)
        self.assertGreaterEqual(len(contigs[0]), 260)


class TestCollectiveAnchorRescue(unittest.TestCase):
    def _add_spaced_errors(self, seq, first, step):
        mutated = list(seq)
        for pos in range(first, len(mutated), step):
            mutated[pos] = "A" if mutated[pos] != "A" else "C"
        return "".join(mutated)

    def test_redundant_noisy_reads_extend_without_relaxing_pairwise_mismatch(self):
        """
        The seed and each tiling read have staggered independent errors.  Their
        real suffix/prefix overlaps are all just above the 5% pairwise mismatch
        ceiling, so ordinary OLC correctly refuses them.  Multiple exact 17-mer
        anchors still place three reads coherently; their shared extension is
        then supported by a pileup rather than a relaxed mismatch threshold.
        """
        rng = random.Random(51)
        truth = rand_seq(rng, 500)
        reads = [
            self._add_spaced_errors(truth[:300], 5, 35),
            self._add_spaced_errors(truth[150:350], 20, 35),
            self._add_spaced_errors(truth[200:400], 20, 35),
            self._add_spaced_errors(truth[250:450], 20, 35),
            self._add_spaced_errors(truth[300:500], 20, 35),
        ]

        baseline = m.assemble_umi(reads, min_ov=20, max_mm=0.05,
                                  min_ctg=100, use_collective_rescue=False)
        rescued = m.assemble_umi(reads, min_ov=20, max_mm=0.05,
                                 min_ctg=100, use_collective_rescue=True)
        self.assertLessEqual(max(map(len, baseline)), 350)
        self.assertGreaterEqual(max(map(len, rescued)), 450)
        self.assertGreater(max(map(len, rescued)), max(map(len, baseline)))

    def test_consumed_reads_can_vote_without_being_consumed_twice(self):
        rng = random.Random(52)
        truth = rand_seq(rng, 400)
        contig = truth[:300]
        pool = [truth[150:400], truth[150:400]]
        index = m._build_pool_kmer_index(pool, m._COLLECTIVE_ANCHOR_K)

        rescued, used = m._collective_anchor_extend(
            contig, pool, index, candidate_set=set(),
            evidence_set={0, 1})

        self.assertEqual(rescued, truth)
        self.assertEqual(used, set())

    def test_evidence_only_extension_cannot_repeat_without_progress(self):
        rng = random.Random(54)
        pool = sorted([rand_seq(rng, 100), rand_seq(rng, 100)],
                      key=lambda seq: (-len(seq), seq))

        def evidence_only(contig, pool_, index, candidates, evidence_set=None):
            if evidence_set is None:
                return None, set()
            return contig + "A", set()

        with mock.patch.object(m, "_collective_anchor_extend",
                               side_effect=evidence_only):
            contig, used = m._extend_one_contig(
                pool, set(range(len(pool))), min_ov=20, max_mm=0.05,
                seed_k=10, use_internal_anchor=False,
                collective_index_holder=[{}],
                collective_evidence_set=set(range(len(pool))),
                boundary_index=m._build_boundary_kmer_index(pool, 20, 10))

        self.assertEqual(contig, pool[0])
        self.assertEqual(used, {0})

    def test_collective_rescue_does_not_rewrite_draft_interior(self):
        rng = random.Random(53)
        truth = rand_seq(rng, 400)
        draft = list(truth[:300])
        draft[100] = "A" if draft[100] != "A" else "C"
        draft = "".join(draft)
        pool = [truth[50:400], truth[50:400]]
        index = m._build_pool_kmer_index(pool, m._COLLECTIVE_ANCHOR_K)

        rescued, _ = m._collective_anchor_extend(
            draft, pool, index, candidate_set={0, 1})

        self.assertEqual(rescued[:len(draft)], draft)
        self.assertEqual(len(rescued), len(truth))


class TestCrossAttemptEvidenceToggle(unittest.TestCase):
    """
    use_cross_attempt_evidence controls whether collective rescue may vote
    using reads already consumed by earlier raw-building attempts for the
    same UMI (True, the wider-evidence behaviour), or only the current
    attempt's own leftover reads (False, the original candidate-only
    semantics). A full 20k real-data run with this True showed both
    lengthened AND shortened main contigs vs the old baseline (lfr.md
    section 21) -- this toggle exists so a dedicated 20k A/B can isolate
    that effect instead of conflating it with unrelated speed fixes.
    """

    def test_defaults_to_on_everywhere(self):
        import inspect
        self.assertTrue(
            inspect.signature(m.assemble_umi).parameters["use_cross_attempt_evidence"].default)
        self.assertTrue(
            inspect.signature(m.configure).parameters["use_cross_attempt_evidence"].default)
        self.assertTrue(m._CFG["use_cross_attempt_evidence"])

    def test_false_restricts_evidence_set_to_none_regardless_of_pool_size(self):
        rng = random.Random(60)
        # small pool (<= _COLLECTIVE_MAX_CROSS_ATTEMPT_READS) would normally
        # get the full-pool evidence set when the toggle is on.
        seqs = [rand_seq(rng, 300) for _ in range(3)]
        captured = []

        def record_and_stop(pool, available, *args, **kwargs):
            captured.append(kwargs.get("collective_evidence_set"))
            seed_idx = min(available)
            return pool[seed_idx], set(available)

        with mock.patch.object(m, "_extend_one_contig", side_effect=record_and_stop):
            m.assemble_umi(seqs, min_ctg=1, use_cross_attempt_evidence=True)
            m.assemble_umi(seqs, min_ctg=1, use_cross_attempt_evidence=False)

        self.assertIsNotNone(captured[0], "on: small pool should get a full-pool evidence set")
        self.assertIsNone(captured[1], "off: must restrict to candidate-only regardless of pool size")

    def test_cli_flag_wired_to_configure(self):
        m.configure(use_cross_attempt_evidence=False)
        self.assertFalse(m._CFG["use_cross_attempt_evidence"])
        m.configure(use_cross_attempt_evidence=True)
        self.assertTrue(m._CFG["use_cross_attempt_evidence"])


class TestCliArgDefaultsEndToEnd(unittest.TestCase):
    """
    Regression test for a real bug: the CLI used to wire
    use_minimizer_dedup=not args.no_minimizer_dedup. --no_minimizer_dedup is
    an opt-OUT flag (action="store_true", default False when absent), so
    running the CLI with NO dedup-related flag at all -- e.g. the exact
    command line used for the real 1.5M production run -- silently computed
    use_minimizer_dedup=True and re-enabled a feature already confirmed (see
    TestMinimizerDedup) to corrupt real 16S assemblies, even though
    assemble_umi/configure/_CFG's own signatures all default to False.
    inspect.signature checks alone (TestMinimizerDedup.test_defaults_to_off_everywhere)
    could not catch this because they never exercise argparse or the
    args-to-configure() mapping -- only an end-to-end parse+map test can.
    Fixed by making minimizer dedup an explicit opt-in flag instead
    (--minimizer_dedup, default False, matching every other default).
    """

    def _configure_kwargs(self, argv):
        captured = {}
        with mock.patch.object(m, "configure", side_effect=lambda **kw: captured.update(kw)):
            args = m._build_arg_parser().parse_args(argv)
            m._configure_from_args(args)
        return captured

    def test_no_flags_matches_every_other_default(self):
        kwargs = self._configure_kwargs(["--sequence_type", "se"])
        self.assertFalse(kwargs["use_minimizer_dedup"],
                          "running with no dedup flag must NOT silently re-enable "
                          "the confirmed-broken minimizer dedup")
        self.assertTrue(kwargs["use_cross_attempt_evidence"])

    def test_explicit_opt_in_flags_flip_correctly(self):
        kwargs = self._configure_kwargs(
            ["--sequence_type", "se", "--minimizer_dedup", "--no_cross_attempt_evidence"])
        self.assertTrue(kwargs["use_minimizer_dedup"])
        self.assertFalse(kwargs["use_cross_attempt_evidence"])


class TestPolishContig(unittest.TestCase):
    def test_corrects_single_read_substitution_error(self):
        rng = random.Random(8)
        truth = rand_seq(rng, 200)
        contig = list(truth)
        contig[100] = "A" if truth[100] != "A" else "C"  # inject 1 error
        contig = "".join(contig)
        # 5 clean reads voting for the true base at that position
        reads = [truth[i:i + 100] for i in range(0, 100, 20)]
        polished = m.polish_contig(contig, reads, min_coverage=3, vote_concordance=0.6)
        self.assertEqual(polished[100], truth[100])


class TestInternalAnchorForward(unittest.TestCase):
    """
    Candidate has a noisy prefix fused onto an otherwise-clean overlap
    with the contig -- suffix_prefix_overlap can't see it (only checks
    literal edges), so the forward internal-anchor fallback exists to
    find the anchor inside the candidate instead.
    """
    def test_rescues_candidate_with_noisy_prefix(self):
        rng = random.Random(9)
        pool_kmer_index = {}
        contig = rand_seq(rng, 400)
        clean_overlap = contig[-80:]
        junk_prefix = rand_seq(rng, 30)
        new_tail = rand_seq(rng, 100)
        candidate = junk_prefix + clean_overlap + new_tail
        kmer_index = m._build_pool_kmer_index([candidate], seed_k=10)
        new_contig, idx = m.internal_anchor_extend_indexed(
            contig, [candidate], kmer_index, {0},
            min_ov=20, max_mm=0.05, seed_k=10, min_verify=60)
        self.assertIsNotNone(new_contig)
        self.assertEqual(new_contig, contig + new_tail)


class TestInternalAnchorReverse(unittest.TestCase):
    """
    Mirror of TestInternalAnchorForward: this time the CONTIG's own edge
    (in practice: the very first seed read's own raw, unverified edge --
    see denovo_seed_olc.py module docstring) carries the noise, not the
    candidate. Added because internal_anchor_extend_indexed's forward-only
    pair structurally could not fix this (it only ever appends new
    candidate content past an anchor; it never replaces contig content
    already committed before that anchor).

    Both directions require the candidate's own content to reach fully
    through to the far end being repaired (contig's end for 3', position
    0 for 5') -- otherwise there's no evidence to tell noise from genuine
    unique tail content, and the fallback correctly declines to guess.
    """
    def test_rescues_seed_with_noisy_3prime_edge(self):
        rng = random.Random(10)
        clean = rand_seq(rng, 600)
        junk_tail = rand_seq(rng, 10)          # seed's own noisy edge
        true_tail = rand_seq(rng, 10)          # what should have been there
        new_seq = rand_seq(rng, 50)            # genuine further extension
        contig = clean + junk_tail
        candidate = clean + true_tail + new_seq

        new_contig, idx = m._internal_anchor_extend_3prime_reverse_indexed(
            contig, [candidate], m._build_pool_kmer_index([candidate], 10), {0},
            min_ov=250, max_mm=0.05, seed_k=10, min_verify=60)
        self.assertIsNotNone(new_contig)
        self.assertEqual(new_contig, clean + true_tail + new_seq)

    def test_rescues_seed_with_noisy_5prime_edge(self):
        rng = random.Random(11)
        clean = rand_seq(rng, 600)
        junk_head = rand_seq(rng, 10)
        true_head = rand_seq(rng, 10)
        new_seq = rand_seq(rng, 50)
        contig = junk_head + clean
        candidate = new_seq + true_head + clean

        new_contig, idx = m._internal_anchor_extend_5prime_reverse_indexed(
            contig, [candidate], m._build_pool_kmer_index([candidate], 10), {0},
            min_ov=250, max_mm=0.05, seed_k=10, min_verify=60)
        self.assertIsNotNone(new_contig)
        self.assertEqual(new_contig, new_seq + true_head + clean)

    def test_declines_when_candidate_does_not_reach_full_end(self):
        """
        Candidate anchors inside the contig's tail but is too short to
        reach all the way to the contig's own current end -- there's no
        way to know if the un-covered remainder is noise or real, so the
        fallback must return None rather than guess.
        """
        rng = random.Random(12)
        clean = rand_seq(rng, 600)
        junk_tail = rand_seq(rng, 10)
        contig = clean + junk_tail
        short_candidate = clean[-100:]  # anchors, but doesn't reach past clean's own end

        new_contig, idx = m._internal_anchor_extend_3prime_reverse_indexed(
            contig, [short_candidate],
            m._build_pool_kmer_index([short_candidate], 10), {0},
            min_ov=250, max_mm=0.05, seed_k=10, min_verify=60)
        self.assertIsNone(new_contig)

    def test_real_barcode_end_to_end_unaffected(self):
        """
        Sanity check that wiring the reverse fallback into
        internal_anchor_extend_indexed didn't change assemble_umi's
        behavior on the case the forward fallback was originally built
        for (candidate-side noise, not seed-side) -- reverse should only
        ever fire as a last resort, after both forward directions fail.
        """
        rng = random.Random(9)
        contig = rand_seq(rng, 400)
        clean_overlap = contig[-80:]
        junk_prefix = rand_seq(rng, 30)
        new_tail = rand_seq(rng, 100)
        candidate = junk_prefix + clean_overlap + new_tail
        contigs = m.assemble_umi([contig, candidate], min_ov=20, max_mm=0.05,
                                 min_ctg=400, use_internal_anchor=True)
        self.assertEqual(len(contigs), 1)
        self.assertEqual(contigs[0], contig + new_tail)


class TestChimericMergeSafety(unittest.TestCase):
    """
    A short sequence shared by two otherwise-unrelated reads must not be
    enough, by itself, to trigger a merge -- min_verify (60bp default)
    is deliberately longer than typical short conserved motifs (e.g.
    bacterial 16S primer sites) specifically to guard against this.
    Covers both the original forward internal-anchor fallback and the
    new reverse one added in this session.
    """
    N_TRIALS = 100

    def test_shared_conserved_motif_forward(self):
        rng = random.Random(20)
        false_merges = 0
        for _ in range(self.N_TRIALS):
            motif = rand_seq(rng, 40)
            a = rand_seq(rng, 200) + motif + rand_seq(rng, 30)
            b = rand_seq(rng, 30) + motif + rand_seq(rng, 200)
            if m.suffix_prefix_overlap(a, b, min_ov=20, max_mm=0.05, seed_k=10) > 0:
                false_merges += 1
        self.assertEqual(false_merges, 0)

    def test_shared_conserved_motif_reverse(self):
        rng = random.Random(21)
        false_3p = false_5p = 0
        for _ in range(self.N_TRIALS):
            motif = rand_seq(rng, 40)
            real_a = rand_seq(rng, 300)
            real_b = rand_seq(rng, 300)

            contig3 = real_a + motif
            cand3 = motif + real_b
            new3, _ = m._internal_anchor_extend_3prime_reverse_indexed(
                contig3, [cand3], m._build_pool_kmer_index([cand3], 10), {0},
                min_ov=20, max_mm=0.05, seed_k=10, min_verify=60)
            if new3 is not None:
                false_3p += 1

            contig5 = motif + real_a
            cand5 = real_b + motif
            new5, _ = m._internal_anchor_extend_5prime_reverse_indexed(
                contig5, [cand5], m._build_pool_kmer_index([cand5], 10), {0},
                min_ov=20, max_mm=0.05, seed_k=10, min_verify=60)
            if new5 is not None:
                false_5p += 1
        self.assertEqual(false_3p, 0)
        self.assertEqual(false_5p, 0)

    def test_shared_polya_tail_is_a_known_preexisting_limitation(self):
        """
        NOT a target for this test suite to enforce as "safe": an exact
        homopolymer run at the boundary is structurally indistinguishable
        from a real overlap using sequence matching alone (no coverage
        depth or quality signal available at this layer). Verified this
        is pre-existing behavior, unchanged by any fix in this file's
        history (same result on the pre-refactor baseline). Recorded here
        so a future change that happens to alter it gets noticed and
        deliberately evaluated, rather than silently drifting.
        """
        rng = random.Random(22)
        false_merges = 0
        for _ in range(self.N_TRIALS):
            a = rand_seq(rng, 200) + "A" * 40
            b = "A" * 40 + rand_seq(rng, 200)
            if m.suffix_prefix_overlap(a, b, min_ov=20, max_mm=0.05, seed_k=10) > 0:
                false_merges += 1
        self.assertEqual(false_merges, self.N_TRIALS)


class TestMinimizerDedup(unittest.TestCase):
    """
    Regression test for a real production bug: _minimizer_dedupe's single
    global canonical minimizer per read is a much weaker signal than real
    minimizer-based dedup tools use, and on 16S data specifically an
    AT-rich conserved motif can easily BE the lexicographically smallest
    k-mer in many otherwise completely different reads that tile unrelated
    true positions and share nothing else. Confirmed at full production
    scale with this dedup enabled by default: total contigs 124,240 ->
    109,039, >=1kb contigs 21,139 -> 17,595, net length change ~-1.06 Mb
    versus the prior baseline. use_minimizer_dedup now defaults to False
    everywhere (assemble_umi, _CFG, configure(), the CLI) until the
    single-minimizer collision risk is fixed (e.g. windowed/multiple
    minimizers, or verifying actual shared sequence length before
    collapsing) and re-validated at full scale.
    """

    def test_conserved_motif_collapses_genuinely_distinct_reads(self):
        rng = random.Random(80)
        conserved_kmer = "A" * 21  # lexicographically minimal -- plausible in real AT-rich conserved stretches
        reads = []
        for _ in range(5):
            unique_part = rand_seq(rng, 500)
            reads.append(unique_part[:250] + conserved_kmer + unique_part[250:])

        minimizers = {m._canonical_minimizer(r, k=21) for r in reads}
        self.assertEqual(len(minimizers), 1,
                          "test setup assumption broken: reads no longer share one minimizer")

        deduped = m._minimizer_dedupe(reads, k=21, keep_n=2)
        self.assertEqual(len(deduped), 2,
                          "documents the known failure mode -- 3 of 5 genuinely distinct "
                          "reads are silently discarded. This is NOT a passing bar to defend; "
                          "it exists so a future redesign fixing this is a deliberate, visible "
                          "change to this assertion, not a silent behavior drift.")

    def test_defaults_to_off_everywhere(self):
        import inspect
        self.assertFalse(inspect.signature(m.assemble_umi).parameters["use_minimizer_dedup"].default)
        self.assertFalse(inspect.signature(m.configure).parameters["use_minimizer_dedup"].default)
        self.assertFalse(m._CFG["use_minimizer_dedup"])


class TestDedupeAndMergeContigs(unittest.TestCase):
    """
    _dedupe_and_merge_contigs used to disable internal-anchor entirely
    (use_internal_anchor=False) on the theory that already-assembled
    contigs have clean, verified boundaries. Real production data (a
    1.5kb-library run) showed that reasoning was incomplete: the vast
    majority of multi-contig UMIs had 200-750bp of genuine shared
    sequence between their separate contigs, sitting hundreds of bp from
    either edge -- not because either contig's edge was noisy, but
    because two independently-grown contigs can have their true
    connection point deep in one or both of them. The default read-scale
    internal-anchor search window (~60bp) structurally cannot reach that
    far, so contig-merging now passes a check_len_override covering the
    whole contig (cheap: at most a handful of contigs per UMI).
    """

    def test_merges_when_connection_point_is_mid_sequence_on_one_side(self):
        """
        Matches the geometry actually observed on real production data:
        the LONGER contig's own matching region sits at ITS OWN edge
        (like a prior verified merge or the original seed's tail), while
        the same shared region sits in the MIDDLE of the shorter
        candidate, which also has genuine unique content on the far side
        of it to contribute. The default read-scale window can't reach a
        match this deep into a several-hundred-bp contig;
        check_len_override (covering the whole contig) can.

        The longer sequence must be the one whose own remaining content
        past the anchor is empty (shared = its exact tail): pool sorts
        longest-first, so the longer one always becomes the seed here,
        and the verification window reaches "as far as possible" toward
        the seed's own far end rather than stopping exactly at the true
        shared-region boundary -- if the seed itself had unrelated
        content past the anchor, that reach would overshoot into it and
        fail verification. This is a real, narrower-than-ideal edge of
        the current fix, not a coincidence of this test.

        NOTE: the current fix does NOT yet handle the fully general case
        where the shared region reaches neither contig's edge on EITHER
        side (verified separately: real barcode "AAAAAAAAAAAAGTT",
        942bp/1045bp with a genuine 383bp shared region, still correctly
        declines to merge rather than guess) -- that needs a proper
        local-alignment splice merge, not just a wider search window.
        """
        rng = random.Random(30)
        shared = rand_seq(rng, 300)
        longer = rand_seq(rng, 300) + shared                                # shared = longer's own tail
        shorter = rand_seq(rng, 50) + shared + rand_seq(rng, 100)           # shared sits in shorter's middle
        merged = m._dedupe_and_merge_contigs([longer, shorter], min_ov=20, max_mm=0.05, seed_k=10)
        self.assertEqual(len(merged), 1)
        self.assertGreater(len(merged[0]), max(len(longer), len(shorter)))

    def test_does_not_silently_discard_a_shorter_contigs_unique_edge(self):
        """
        Regression test for a real bug found while validating the
        check_len_override fix above: forward-direction internal anchor
        can find a match between a candidate's OWN content and the
        middle of a longer contig where the match reaches all the way to
        the candidate's own end (zero new trailing content). The old
        code still accepted this as a "successful merge", consuming the
        candidate's pool slot and silently discarding whatever unique
        content it had BEFORE the anchor -- on real data this made an
        entire genuine 584bp contig vanish with no corresponding growth
        anywhere else. Fixed by requiring genuine new content before
        accepting a forward-direction match.
        """
        rng = random.Random(31)
        unique_prefix = rand_seq(rng, 100)   # exists ONLY in the shorter contig
        shared_middle = rand_seq(rng, 300)
        shorter = unique_prefix + shared_middle          # shared sits at shorter's own END
        longer = rand_seq(rng, 200) + shared_middle + rand_seq(rng, 100)  # shared sits in longer's MIDDLE

        merged = m._dedupe_and_merge_contigs([longer, shorter], min_ov=20, max_mm=0.05, seed_k=10)
        total_len = sum(len(c) for c in merged)
        # whatever the merge decides to do, unique_prefix's content must
        # not simply vanish -- either `shorter` survives untouched (safe,
        # conservative fallback) or some contig grew enough to account
        # for it. What must NOT happen: `longer` reappears completely
        # unchanged as the ONLY output with `shorter` gone and no growth.
        self.assertFalse(
            len(merged) == 1 and merged[0] == longer,
            "shorter contig's unique prefix was silently discarded with zero growth"
        )


class TestMultiprocessingConfigPropagation(unittest.TestCase):
    """
    Regression test for a real bug found while reviewing this file:
    _init_pool_worker only forwarded meta_data1/meta_data2/lock, not
    _CFG. Under 'fork' (Linux default) this silently worked because a
    forked worker inherits the parent's already-configure()'d _CFG via
    copy-on-write memory -- but under 'spawn' (macOS/Windows default),
    each worker re-imports the module fresh, resetting _CFG to its
    hardcoded defaults and discarding every configure() call the parent
    made. Reproduced concretely before the fix: with configure() called
    in the parent only, spawned workers tried to write
    'denovo/final_contigs_0.fa' (the hardcoded default) instead of the
    configured path. This test runs actual multiprocessing (not just a
    syntax check) and would fail again if that regressed.
    """
    def test_single_vs_multi_process_output_identical(self):
        meta = {
            "BC_A": [">r0", "ACGT" * 100, ">r1", "ACGT" * 100],
            "BC_B": [">r0", "TTTT" * 100, ">r1", "TTTT" * 100],
        }

        tmp1 = tempfile.mkdtemp()
        out1 = os.path.join(tmp1, "out_{id}.fa")
        m.configure(min_ctg_len=100, out_id=0, out_file=out1)
        m._process_se_metadata(dict(meta), num_processes=1)
        with open(out1.format(id=0)) as f:
            single_process_output = f.read()

        tmp2 = tempfile.mkdtemp()
        out2 = os.path.join(tmp2, "out_{id}.fa")
        m.configure(min_ctg_len=100, out_id=0, out_file=out2)
        m._process_se_metadata(dict(meta), num_processes=2)
        with open(out2.format(id=0)) as f:
            multi_process_output = f.read()

        self.assertEqual(sorted(single_process_output.split()),
                          sorted(multi_process_output.split()))

    def test_configured_min_ctg_propagates_to_spawned_workers(self):
        """
        An absurdly high min_ctg means every contig should be filtered
        out -- if a worker silently fell back to the default min_ctg
        (400), the short test contigs here would pass and get written,
        proving the parent's configure() call was lost.
        """
        meta = {
            "BC_A": [">r0", "ACGT" * 100, ">r1", "ACGT" * 100],
            "BC_B": [">r0", "TTTT" * 100, ">r1", "TTTT" * 100],
        }
        tmp = tempfile.mkdtemp()
        out_file = os.path.join(tmp, "out_{id}.fa")
        m.configure(min_ctg_len=999999, out_id=0, out_file=out_file)
        m._process_se_metadata(dict(meta), num_processes=2)
        path = out_file.format(id=0)
        content = open(path).read() if os.path.exists(path) else ""
        self.assertEqual(content.strip(), "")


if __name__ == "__main__":
    unittest.main(verbosity=2)
