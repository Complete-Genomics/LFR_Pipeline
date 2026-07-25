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
            contig, [candidate], {0},
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
            contig, [candidate], {0},
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
            contig, [short_candidate], {0},
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
                contig3, [cand3], {0}, min_ov=20, max_mm=0.05, seed_k=10, min_verify=60)
            if new3 is not None:
                false_3p += 1

            contig5 = motif + real_a
            cand5 = real_b + motif
            new5, _ = m._internal_anchor_extend_5prime_reverse_indexed(
                contig5, [cand5], {0}, min_ov=20, max_mm=0.05, seed_k=10, min_verify=60)
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
