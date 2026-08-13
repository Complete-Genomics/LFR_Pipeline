#!/usr/bin/env python3
"""Unit tests for denovo_read_features.py.

Run:
    python3 -m unittest test_denovo_read_features -v

Pure stdlib (unittest + random + multiprocessing), matching the convention in
test_denovo_seed_olc.py.
"""
import gzip
import os
import random
import subprocess
import sys
import tempfile
import unittest
from collections import deque

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import denovo_read_features as m


def _reference_minimizer_gaps(seq, k=m.K_MINI, w=m.W_MINI):
    """The straightforward O(n*w) definition the optimized version replaced.

    Kept here as the test oracle: the fast version's whole justification is
    that it computes exactly this, so the property worth asserting is
    equivalence to the definition, not any hand-picked expected output.
    """
    n = len(seq)
    if n < k:
        return []
    kmers = [(seq[i:i + k], i) for i in range(n - k + 1)]
    if len(kmers) < w:
        return []
    positions = []
    last = None
    for i in range(0, len(kmers) - w + 1):
        _, mpos = min(kmers[i:i + w], key=lambda x: x[0])
        if mpos != last:
            positions.append(mpos)
            last = mpos
    positions = sorted(set(positions))
    return [positions[i + 1] - positions[i] for i in range(len(positions) - 1)]


def _write_fastq(path, records):
    with open(path, "w") as fh:
        for rid, seq, qual in records:
            fh.write("@{}\tBX:Z:{}\n{}\n+\n{}\n".format(
                rid, rid.split("#")[1].split("/")[0], seq, qual))


def _rand_seq(rng, n):
    return "".join(rng.choice("ACGT") for _ in range(n))


class TestMinimizerGaps(unittest.TestCase):
    def test_matches_naive_definition_on_random_sequence(self):
        rng = random.Random(0)
        for _ in range(200):
            seq = _rand_seq(rng, rng.randint(1, 400))
            self.assertEqual(m.minimizer_gaps(seq),
                             _reference_minimizer_gaps(seq))

    def test_matches_naive_definition_on_homopolymer_runs(self):
        # Ties everywhere -- the deque's strict '>' pop is what keeps the
        # leftmost of equal k-mers, matching min()'s first-wins behavior. A
        # '>=' there would silently pick a different position on exactly this
        # input, which is common in real reads (denovo.md sec 40/41).
        rng = random.Random(1)
        for _ in range(100):
            seq = "".join(rng.choice("AACCGGTT") for _ in range(rng.randint(20, 200)))
            self.assertEqual(m.minimizer_gaps(seq),
                             _reference_minimizer_gaps(seq))
        self.assertEqual(m.minimizer_gaps("A" * 100),
                         _reference_minimizer_gaps("A" * 100))

    def test_too_short_returns_empty(self):
        self.assertEqual(m.minimizer_gaps("ACGT"), [])
        # long enough for a k-mer but not for a full window
        self.assertEqual(m.minimizer_gaps("A" * (m.K_MINI + 2)), [])


class TestPoolFeatures(unittest.TestCase):
    def test_popular_frac_needs_a_second_read_to_share_the_kmer(self):
        seq = "ACGT" * 30
        alone = m.features_for_pool("BC", [("r1#BC/2", seq, "I" * len(seq))])
        self.assertEqual(alone[0].split("\t")[m.FIELDS.index("pool_kmer_popular_frac")],
                         "0.0000")
        shared = m.features_for_pool("BC", [("r1#BC/2", seq, "I" * len(seq)),
                                            ("r2#BC/2", seq, "I" * len(seq))])
        self.assertEqual(shared[0].split("\t")[m.FIELDS.index("pool_kmer_popular_frac")],
                         "1.0000")

    def test_pool_size_reflects_the_whole_group(self):
        reads = [("r%d#BC/2" % i, "ACGT" * 30, "I" * 120) for i in range(5)]
        rows = m.features_for_pool("BC", reads)
        self.assertEqual(len(rows), 5)
        for row in rows:
            self.assertEqual(row.split("\t")[m.FIELDS.index("pool_size")], "5")

    def test_row_has_exactly_the_declared_columns(self):
        rows = m.features_for_pool("BC", [("r1#BC/2", "ACGT" * 30, "I" * 120)])
        self.assertEqual(len(rows[0].split("\t")), len(m.FIELDS))


class TestFastqGrouping(unittest.TestCase):
    def test_read_id_stops_at_the_tab(self):
        """The sgrep-style header carries a trailing '\\tBX:Z:<barcode>'.
        Taking the whole line as the read id shifted every column in the
        prototype's output -- an actual bug, hence an explicit test.
        """
        tmp = tempfile.mkdtemp()
        path = os.path.join(tmp, "in.fastq")
        _write_fastq(path, [("r1#AAAAAAAAAAAAAAA/2", "ACGT" * 30, "I" * 120)])
        groups = list(m.iter_fastq_barcode_groups(path))
        self.assertEqual(len(groups), 1)
        bc, reads = groups[0]
        self.assertEqual(bc, "AAAAAAAAAAAAAAA")
        self.assertEqual(reads[0][0], "r1#AAAAAAAAAAAAAAA/2")

    def test_groups_are_split_on_barcode_change(self):
        tmp = tempfile.mkdtemp()
        path = os.path.join(tmp, "in.fastq")
        _write_fastq(path, [
            ("r1#AAAAAAAAAAAAAAA/2", "ACGT" * 30, "I" * 120),
            ("r2#AAAAAAAAAAAAAAA/2", "ACGT" * 30, "I" * 120),
            ("r3#CCCCCCCCCCCCCCC/2", "ACGT" * 30, "I" * 120),
        ])
        groups = list(m.iter_fastq_barcode_groups(path))
        self.assertEqual([(bc, len(r)) for bc, r in groups],
                         [("AAAAAAAAAAAAAAA", 2), ("CCCCCCCCCCCCCCC", 1)])

    def test_reads_gzipped_input(self):
        tmp = tempfile.mkdtemp()
        plain = os.path.join(tmp, "in.fastq")
        _write_fastq(plain, [("r1#AAAAAAAAAAAAAAA/2", "ACGT" * 30, "I" * 120)])
        gz = plain + ".gz"
        with open(plain, "rb") as src, gzip.open(gz, "wb") as dst:
            dst.write(src.read())
        self.assertEqual(len(list(m.iter_fastq_barcode_groups(gz))), 1)


class TestParallelEquivalence(unittest.TestCase):
    """The model consumes these columns positionally, so a worker-count-
    dependent difference (row order, or a value) would silently change
    predictions rather than fail loudly. Assert byte-identical output.
    """

    def _run(self, fastq, out, nproc):
        script = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                              "denovo_read_features.py")
        subprocess.check_call(
            [sys.executable, script, "--fastq", fastq, "--out", out,
             "--num_processes", str(nproc), "--batch-barcodes", "7"],
            stderr=subprocess.DEVNULL)

    def test_single_and_multi_process_output_identical(self):
        rng = random.Random(7)
        tmp = tempfile.mkdtemp()
        fastq = os.path.join(tmp, "in.fastq")
        records = []
        # deliberately uneven pool depths, and >batch-barcodes groups, so the
        # run crosses a batch boundary with a heavy-tailed cost distribution
        for b in range(20):
            bc = "".join(rng.choice("ACGT") for _ in range(15))
            for i in range(rng.randint(1, 12)):
                seq = _rand_seq(rng, rng.randint(60, 300))
                records.append(("r%d_%d#%s/2" % (b, i, bc), seq, "I" * len(seq)))
        _write_fastq(fastq, records)

        one = os.path.join(tmp, "p1.tsv")
        four = os.path.join(tmp, "p4.tsv")
        self._run(fastq, one, 1)
        self._run(fastq, four, 4)
        with open(one) as a, open(four) as b:
            self.assertEqual(a.read(), b.read())


if __name__ == "__main__":
    unittest.main(verbosity=2)
