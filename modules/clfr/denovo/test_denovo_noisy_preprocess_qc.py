#!/usr/bin/env python3
"""Tests for denovo_noisy_preprocess_qc.py."""

import os
import tempfile
import unittest

import denovo_noisy_preprocess_qc as qc


def reads(depths, short_by_barcode=None):
    short_by_barcode = short_by_barcode or {}
    for i, depth in enumerate(depths):
        barcode = "BC{:04d}".format(i)
        n_short = short_by_barcode.get(i, 0)
        for j in range(depth):
            yield barcode, ("A" * (200 if j < n_short else 500))


def zymo_like_baseline():
    return {
        "depth": {
            "uppers": qc.DEPTH_BIN_UPPERS,
            "counts": [0, 0, 100, 0, 0, 0, 0, 0, 0],
            "top1pct_read_fraction": 0.03,
        },
        "length": {
            "uppers": qc.LENGTH_BIN_UPPERS,
            "counts": [0, 0, 600, 0, 0, 5400, 0, 0],
            "short_read_fraction": 0.10,
        },
    }


class NoisyPreprocessQCTest(unittest.TestCase):
    def test_healthy_sample_is_pass_through_candidate(self):
        result = qc.summarize(reads([60] * 100),
                              baseline=zymo_like_baseline())
        self.assertEqual(result["candidate_verdict"], "pass_through")
        self.assertFalse(result["depth_adverse_drift"])
        self.assertFalse(result["read_length_adverse_drift"])
        self.assertEqual(result["projected_depth_median"], 60)

    def test_hs1_like_sample_is_salvage_candidate(self):
        depths = [5000] + [1000] * 9 + [60] * 90
        short = {0: 4500}
        short.update({i: 900 for i in range(1, 10)})
        result = qc.summarize(reads(depths, short),
                              baseline=zymo_like_baseline())
        self.assertEqual(result["candidate_verdict"], "salvage_candidate")
        self.assertTrue(result["depth_adverse_drift"])
        self.assertTrue(result["read_length_adverse_drift"])
        self.assertTrue(result["projected_postfilter_safe"])

    def test_short_low_depth_sample_is_not_salvaged(self):
        depths = [30] * 100
        short = {i: 10 for i in range(100)}
        result = qc.summarize(reads(depths, short),
                              baseline=zymo_like_baseline())
        self.assertEqual(result["candidate_verdict"], "drift_no_salvage")
        self.assertFalse(result["depth_adverse_drift"])
        self.assertTrue(result["read_length_adverse_drift"])
        self.assertFalse(result["projected_postfilter_safe"])

    def test_hs2_like_higher_depth_is_not_adverse_drift(self):
        short = {i: 9 for i in range(100)}
        result = qc.summarize(reads([75] * 100, short),
                              baseline=zymo_like_baseline())
        self.assertTrue(result["depth_distribution_drift"])
        self.assertFalse(result["depth_adverse_drift"])
        self.assertFalse(result["read_length_adverse_drift"])
        self.assertEqual(result["candidate_verdict"], "pass_through")

    def test_report_is_explicitly_non_acting(self):
        result = qc.summarize(reads([60] * 10),
                              baseline=zymo_like_baseline())
        with tempfile.TemporaryDirectory() as tmp:
            out = os.path.join(tmp, "report.tsv")
            qc.write_report(out, result, "reads.tsv", "tsv", 300, 300,
                            mode="report", action="report_only")
            with open(out) as fh:
                report = dict(line.rstrip("\n").split("\t", 1)
                              for line in fh if "\t" in line)
        self.assertEqual(report["mode"], "report")
        self.assertEqual(report["action"], "report_only")

    def test_empty_input_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "no reads"):
            qc.summarize(iter(()))

    def test_salvage_keeps_only_long_reads_and_caps_each_umi(self):
        with tempfile.TemporaryDirectory() as tmp:
            source = os.path.join(tmp, "reads.tsv")
            output = os.path.join(tmp, "salvaged.tsv")
            with open(source, "w") as out:
                for i in range(4):
                    length = 200 if i == 0 else 500
                    out.write("BX:Z:A @r{}\t{}\n".format(i, "A" * length))
                for i in range(2):
                    out.write("BX:Z:B @s{}\t{}\n".format(i, "C" * 500))
            kept = qc.filter_tsv(source, output, 300, 2)
            with open(output) as inp:
                lines = inp.readlines()
        self.assertEqual(kept, 4)
        self.assertEqual(len(lines), 4)
        self.assertTrue(all(len(line.rstrip().split("\t")[1]) >= 300
                            for line in lines))


if __name__ == "__main__":
    unittest.main()
