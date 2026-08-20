#!/usr/bin/env python3
"""Unit tests for denovo_candidate_select.py.

This is the component that decides which candidate contig a UMI actually
delivers, so the tests below focus on what matters for a delivery decision:
the default mode must reproduce the historical "always k41_0" behaviour
byte-for-byte, gated_switch must only ever switch into a candidate that
genuinely passes the gate, and a candidate with no QC row (a genuinely
missing signal, not a signal that failed) must never be treated as a pass.

Run:
    python3 -m unittest test_denovo_candidate_select -v
"""
import csv
import os
import subprocess
import sys
import tempfile
import unittest

SCRIPT = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "denovo_candidate_select.py")


class SelectCase(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()

    def _write_fasta(self, barcodes_to_seqs):
        """barcodes_to_seqs: {barcode: [seq_rank0, seq_rank1, ...]}."""
        path = os.path.join(self.tmp, "contigs.fa")
        with open(path, "w") as fh:
            for barcode, seqs in barcodes_to_seqs.items():
                for rank, seq in enumerate(seqs):
                    fh.write(f">{barcode}>k41_{rank}\n{seq}\n")
        return path

    def _write_qc(self, rows):
        """rows: list of dicts with barcode, k41_rank, span_cov_ratio,
        placed_reads (contig_len filled in if absent)."""
        path = os.path.join(self.tmp, "qc.tsv")
        fields = ["barcode", "header", "k41_rank", "contig_len",
                  "placed_reads", "span_cov_ratio", "min_local_span_ratio"]
        with open(path, "w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
            writer.writeheader()
            for row in rows:
                full = {"header": f"{row['barcode']}>k41_{row['k41_rank']}",
                        "contig_len": 100, "min_local_span_ratio": 0.1}
                full.update(row)
                writer.writerow(full)
        return path

    def _run(self, contigs, mode, candidate_qc=None, extra=None):
        out_fasta = os.path.join(self.tmp, "out.fa")
        out_report = os.path.join(self.tmp, "report.tsv")
        out_decision = os.path.join(self.tmp, "decision.tsv")
        cmd = [sys.executable, SCRIPT,
               "--contigs", contigs, "--mode", mode,
               "--out-fasta", out_fasta, "--out-report", out_report,
               "--out-decision", out_decision]
        if candidate_qc:
            cmd += ["--candidate-qc", candidate_qc]
        cmd += extra or []
        proc = subprocess.run(cmd, capture_output=True, text=True)
        self.assertEqual(proc.returncode, 0, proc.stderr)

        fasta = {}
        header = None
        with open(out_fasta) as fh:
            for line in fh:
                line = line.rstrip("\n")
                if line.startswith(">"):
                    header = line[1:]
                else:
                    fasta[header] = line

        report = {}
        with open(out_report) as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                report[row["barcode"]] = row

        return fasta, report


class TestLongestMode(SelectCase):
    def test_always_picks_rank0_even_with_no_qc(self):
        contigs = self._write_fasta({"BC1": ["AAAA", "CCCCCC"]})
        fasta, report = self._run(contigs, mode="longest")
        self.assertEqual(fasta["BC1>k41_0"], "AAAA")
        self.assertEqual(report["BC1"]["chosen_rank"], "0")
        self.assertEqual(report["BC1"]["switched"], "0")

    def test_matches_historical_awk_behavior_on_multiple_barcodes(self):
        contigs = self._write_fasta({
            "BC1": ["AAAA", "CCCC"],
            "BC2": ["GGGGGGGG"],
        })
        fasta, report = self._run(contigs, mode="longest")
        self.assertEqual(set(fasta), {"BC1>k41_0", "BC2>k41_0"})
        self.assertEqual(report["BC1"]["switched"], "0")
        self.assertEqual(report["BC2"]["switched"], "0")


class TestGatedSwitchMode(SelectCase):
    def test_switches_to_first_qualifying_candidate(self):
        contigs = self._write_fasta({"BC1": ["AAAA", "CCCC", "GGGG"]})
        qc = self._write_qc([
            {"barcode": "BC1", "k41_rank": 0, "span_cov_ratio": 0.1, "placed_reads": 5},
            {"barcode": "BC1", "k41_rank": 1, "span_cov_ratio": 0.4, "placed_reads": 3},
            {"barcode": "BC1", "k41_rank": 2, "span_cov_ratio": 0.9, "placed_reads": 9},
        ])
        fasta, report = self._run(contigs, mode="gated_switch", candidate_qc=qc)
        # rank 0 fails the span gate; rank 1 is the first to pass -> switch there,
        # not all the way to rank 2 which is also a pass but later in rank order.
        self.assertEqual(fasta["BC1>k41_1"], "CCCC")
        self.assertEqual(report["BC1"]["chosen_rank"], "1")
        self.assertEqual(report["BC1"]["switched"], "1")

    def test_falls_back_to_primary_when_nothing_passes_gate(self):
        contigs = self._write_fasta({"BC1": ["AAAA", "CCCC"]})
        qc = self._write_qc([
            {"barcode": "BC1", "k41_rank": 0, "span_cov_ratio": 0.1, "placed_reads": 1},
            {"barcode": "BC1", "k41_rank": 1, "span_cov_ratio": 0.1, "placed_reads": 1},
        ])
        fasta, report = self._run(contigs, mode="gated_switch", candidate_qc=qc)
        self.assertEqual(fasta["BC1>k41_0"], "AAAA")
        self.assertEqual(report["BC1"]["chosen_rank"], "0")
        self.assertEqual(report["BC1"]["switched"], "0")
        self.assertEqual(report["BC1"]["reason"], "no_candidate_passed_gate")

    def test_primary_itself_passing_gate_is_not_counted_as_switched(self):
        contigs = self._write_fasta({"BC1": ["AAAA", "CCCC"]})
        qc = self._write_qc([
            {"barcode": "BC1", "k41_rank": 0, "span_cov_ratio": 0.9, "placed_reads": 9},
            {"barcode": "BC1", "k41_rank": 1, "span_cov_ratio": 0.9, "placed_reads": 9},
        ])
        fasta, report = self._run(contigs, mode="gated_switch", candidate_qc=qc)
        self.assertEqual(report["BC1"]["chosen_rank"], "0")
        self.assertEqual(report["BC1"]["switched"], "0")

    def test_candidate_missing_a_qc_row_is_never_selected(self):
        # rank 1 has no QC row at all (e.g. denovo_junction_qc.py's analyze()
        # returned None for it, zero placed reads) -- must be treated as a
        # gate failure, never as an implicit pass.
        contigs = self._write_fasta({"BC1": ["AAAA", "CCCC", "GGGG"]})
        qc = self._write_qc([
            {"barcode": "BC1", "k41_rank": 0, "span_cov_ratio": 0.1, "placed_reads": 1},
            {"barcode": "BC1", "k41_rank": 2, "span_cov_ratio": 0.9, "placed_reads": 9},
        ])
        fasta, report = self._run(contigs, mode="gated_switch", candidate_qc=qc)
        self.assertEqual(report["BC1"]["chosen_rank"], "2")

    def test_respects_custom_thresholds(self):
        contigs = self._write_fasta({"BC1": ["AAAA", "CCCC"]})
        qc = self._write_qc([
            {"barcode": "BC1", "k41_rank": 0, "span_cov_ratio": 0.3, "placed_reads": 2},
            {"barcode": "BC1", "k41_rank": 1, "span_cov_ratio": 0.3, "placed_reads": 2},
        ])
        # default thresholds (span>=0.25, placed>=2) would pass rank 0 itself
        fasta, report = self._run(contigs, mode="gated_switch", candidate_qc=qc)
        self.assertEqual(report["BC1"]["chosen_rank"], "0")
        # a stricter min-placed-reads should now reject both, falling back
        fasta, report = self._run(contigs, mode="gated_switch", candidate_qc=qc,
                                   extra=["--min-placed-reads", "5"])
        self.assertEqual(report["BC1"]["reason"], "no_candidate_passed_gate")

    def test_requires_candidate_qc(self):
        contigs = self._write_fasta({"BC1": ["AAAA"]})
        out_fasta = os.path.join(self.tmp, "out.fa")
        cmd = [sys.executable, SCRIPT, "--contigs", contigs, "--mode", "gated_switch",
               "--out-fasta", out_fasta,
               "--out-report", os.path.join(self.tmp, "r.tsv"),
               "--out-decision", os.path.join(self.tmp, "d.tsv")]
        proc = subprocess.run(cmd, capture_output=True, text=True)
        self.assertNotEqual(proc.returncode, 0)
        self.assertIn("--candidate-qc", proc.stderr)


class TestDecisionRecord(SelectCase):
    def test_decision_records_settings_and_switch_count(self):
        contigs = self._write_fasta({
            "BC1": ["AAAA", "CCCC"],
            "BC2": ["GGGG", "TTTT"],
        })
        qc = self._write_qc([
            {"barcode": "BC1", "k41_rank": 0, "span_cov_ratio": 0.1, "placed_reads": 1},
            {"barcode": "BC1", "k41_rank": 1, "span_cov_ratio": 0.9, "placed_reads": 9},
            {"barcode": "BC2", "k41_rank": 0, "span_cov_ratio": 0.9, "placed_reads": 9},
            {"barcode": "BC2", "k41_rank": 1, "span_cov_ratio": 0.9, "placed_reads": 9},
        ])
        self._run(contigs, mode="gated_switch", candidate_qc=qc)
        decision = {}
        with open(os.path.join(self.tmp, "decision.tsv")) as fh:
            fh.readline()
            for line in fh:
                k, _, v = line.rstrip("\n").partition("\t")
                decision[k] = v
        self.assertEqual(decision["mode"], "gated_switch")
        self.assertEqual(decision["barcodes_total"], "2")
        self.assertEqual(decision["barcodes_switched"], "1")


if __name__ == "__main__":
    unittest.main()
