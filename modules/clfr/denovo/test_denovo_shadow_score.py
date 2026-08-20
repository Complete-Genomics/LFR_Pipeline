import csv
import importlib.util
import sys
import tempfile
import types
import unittest
from pathlib import Path
from unittest.mock import patch


MODULE_PATH = Path(__file__).with_name("denovo_shadow_score.py")
SPEC = importlib.util.spec_from_file_location("denovo_shadow_score", MODULE_PATH)
shadow_score = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(shadow_score)


class ShadowScoreTest(unittest.TestCase):
    def test_scores_all_candidates_and_reports_rule_disagreement(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            candidate_qc = tmp / "candidate_qc.tsv"
            selection = tmp / "selection.tsv"
            candidate_out = tmp / "shadow" / "candidate_scores.tsv"
            summary_out = tmp / "shadow" / "umi_summary.tsv"
            candidate_qc.write_text(
                "barcode\theader\tk41_rank\tcontig_len\tplaced_reads\t"
                "span_cov_ratio\tmin_local_span_ratio\n"
                "bc1\tbc1>k41_0\t0\t1000\t10\t0.30\t0.10\n"
                "bc1\tbc1>k41_1\t1\t800\t8\t0.40\t0.20\n"
                "bc2\tbc2>k41_0\t0\t900\t6\t0.35\t0.15\n"
            )
            selection.write_text(
                "barcode\tchosen_rank\tn_candidates\tswitched\treason\n"
                "bc1\t0\t2\t0\tmode=longest\n"
                "bc2\t0\t1\t0\tmode=longest\n"
            )

            class FakeBooster:
                def __init__(self, model_file):
                    self.model_file = model_file

                def predict(self, features):
                    return [0.9 if int(row[5]) == 0 else 0.1 for row in features]

            fake_lgb = types.SimpleNamespace(Booster=FakeBooster)
            argv = [
                "denovo_shadow_score.py",
                "--candidate-qc", str(candidate_qc),
                "--selection", str(selection),
                "--model", "fake.lgb",
                "--out-candidates", str(candidate_out),
                "--out-summary", str(summary_out),
            ]
            with patch.dict(sys.modules, {"lightgbm": fake_lgb}), \
                    patch.object(sys, "argv", argv):
                shadow_score.main()

            with candidate_out.open() as fh:
                candidate_rows = list(csv.DictReader(fh, delimiter="\t"))
            with summary_out.open() as fh:
                summary_rows = list(csv.DictReader(fh, delimiter="\t"))

            self.assertEqual(len(candidate_rows), 3)
            self.assertEqual(candidate_rows[1]["selected_by_rule"], "0")
            self.assertEqual(summary_rows[0]["lowest_p_chimera_rank"], "1")
            self.assertEqual(summary_rows[0]["gbdt_disagrees_with_rule"], "1")
            self.assertEqual(summary_rows[1]["gbdt_disagrees_with_rule"], "0")


if __name__ == "__main__":
    unittest.main()
