"""Focused tests for the shadow-only cLFR resource policy."""

from __future__ import annotations

import csv
import sys
import tempfile
import unittest
from pathlib import Path


MODULE_DIR = Path(__file__).resolve().parents[1] / "modules" / "shared" / "resource_allocation"
sys.path.insert(0, str(MODULE_DIR))
import resource_allocation as allocation  # noqa: E402


FIELDS = [
    "task_key", "run_id", "rule", "input_size_bytes", "static_memory_gb",
    "static_runtime_min", "max_rss_bytes", "runtime_sec", "status",
    "workflow_version", "node_class",
]


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str] = FIELDS) -> None:
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


class ResourceAllocationTest(unittest.TestCase):
    def test_rule_bucket_upper_bound_can_reduce_only_with_explicit_floor(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            history_path = root / "history.tsv"
            candidates_path = root / "candidates.tsv"
            output_path = root / "shadow.tsv"
            history = [
                {
                    "task_key": f"old-{index}", "run_id": f"r{index}", "rule": "map_reads_minimap",
                    "input_size_bytes": str(2 ** 30), "static_memory_gb": "96", "static_runtime_min": "600",
                    "max_rss_bytes": str((20 + index) * allocation.GIB), "runtime_sec": str((100 + index) * 60),
                    "status": "success", "workflow_version": "v1", "node_class": "c40",
                }
                for index in range(3)
            ]
            write_tsv(history_path, history)
            candidate = dict(history[0])
            candidate.update({"task_key": "new", "run_id": "new-run", "memory_floor_gb": "32", "runtime_floor_min": "60"})
            write_tsv(candidates_path, [candidate], FIELDS + ["memory_floor_gb", "runtime_floor_min"])

            allocation.shadow(
                history_path, candidates_path, output_path, 0.95, 1.30, 1.25, 3, 3
            )
            row = allocation.read_tsv(output_path)[0]
            self.assertEqual(row["reason"], "quantile_upper_bound")
            self.assertEqual(row["stratum"], "rule_size_2pow_30")
            self.assertLess(float(row["suggested_memory_gb"]), 96)

    def test_unknown_context_and_ood_input_fall_back_to_static(self):
        history = [
            {
                "task_key": "old", "run_id": "r1", "rule": "get_consensus_fasta",
                "input_size_bytes": "1000", "static_memory_gb": "64", "static_runtime_min": "240",
                "max_rss_bytes": str(10 * allocation.GIB), "runtime_sec": "1200", "status": "success",
                "workflow_version": "v1", "node_class": "c40",
            }
        ]
        candidate = dict(history[0])
        candidate.update({"task_key": "new", "workflow_version": "v2"})
        result = allocation.recommend_one(candidate, history, 0.95, 1.3, 1.25, 1, 1)
        self.assertEqual(result["reason"], "static_fallback_unknown_context")

        candidate["workflow_version"] = "v1"
        candidate["input_size_bytes"] = "2000"
        result = allocation.recommend_one(candidate, history, 0.95, 1.3, 1.25, 1, 1)
        self.assertEqual(result["reason"], "static_fallback_ood_input_size")

    def test_collect_reads_snakemake_benchmark(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            benchmark = root / "job.txt"
            write_tsv(benchmark, [{"s": "30.5", "max_rss": "2048"}], ["s", "max_rss"])
            manifest = root / "manifest.tsv"
            row = {
                "task_key": "r1:map", "run_id": "r1", "rule": "map_reads_minimap",
                "input_size_bytes": "100", "static_memory_gb": "64", "static_runtime_min": "20",
                "benchmark_path": str(benchmark),
            }
            write_tsv(manifest, [row], list(row))
            output = root / "history.tsv"
            allocation.collect_history(manifest, output, "kb")
            parsed = allocation.read_tsv(output)[0]
            self.assertEqual(parsed["max_rss_bytes"], str(2048 * 1024))
            self.assertEqual(parsed["runtime_sec"], "30.5")

    def test_prepare_candidate_reads_metadata_only(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            input_path = root / "input.bam"
            input_path.write_bytes(b"a" * 17)
            output = root / "candidate.tsv"
            allocation.prepare_candidate(
                "r1:sort", "r1", "sort_reformated_bam", [str(input_path)], output,
                32, 60, workflow_version="v1", node_class="c40", umi_groups=5,
            )
            row = allocation.read_tsv(output)[0]
            self.assertEqual(row["input_size_bytes"], "17")
            self.assertEqual(row["memory_floor_gb"], "32")
            self.assertEqual(row["umi_groups"], "5")


if __name__ == "__main__":
    unittest.main()
