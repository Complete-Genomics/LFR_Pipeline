"""Build cLFR resource history and emit conservative shadow recommendations.

This intentionally does not alter Snakemake resources.  It operates on completed
benchmark records and a pre-submit candidate manifest so that recommendations can
be audited before any canary is considered.
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path
from typing import Iterable


GIB = 1024 ** 3
SUCCESS_STATES = {"", "success", "completed", "complete"}
CONTEXT_COLUMNS = ("workflow_version", "tool_version", "reference_id", "node_class")
REQUIRED_MANIFEST_COLUMNS = {
    "task_key",
    "run_id",
    "rule",
    "input_size_bytes",
    "static_memory_gb",
    "static_runtime_min",
}


def read_tsv(path: str | Path) -> list[dict[str, str]]:
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: str | Path, rows: Iterable[dict[str, object]], fieldnames: list[str]) -> None:
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    with open(destination, "w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def as_float(row: dict[str, str], field: str, default: float | None = None) -> float | None:
    value = row.get(field, "")
    if value in (None, ""):
        return default
    return float(value)


def runtime_seconds(row: dict[str, str]) -> float | None:
    seconds = as_float(row, "s")
    if seconds is not None:
        return seconds
    value = row.get("h:m:s", "")
    if not value:
        return None
    hours, minutes, seconds_text = value.split(":")
    return int(hours) * 3600 + int(minutes) * 60 + float(seconds_text)


def rss_multiplier(unit: str) -> float:
    return {"bytes": 1.0, "kb": 1024.0, "mb": 1024.0 ** 2, "gb": float(GIB)}[unit]


def validate_manifest(rows: list[dict[str, str]]) -> None:
    if not rows:
        raise ValueError("resource manifest has no rows")
    missing = REQUIRED_MANIFEST_COLUMNS.difference(rows[0])
    if missing:
        raise ValueError("resource manifest is missing columns: " + ", ".join(sorted(missing)))


def collect_history(
    manifest_path: str | Path,
    output_path: str | Path,
    max_rss_unit: str = "bytes",
) -> None:
    manifest = read_tsv(manifest_path)
    validate_manifest(manifest)
    rows: list[dict[str, object]] = []
    multiplier = rss_multiplier(max_rss_unit)

    for metadata in manifest:
        benchmark_path = metadata.get("benchmark_path", "")
        if not benchmark_path:
            raise ValueError(f"task_key={metadata['task_key']} is missing benchmark_path")
        benchmark_rows = read_tsv(benchmark_path)
        if len(benchmark_rows) != 1:
            raise ValueError(f"benchmark must contain exactly one row: {benchmark_path}")
        benchmark = benchmark_rows[0]
        max_rss = as_float(benchmark, "max_rss")
        elapsed = runtime_seconds(benchmark)
        if max_rss is None or elapsed is None:
            raise ValueError(f"benchmark lacks max_rss or runtime: {benchmark_path}")
        merged: dict[str, object] = dict(metadata)
        merged["max_rss_bytes"] = int(max_rss * multiplier)
        merged["runtime_sec"] = elapsed
        rows.append(merged)

    fieldnames = list(dict.fromkeys(list(manifest[0]) + ["max_rss_bytes", "runtime_sec"]))
    write_tsv(output_path, rows, fieldnames)


def prepare_candidate(
    task_key: str,
    run_id: str,
    rule: str,
    input_paths: list[str],
    output_path: str | Path,
    static_memory_gb: float,
    static_runtime_min: float,
    memory_floor_gb: float | None = None,
    runtime_floor_min: float | None = None,
    workflow_version: str = "",
    tool_version: str = "",
    reference_id: str = "",
    node_class: str = "",
    input_records: int | None = None,
    umi_groups: int | None = None,
    qualifying_umi: int | None = None,
) -> None:
    paths = [Path(path) for path in input_paths]
    missing = [str(path) for path in paths if not path.exists()]
    if missing:
        raise FileNotFoundError("candidate inputs do not exist: " + ", ".join(missing))
    sizes = [path.stat().st_size for path in paths]
    row: dict[str, object] = {
        "task_key": task_key,
        "run_id": run_id,
        "rule": rule,
        "input_size_bytes": sum(sizes),
        "input_file_count": len(paths),
        "largest_input_bytes": max(sizes),
        "static_memory_gb": static_memory_gb,
        "static_runtime_min": static_runtime_min,
        "memory_floor_gb": memory_floor_gb if memory_floor_gb is not None else static_memory_gb,
        "runtime_floor_min": runtime_floor_min if runtime_floor_min is not None else static_runtime_min,
        "workflow_version": workflow_version,
        "tool_version": tool_version,
        "reference_id": reference_id,
        "node_class": node_class,
        "input_records": input_records if input_records is not None else "",
        "umi_groups": umi_groups if umi_groups is not None else "",
        "qualifying_umi": qualifying_umi if qualifying_umi is not None else "",
    }
    write_tsv(output_path, [row], list(row))


def quantile(values: list[float], q: float) -> float:
    if not values:
        raise ValueError("cannot calculate a quantile of an empty list")
    ordered = sorted(values)
    if len(ordered) == 1:
        return ordered[0]
    position = (len(ordered) - 1) * q
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return ordered[lower]
    return ordered[lower] + (ordered[upper] - ordered[lower]) * (position - lower)


def size_bucket(input_size_bytes: float) -> int:
    return int(math.floor(math.log2(max(input_size_bytes, 1.0))))


def completed_history(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    result = []
    for row in rows:
        if row.get("status", "").strip().lower() not in SUCCESS_STATES:
            continue
        if as_float(row, "max_rss_bytes") is None or as_float(row, "runtime_sec") is None:
            continue
        result.append(row)
    return result


def context_is_known(candidate: dict[str, str], history: list[dict[str, str]]) -> bool:
    for column in CONTEXT_COLUMNS:
        value = candidate.get(column, "").strip()
        known = {row.get(column, "").strip() for row in history if row.get(column, "").strip()}
        if value and known and value not in known:
            return False
    return True


def static_recommendation(candidate: dict[str, str], reason: str) -> dict[str, object]:
    return {
        "task_key": candidate["task_key"],
        "run_id": candidate["run_id"],
        "rule": candidate["rule"],
        "input_size_bytes": candidate["input_size_bytes"],
        "suggested_memory_gb": candidate["static_memory_gb"],
        "suggested_runtime_min": candidate["static_runtime_min"],
        "sample_count": 0,
        "stratum": "static",
        "reason": reason,
        "mode": "shadow",
    }


def recommend_one(
    candidate: dict[str, str],
    history: list[dict[str, str]],
    quantile_level: float,
    memory_margin: float,
    runtime_margin: float,
    min_rule_samples: int,
    min_bucket_samples: int,
) -> dict[str, object]:
    rule_history = [row for row in history if row.get("rule") == candidate["rule"]]
    if len(rule_history) < min_rule_samples:
        return static_recommendation(candidate, "static_fallback_insufficient_rule_history")
    if not context_is_known(candidate, rule_history):
        return static_recommendation(candidate, "static_fallback_unknown_context")

    input_size = as_float(candidate, "input_size_bytes")
    assert input_size is not None
    sizes = [as_float(row, "input_size_bytes") for row in rule_history]
    if any(size is None for size in sizes):
        return static_recommendation(candidate, "static_fallback_invalid_history")
    numeric_sizes = [size for size in sizes if size is not None]
    if input_size < min(numeric_sizes) or input_size > max(numeric_sizes):
        return static_recommendation(candidate, "static_fallback_ood_input_size")

    bucket = size_bucket(input_size)
    bucket_history = [
        row for row in rule_history
        if size_bucket(as_float(row, "input_size_bytes") or 0) == bucket
    ]
    if len(bucket_history) >= min_bucket_samples:
        selected = bucket_history
        stratum = f"rule_size_2pow_{bucket}"
    else:
        selected = rule_history
        stratum = "rule"

    mem_gb = [as_float(row, "max_rss_bytes", 0.0) / GIB for row in selected]
    runtime_min = [as_float(row, "runtime_sec", 0.0) / 60 for row in selected]
    predicted_memory = quantile([value for value in mem_gb if value is not None], quantile_level) * memory_margin
    predicted_runtime = quantile([value for value in runtime_min if value is not None], quantile_level) * runtime_margin

    static_memory = as_float(candidate, "static_memory_gb", 0.0) or 0.0
    static_runtime = as_float(candidate, "static_runtime_min", 0.0) or 0.0
    memory_floor = as_float(candidate, "memory_floor_gb", static_memory) or static_memory
    runtime_floor = as_float(candidate, "runtime_floor_min", static_runtime) or static_runtime
    memory_ceiling = as_float(candidate, "memory_ceiling_gb", float("inf")) or float("inf")
    runtime_ceiling = as_float(candidate, "runtime_ceiling_min", float("inf")) or float("inf")
    suggested_memory = min(memory_ceiling, max(memory_floor, predicted_memory))
    suggested_runtime = min(runtime_ceiling, max(runtime_floor, predicted_runtime))

    return {
        "task_key": candidate["task_key"],
        "run_id": candidate["run_id"],
        "rule": candidate["rule"],
        "input_size_bytes": candidate["input_size_bytes"],
        "suggested_memory_gb": f"{suggested_memory:.3f}",
        "suggested_runtime_min": f"{suggested_runtime:.3f}",
        "sample_count": len(selected),
        "stratum": stratum,
        "reason": "quantile_upper_bound",
        "mode": "shadow",
    }


def shadow(
    history_path: str | Path,
    candidates_path: str | Path,
    output_path: str | Path,
    quantile_level: float,
    memory_margin: float,
    runtime_margin: float,
    min_rule_samples: int,
    min_bucket_samples: int,
) -> None:
    history = completed_history(read_tsv(history_path))
    candidates = read_tsv(candidates_path)
    validate_manifest(candidates)
    rows = [
        recommend_one(
            candidate,
            history,
            quantile_level,
            memory_margin,
            runtime_margin,
            min_rule_samples,
            min_bucket_samples,
        )
        for candidate in candidates
    ]
    fields = [
        "task_key", "run_id", "rule", "input_size_bytes", "suggested_memory_gb",
        "suggested_runtime_min", "sample_count", "stratum", "reason", "mode",
    ]
    write_tsv(output_path, rows, fields)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    collect = subparsers.add_parser("collect", help="join completed benchmarks to a task manifest")
    collect.add_argument("--manifest", required=True)
    collect.add_argument("--output", required=True)
    collect.add_argument("--max-rss-unit", choices=("bytes", "kb", "mb", "gb"), default="bytes")

    candidate = subparsers.add_parser("prepare-candidate", help="derive pre-submit size features from input metadata")
    candidate.add_argument("--task-key", required=True)
    candidate.add_argument("--run-id", required=True)
    candidate.add_argument("--rule", required=True)
    candidate.add_argument("--input", action="append", required=True)
    candidate.add_argument("--output", required=True)
    candidate.add_argument("--static-memory-gb", required=True, type=float)
    candidate.add_argument("--static-runtime-min", required=True, type=float)
    candidate.add_argument("--memory-floor-gb", type=float)
    candidate.add_argument("--runtime-floor-min", type=float)
    candidate.add_argument("--workflow-version", default="")
    candidate.add_argument("--tool-version", default="")
    candidate.add_argument("--reference-id", default="")
    candidate.add_argument("--node-class", default="")
    candidate.add_argument("--input-records", type=int)
    candidate.add_argument("--umi-groups", type=int)
    candidate.add_argument("--qualifying-umi", type=int)

    shadow_parser = subparsers.add_parser("shadow", help="write non-binding conservative recommendations")
    shadow_parser.add_argument("--history", required=True)
    shadow_parser.add_argument("--candidates", required=True)
    shadow_parser.add_argument("--output", required=True)
    shadow_parser.add_argument("--quantile", type=float, default=0.95)
    shadow_parser.add_argument("--memory-margin", type=float, default=1.30)
    shadow_parser.add_argument("--runtime-margin", type=float, default=1.25)
    shadow_parser.add_argument("--min-rule-samples", type=int, default=20)
    shadow_parser.add_argument("--min-bucket-samples", type=int, default=10)
    return parser


def main() -> None:
    args = build_parser().parse_args()
    if args.command == "collect":
        collect_history(args.manifest, args.output, args.max_rss_unit)
    elif args.command == "prepare-candidate":
        prepare_candidate(
            args.task_key,
            args.run_id,
            args.rule,
            args.input,
            args.output,
            args.static_memory_gb,
            args.static_runtime_min,
            args.memory_floor_gb,
            args.runtime_floor_min,
            args.workflow_version,
            args.tool_version,
            args.reference_id,
            args.node_class,
            args.input_records,
            args.umi_groups,
            args.qualifying_umi,
        )
    else:
        if not 0 < args.quantile <= 1:
            raise ValueError("--quantile must be in (0, 1]")
        if args.min_rule_samples < 1 or args.min_bucket_samples < 1:
            raise ValueError("minimum sample counts must be positive")
        shadow(
            args.history,
            args.candidates,
            args.output,
            args.quantile,
            args.memory_margin,
            args.runtime_margin,
            args.min_rule_samples,
            args.min_bucket_samples,
        )


if __name__ == "__main__":
    main()
