"""Record pipeline code and configuration provenance for each run."""

import datetime
import hashlib
import json
import os
import shlex
import subprocess
import sys
from pathlib import Path


def _git(pipeline_root, *args):
    try:
        result = subprocess.run(
            ["git", "-C", str(pipeline_root)] + list(args),
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True,
            check=False,
        )
    except OSError:
        return ""
    if result.returncode != 0:
        return ""
    return result.stdout.rstrip("\n")


def record_run_metadata(pipeline_root, effective_config):
    """Write one immutable metadata directory for the current Snakemake run."""
    pipeline_root = Path(pipeline_root).resolve()
    timestamp = datetime.datetime.utcnow().strftime("%Y%m%dT%H%M%SZ")
    run_id = "{}_{}".format(timestamp, os.getpid())
    status = _git(pipeline_root, "status", "--porcelain=v1", "--untracked-files=all")
    diff = _git(pipeline_root, "diff", "--binary", "HEAD", "--")

    run_dir = Path.cwd() / "run_metadata" / run_id
    run_dir.mkdir(parents=True, exist_ok=False)
    config_text = json.dumps(
        effective_config,
        indent=2,
        sort_keys=True,
        default=str,
    ) + "\n"
    (run_dir / "effective_config.json").write_text(config_text)
    (run_dir / "git_status.txt").write_text(status + ("\n" if status else ""))
    (run_dir / "git_diff.patch").write_text(diff + ("\n" if diff else ""))

    branch = _git(pipeline_root, "symbolic-ref", "--short", "-q", "HEAD")
    metadata = {
        "run_id": run_id,
        "start_time_utc": timestamp,
        "working_directory": str(Path.cwd()),
        "pipeline_root": str(pipeline_root),
        "command": " ".join(shlex.quote(arg) for arg in sys.argv),
        "git_commit": _git(pipeline_root, "rev-parse", "HEAD") or "unknown",
        "git_branch": branch or "detached",
        "git_describe": _git(pipeline_root, "describe", "--always", "--tags", "--dirty") or "unknown",
        "git_dirty": "true" if status else "false",
        "config_sha256": hashlib.sha256(config_text.encode("utf-8")).hexdigest(),
    }
    with (run_dir / "metadata.tsv").open("w") as handle:
        for key, value in metadata.items():
            handle.write("{}\t{}\n".format(key, value))

    print("Run metadata: {}".format(run_dir), file=sys.stderr)
    if status:
        print(
            "WARNING: pipeline checkout is dirty; tracked changes are saved in git_diff.patch.",
            file=sys.stderr,
        )
    return run_dir


def record_run_metadata_safely(pipeline_root, effective_config):
    """Keep provenance diagnostics from blocking the scientific workflow."""
    try:
        return record_run_metadata(pipeline_root, effective_config)
    except Exception as exc:
        print("WARNING: could not record run metadata: {}".format(exc), file=sys.stderr)
        return None
