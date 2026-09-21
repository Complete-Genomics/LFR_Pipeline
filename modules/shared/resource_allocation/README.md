# cLFR conservative resource shadow

This module is derived from the safety contract of
`production_agenticML/subprojects/Adaptive_Resource_Allocation`, but stays
dependency-free and Snakemake-native. It does not train a point-prediction
model and never changes a workflow resource request.

## Data contract

Create one manifest row for each completed Snakemake job. Required columns are:

```text
task_key  run_id  rule  input_size_bytes  static_memory_gb  static_runtime_min  benchmark_path
```

`task_key` must include `run_id + rule + wildcard/shard + attempt`. Add
`workflow_version`, `tool_version`, `reference_id`, `node_class`, and workload
features such as `input_records`, `umi_groups`, and `qualifying_umi` when
available before submission. Completed rows may include `status=success`.

Build a history table after a run:

```bash
python modules/shared/resource_allocation/resource_allocation.py collect \
  --manifest completed_manifest.tsv \
  --output resource_allocation/history.tsv
```

The default assumes Snakemake's `max_rss` is bytes. Pass `--max-rss-unit kb`
only when the accounting export is known to use KiB.

Create a pre-submit candidate manifest with the same required fields, plus
optional `memory_floor_gb`, `runtime_floor_min`, `memory_ceiling_gb`, and
`runtime_ceiling_min`. An omitted floor defaults to the static request, so the
policy cannot reduce resources accidentally.

For a single candidate, derive file-size features without reading file content:

```bash
python modules/shared/resource_allocation/resource_allocation.py prepare-candidate \
  --task-key run42:map_reads_minimap:attempt1 \
  --run-id run42 --rule map_reads_minimap \
  --input data/split_read_2_trimmed.fastq.gz \
  --static-memory-gb 64 --static-runtime-min 720 \
  --workflow-version <git-commit> --node-class c40 \
  --output resource_allocation/candidates/map.tsv
```

Concatenate one-row candidate manifests only when their `task_key` values are
unique. Add already-computed UMI fields for consensus/denovo shards; never make
the pre-submit command scan the full input to obtain them.

## Policy and integration

For a candidate with sufficient same-rule history, the policy selects a same
`input_size_bytes` power-of-two bucket when it has enough rows; otherwise it
uses the whole rule. It applies a high quantile plus margins to RSS and runtime.
Unknown rule/context, input size outside historical range, or insufficient data
emits the original static request with a fallback reason.

Enable only after history is curated:

```yaml
modules:
  resource_allocation: true
resource_allocation:
  history: resource_allocation/history.tsv
  candidates: resource_allocation/candidates.tsv
```

This adds `resource_allocation/shadow.tsv` to the DAG as an audit artifact only.
It does not feed `threads`, `resources`, or shell commands. Evaluate grouped,
time-held-out coverage and scheduler-confirmed OOM/timeout outcomes before any
canary.
