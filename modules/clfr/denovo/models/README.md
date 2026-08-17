# models — shipped model artifacts

Canonical, git-tracked home for small (<a few MB) trained models this module
consumes at runtime. Large training workspaces (raw reads, label chunks,
retraining scripts) do NOT belong here — they stay in whatever scratch/working
area produced the model; this directory only holds the artifact itself plus
enough provenance to trace it back.

## Zymo noisy-preprocess reference distributions

`zymo_reads_per_umi_distribution.tsv` and
`zymo_read_length_distribution.tsv` are the two runtime references used by
`denovo_noisy_preprocess_qc.py`.  They were built from
`test/olc/ml_pre_filter/zymo_train_eval_5000bc.fastq.gz` (221,400 reads,
5,000 UMIs); each file records the source SHA-256 and fixed analysis bins.
Regenerate them with:

    python modules/clfr/denovo/denovo_noisy_preprocess_qc.py \
      --r2 test/olc/ml_pre_filter/zymo_train_eval_5000bc.fastq.gz \
      --r2-format fastq --build-baseline-dir modules/clfr/denovo/models

The automatic gate requires PSI >=0.25 for both distributions.  PSI alone is
not enough: the top-1%-UMI read share and the <300-bp read fraction must each
be at least 0.10 above Zymo, and projected post-filter depth must remain safe.
This directional guard keeps a uniformly deeper sample from being mistaken
for the hs1 failure mode.

## model_identity.lgb

LightGBM regressor predicting per-read identity from `denovo_read_features.py`'s
13-column feature set. Consumed by `denovo_read_filter.py`'s `--ml-model` flag
as a conflict-graph tie-break (see that file's module docstring and
`find_contaminants()`'s `ml_scores` parameter).

**Feature contract** (order matters — this is a positional contract, not named
columns):

    length, qual_mean, qual_min, qual_head, qual_tail, qual_trend,
    hp_frac, hp_max_run, mini_gap_mean, mini_gap_max, mini_gap_var,
    pool_size, pool_kmer_popular_frac

`denovo_read_features.py`'s `FIELDS` list is the authority for this; if it
changes, this model is invalidated until retrained.

**Provenance**: trained on the ZymoBIOMICS mock community (vsearch-labeled
per-read identity against the known reference), reproduced 2026-08-13 after
the original training artifacts were lost with a `/tmp` wipe. Full retraining
history, the reproduction match table, and the training workspace (raw reads,
feature/label TSVs, `train_model.py`, `run_chunked_label.sh`) live in
`subprojects/olc/mlpf/` (not git-tracked — that directory is scratch/working
space, not a release artifact). See `denovo.md` sec 74-76 for the full
narrative, including two real mistakes made and corrected while reproducing
it.

**Validation status**: the conflict-graph tie-break this model drives was
validated (+0.03..+0.07 identity vs the plain filter across 4 samples,
denovo.md sec 51/59/62/67) using the *original* model, which was lost before
this copy was reproduced. `random_inspection.smk` exists to re-run that same
validation on demand against any model file, including this one — check
denovo.md for whether that re-validation has been completed with this exact
file before trusting the original numbers apply here unchanged.
