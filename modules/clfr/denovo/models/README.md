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

## model_candidate_chimera.lgb

LightGBM classifier predicting P(confident_chimera) for one assembled
candidate contig (a UMI's k41_0..k41_N OLC output), from 7 features computed
by `denovo_junction_qc.py --all-candidates`. This is a different model from
`model_identity.lgb` above (that one scores individual reads pre-assembly;
this one scores whole assembled candidate contigs post-assembly) and is
consumed by `denovo_shadow_score.py`, wired through the optional
`denovo_shadow.smk` sidecar via `frag_de_novo.shadow_score_model`. It writes
`denovo/shadow/candidate_scores.tsv` and `denovo/shadow/umi_summary.tsv` only
when those shadow targets are explicitly requested; neither output is an
input to the production delivery target.

**Production role: SHADOW ONLY.** This model's score must never be allowed
to change which candidate gets delivered. `denovo_candidate_select.py`'s
`gated_switch` mode (the validated, rule-only candidate switch) does not use
it at all. Every active-selection deployment shape tried for this model
(argmin-in-gate ranking, length/support guardrails on top of it, an
asymmetric veto) either reproduced or failed to fix a severe-loss regression
too large to ship, and all training/eval was on the ZymoBIOMICS control, so
there is currently no way to verify it against ground truth on a real
sample. See `denovo.md` sec 106-120 for the full narrative (sec 118 has the
deployment decision; sec 120 has the most recent revisit of it) before ever
reconsidering a more active role.

**Feature contract** (read by name via pandas/the training script, not
positional):

    span_cov_ratio, min_local_span_ratio, placed_reads, contig_len,
    len_ratio, k41_rank, n_candidates

`len_ratio` and `n_candidates` are derived per barcode from the full
candidate set (not columns `denovo_junction_qc.py` writes directly) — see
`denovo_shadow_score.py`'s `build_features()` for the exact derivation used
at inference time, and `denovo_train_candidate_model.py`'s
`add_derived_features()` for training time.

**Provenance**: trained 2026-08-20 on the ZymoBIOMICS mock community,
reproducing the exact frozen recipe validated in denovo.md sec 107-112
(`LGBMClassifier(n_estimators=300, max_depth=5, learning_rate=0.05,
num_leaves=31, min_child_samples=30)`, quarter-split ground-truth labels,
57,840 decidable candidates pooled from the original 2,000-barcode training
set plus the sec109 20,000-barcode expansion, 5-fold GroupKFold-by-barcode
CV accuracy 0.9742). Retraining script: `denovo_train_candidate_model.py`.
The training workspace (filtered candidate-qc/label TSVs) is scratch, not
git-tracked; the underlying labeled data lives in
`salvage/2026-08-19_ml_draft_gain_optionB/` and
`salvage/2026-08-19_n5k_20k_case_expansion/`.

**Validation status**: the model's *ranking* ability was validated on a
one-time held-out set (tail_raw, sec 112: within-UMI pairwise accuracy 0.85
vs 0.76 for the best single rule feature — genuine generalization, not
overfitting). Its *deployment* as an active decision-maker was NOT
validated — see the shadow-only note above. If this file is retrained, the
sec 112 tail_raw numbers no longer apply to it until re-verified (that
held-out set has already been spent once and should not be re-used
casually — see denovo.md sec 112's own note on this).
