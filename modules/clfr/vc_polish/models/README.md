# models — shipped model artifacts

Deployed copy of the trained Claim A model `vc_polish.smk` /
`04_apply_rescore.py` consume at runtime. Training happens entirely in the
separate `cLFR_eval` repo (`src/01_make_candidates.py` / `03_train_eval.py`,
GIAB-truth-labeled) — nothing here trains anything. This directory only
holds the deployable artifact (`.lgb.txt`) plus enough provenance to trace it
back to that source repo. Full diagnostics (feature_importance/calibration
CSVs, all 22 HG002 + 18 HG004 LOCO fold models, the raw features TSVs) live
in `cLFR_eval/models/` and `cLFR_eval/out/genomewide/` — not copied here to
keep this deployed module lean; see `cLFR_eval/models/README.md` for the
complete version of this document.

## claimA_snv_confidence.hg002_chr20holdout.2026-08-25.lgb.txt

Filename carries the training date (`2026-08-25`) as the version marker --
retraining the same fold (or a different held-out chrom) produces a
differently-dated file rather than overwriting this one silently.

LightGBM binary classifier: P(true variant) for a candidate SNP, from the 25
molecule-linkage + per-read + pileup features `02_extract_features.py`
produces (`ALL_FEATURES` in `03_train_eval.py` / `04_apply_rescore.py` — order
doesn't matter for LightGBM's `.txt` format, it stores feature names, but the
25-column contract itself does: retraining that changes `FEATURES` in
`02_extract_features.py` invalidates this model until retrained).

**Provenance**: one fold of the HG002 22-chromosome genome-wide
leave-one-chromosome-out CV run 2026-08-25 (`--feature-set all`, cLFR_eval's
`out/genomewide/loco_hg002/chr20_all/`). Trained on chr1-19,21,22 (21
chromosomes, 15,070,162 fit + 3,765,259 calibration rows), held out **chr20**
for validation. chr20 was chosen as the held-out fold specifically *because*
it's the chromosome most extensively characterized as "clean" across this
project's iteration (no repeat/segdup concentration like chr22, no depth
anomaly like chr19 — see cLFR_eval's `memory/plan.md`), so its held-out score is a
trustworthy, non-cherry-picked read of accuracy. This is **one of 22
statistically equivalent LOCO fold models** (each trained on ~95% of the
genome-wide data, differing only in which chromosome was excluded) — chr20
was picked for the *validation number attached to this shipped copy*, not
because this particular fold's weights are special.

**Validation status (held-out chr20, HG002)**:
- PR-AUC = 0.7954, ROC-AUC = 0.9969, Brier = 0.0015 (well-calibrated)
- UMI/molecule contribution (PR-AUC(all) − PR-AUC(no_molecule)) = **+0.0349**
- Full 22-fold genome-wide picture (this fold is representative, not an
  outlier): mean PR-AUC 0.7578 ± 0.0666, mean UMI delta +0.0258 ± 0.0223,
  2/22 folds negative (chr15, chr19 — both HG002-specific depth-covariate
  anomalies, diagnosed and NOT present on the same chromosomes in HG004; see
  cLFR_eval's `memory/plan.md`).
- Cross-sample replication (HG004, 18/22 chromosomes so far): mean PR-AUC
  0.7522 ± 0.0364, mean UMI delta +0.0236 ± 0.0102, **0/18 folds negative**.
- VAF-stratified: UMI contribution follows an inverted-U, peaking at mid VAF
  (0.05–0.5) at 3–4× the aggregate delta on both samples — see plan.md for the
  full per-bin table. Report VAF-stratified numbers, not one headline AUC.

**What this model is NOT validated for**:
- **Not titration/low-AF-sensitivity tested** (Claim B). This is a Claim A
  model: error suppression + calibrated confidence on germline-AF (~50/100%)
  HG002 truth. Do not claim low-AF (ctDNA-style) sensitivity from this model
  without a titration validation pass.
- **Not validated on genuine hom-ref-vs-error at production inference time**
  — training/validation candidates came from `01_make_candidates.py`
  (GIAB-confident-BED + truth-VCF labeled), not from a real called VCF's SNP
  list the way `vc_polish.smk` feeds it in production. The feature
  *extraction* code path is identical (`02_extract_features.py` doesn't care
  where candidate positions came from), but the candidate *composition*
  (a real caller's SNP calls vs. a labeled train/test split) hasn't been
  cross-checked for distribution shift.
- **DNA (HG002 WGS-style GIAB truth) trained, not RNA-recalibrated.** If
  applied to cLFR RNA/isoform consensus SNPs (Step 4's use case), recalibrate
  on ERCC (error-only truth) first — see cLFR_eval's `README.md` Step 4 caveats.
- **HG004 gap**: HG004 genome-wide LOCO is missing chr2/chr5/chr10/chr15
  (feature extraction not finished on the server) — HG004 replication is
  strong but not yet complete across all 22 chromosomes.

**Feature contract** (25 columns, `FEATURES` in `02_extract_features.py` /
`ALL_FEATURES` in `03_train_eval.py` and `04_apply_rescore.py` — read by name
via pandas, not positional):

    dp, alt_reads, ref_reads, vaf,
    n_mol_total, n_mol_alt, n_mol_ref, mol_alt_fraction,
    alt_reads_per_alt_mol_mean, within_mol_alt_agreement_mean,
    alt_bq_mean, alt_bq_min, ref_bq_mean,
    alt_mapq_mean, alt_mapq_min, ref_mapq_mean,
    alt_strand_balance, alt_softclip_frac_mean,
    alt_readpos_fromend_mean, alt_readpos_fromend_min,
    alt_indel_near_frac, alt_nm_mean, homopolymer_run,
    alt_clip_frac_mean, alt_supplementary_frac

**Do not use `--feature-set no_abs` as a fix for the chr15/chr19 anomaly.**
Tested genome-wide 2026-08-25: it fixes chr15 but breaks
chr1/chr7/chr12/chr22 and worsens chr19; net negative across 22 folds (mean
PR-AUC 0.7482 vs 0.7578, 5/22 negative folds vs 2/22). `all` (this model) is
the validated default. Full experiment and the "why the earlier chr22
finding didn't generalize" lesson: cLFR_eval's `memory/plan.md`.

This deployed copy intentionally does **not** include
`claimA_snv_confidence.hg002_chr20holdout.2026-08-25.feature_importance.csv`
or `.calibration.csv` -- the numbers are already in this README's prose above.
The raw CSVs live alongside the full-precision model copy in cLFR_eval's
`models/` directory if you need them.

**Regenerate / pick a different fold**: in the `cLFR_eval` repo,
`out/genomewide/loco_hg002/` and `out/genomewide/loco_hg004/` hold all 22
(HG002) / 18 (HG004) LOCO fold models — any `<chrom>_all/model.txt` is a
drop-in alternative (copy it here and update the `vc_polish:` config's
`model` path, or overwrite this file keeping the same name). Retrain from
scratch with `cLFR_eval/src/run.sh` (single train/test split) or the
`out/genomewide/run_loco_hg002.sh` pattern (full LOCO sweep).
