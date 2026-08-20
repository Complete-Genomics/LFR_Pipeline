# Unattended model retraining on a customer's own mock-community control.
#
# WHY THIS EXISTS. The shipped model (modules/clfr/denovo/models/) was trained
# on our ZymoBIOMICS control, on our library prep, our operator, our reagent
# lots and our instrument. A customer running the same protocol still differs
# in all of those, and what this model learns is precisely
# library/sequencing-artifact structure (read length, quality shape,
# homopolymers, within-pool k-mer popularity) rather than anything about which
# organisms are present -- so it is exactly the kind of model that can drift
# between labs. Retraining on the customer's own control is therefore worth
# doing, and since customers cannot be expected to make an ML judgement call,
# it has to happen without one.
#
# THREE STAGES, and it matters which is which -- they answer different
# questions and cannot substitute for each other:
#
#   ADMISSION -- CONTROL QC (controlQC_retrain). Is this control healthy
#      enough to be training data at all? denovo.md sec 61 found a real
#      production sample at 52.5% 16S content in a batch whose other samples
#      were fine; training on something like that would teach the model that
#      noise is normal. Failing aborts the retrain, keeps the incumbent, and
#      names the failing check -- an actionable message about the wet lab,
#      not an ML question.
#
#   TRIGGER -- DRIFT (driftReport_retrain, target auto_retrain/drift_report.tsv).
#      Has the feature distribution moved away from what the incumbent was
#      trained on, and if so, WHICH WAY? Three verdicts: `no_drift`,
#      `drift` (moved, quality intact -- a protocol change, retrain), and
#      `degradation` (moved, quality features moved the wrong way -- do NOT
#      retrain; the model would learn to treat the decay as normal and the
#      wet-lab problem would disappear into the product). PSI alone is
#      symmetric and cannot separate the last two, which is why
#      denovo_feature_drift.py tracks direction as well as distance.
#      Run this cheap target first (sort + features + PSI only, no assembly)
#      and only invoke the full retrain when the verdict is `drift`.
#
#      What one control CANNOT tell you is whether a degradation is a bad run
#      or a permanently changed process -- that needs several consecutive
#      controls, i.e. history this workflow does not keep. The directional
#      check is the single-sample approximation, and it errs toward refusing
#      to retrain.
#
#   GATE -- PROMOTION (promotionGate_retrain). Is the candidate actually
#      BETTER? Drift says the input moved, not that the new model wins, and
#      training metrics do not decide either: sec 50 is the standing
#      counterexample of a model with good MAE that made a worse filter. So
#      both models actually filter and assemble the control, the contigs are
#      scored against the control's KNOWN reference, and the candidate must
#      beat the incumbent by a margin without an excessive severe-loss tail or
#      a loss of assembled/read-backed yield. This is a true A/B on real ground
#      truth, unlike random_inspection.smk's greengenes proxy for field
#      samples (which compresses effect size 3-5x and cannot rank close arms
#      -- sec 55/59 -- so it is the wrong instrument for an automatic
#      promote/reject).
#
# Training data is ALWAYS the control, never field samples: an unlabeled
# production sample offers no way to check whether the resulting model got
# better or worse.
#
# Output contract: models/in_use.lgb is what denovo_olc.smk should be pointed
# at (frag_de_novo.read_filter_ml_model). It is written on every run, whether
# the candidate was promoted or rejected, so the production config never has
# to change. retrain_decision.tsv records which model won and why; the
# candidate is kept on disk either way, so a rejected model can be inspected
# and a promoted one rolled back by hand.
#
# Usage -- cheap trigger check first, full retrain only if it says drift.
# Note the exact match on `drift`: a `degradation` verdict must not start a
# retrain, so a substring test would be wrong here.
#   snakemake -s auto_retrain.smk --configfile <config.yaml> -j <N> \
#       auto_retrain/drift_report.tsv
#   awk -F'\t' '$1=="#verdict"{exit $2!="drift"}' auto_retrain/drift_report.tsv && \
#   snakemake -s auto_retrain.smk --configfile <config.yaml> -j <N> \
#       auto_retrain/retrain_decision.tsv
#
# Config (under `auto_retrain:`):
#   control_fastq     barcode-tagged FASTQ(.gz) of the customer's mock control
#   control_reference FASTA of that mock's known 16S sequences
#   outdir            default "auto_retrain"
#   incumbent_model   default modules/clfr/denovo/models/model_identity.lgb
#   label_id          vsearch identity floor for LABELING, default 0.50.
#                     NOT a quality threshold -- it is deliberately permissive
#                     so noisy reads still receive a (low) identity label,
#                     which is the signal the model is trained to predict.
#                     Raising it silently truncates the training distribution:
#                     0.70 dropped label coverage 94.4% -> 86.8% and pushed
#                     mean label identity 91.2 -> 94.0 (denovo.md sec 75).
#   min_16s_rate      control QC: min fraction of contigs matching the
#                     reference, default 0.75
#   min_barcodes      control QC: min assembled barcodes, default 1000
#   baseline_profile  drift: feature profile the incumbent was trained on,
#                     default <incumbent_model minus .lgb>.baseline.json
#   psi_threshold     drift: per-feature PSI counted as shifted, default 0.25
#   min_flagged       drift: shifted features constituting drift, default 2
#   degradation_tolerance
#                     drift: relative worsening of a shifted quality feature
#                     that turns the verdict into `degradation`, default 0.05
#   min_improvement   promotion margin in identity points, default 0.05
#   severe_loss_points promotion: paired identity drop counted as severe,
#                     default 5.0
#   max_severe_loss_rate
#                     promotion: maximum fraction of severe regressions,
#                     default 0.01
#   min_primary_len    promotion: length defining assembled-contig yield,
#                     default 1000
#   min_yield_ratio, min_readback_bp_ratio
#                     promotion: candidate must not reduce >=min_primary_len
#                     yield or raw-read-supported bases; both default 1.0
#   vsearch, samtools, num_processes

import os

AR = config.get('auto_retrain', {})
AR_OUT = AR.get('outdir', 'auto_retrain')
AR_SRC = config['params']['src_dir'] + "/modules/clfr/denovo"
AR_PYTHON = config['params']['general_python']
AR_INCUMBENT = AR.get('incumbent_model', config['params']['src_dir'] +
                       "/modules/clfr/denovo/models/model_identity.lgb")
AR_CONTROL_FQ = AR.get('control_fastq')
AR_CONTROL_REF = AR.get('control_reference')
AR_VSEARCH = AR.get('vsearch', 'vsearch')
AR_SAMTOOLS = AR.get('samtools', 'samtools')
AR_NPROC = AR.get('num_processes', 4)
AR_LABEL_ID = AR.get('label_id', 0.50)
AR_MIN_16S = AR.get('min_16s_rate', 0.75)
AR_MIN_BC = AR.get('min_barcodes', 1000)
AR_MIN_IMPROVEMENT = AR.get('min_improvement', 0.05)
AR_SEVERE_LOSS_POINTS = AR.get('severe_loss_points', 5.0)
AR_MAX_SEVERE_LOSS_RATE = AR.get('max_severe_loss_rate', 0.01)
AR_MIN_PRIMARY_LEN = AR.get('min_primary_len', 1000)
AR_MIN_YIELD_RATIO = AR.get('min_yield_ratio', 1.0)
AR_MIN_READBACK_BP_RATIO = AR.get('min_readback_bp_ratio', 1.0)
AR_BASELINE = AR.get('baseline_profile',
                      AR_INCUMBENT[:-len(".lgb")] + ".baseline.json"
                      if AR_INCUMBENT.endswith(".lgb") else AR_INCUMBENT + ".baseline.json")
AR_PSI = AR.get('psi_threshold', 0.25)
AR_MIN_FLAGGED = AR.get('min_flagged', 2)
AR_DEG_TOL = AR.get('degradation_tolerance', 0.05)

wildcard_constraints:
    model = "incumbent|candidate"

# See random_inspection.smk for why .strip() is required here.
_PRIMARY_AWK = r"""
    awk '
        /^>/ {{
            if (seq != "" && keep) print header"\n"seq
            header=$0
            keep = (header ~ /k41_0$/)
            seq=""
            next
        }}
        {{ seq = seq $0 }}
        END {{ if (seq != "" && keep) print header"\n"seq }}
    '
""".strip()


rule auto_retrain_all:
    input:
        f"{AR_OUT}/retrain_decision.tsv"


rule sortControl_retrain:
    input:
        AR_CONTROL_FQ
    output:
        f"{AR_OUT}/control_sorted.fastq.gz"
    benchmark:
        "Benchmarks/auto_retrain.sortControl.txt"
    shell:
        "bash {AR_SRC}/denovo_sort_fastq.sh {input} {output} {AR_NPROC} {AR_SAMTOOLS}"


rule featuresControl_retrain:
    input:
        f"{AR_OUT}/control_sorted.fastq.gz"
    output:
        f"{AR_OUT}/control_features.tsv"
    benchmark:
        "Benchmarks/auto_retrain.featuresControl.txt"
    shell:
        "{AR_PYTHON} {AR_SRC}/denovo_read_features.py "
        "--fastq {input} --out {output} --num_processes {AR_NPROC}"


## TRIGGER. Deliberately cheap -- depends only on sort + features, no
## assembly and no vsearch, so a scheduler can poll it often. Measured on
## real data (denovo.md sec 78): a soil field sample against the shipped
## Zymo-trained baseline scores max PSI 0.0775 -> no_drift, confirming the
## detector does not fire merely because the sample type changed; the same
## features with quality degradation injected score PSI 5-11 on exactly the
## degraded columns -> degradation; a simulated chemistry change (longer
## reads, deeper pools, quality intact) -> drift.
rule driftReport_retrain:
    input:
        features = f"{AR_OUT}/control_features.tsv",
        baseline = AR_BASELINE
    output:
        f"{AR_OUT}/drift_report.tsv"
    benchmark:
        "Benchmarks/auto_retrain.driftReport.txt"
    shell:
        "{AR_PYTHON} {AR_SRC}/denovo_feature_drift.py "
        "--baseline {input.baseline} --features {input.features} --out {output} "
        "--psi-threshold {AR_PSI} --min-flagged {AR_MIN_FLAGGED} "
        "--degradation-tolerance {AR_DEG_TOL}"


## Refuses to proceed into training when the trigger said `degradation`.
## Belt and braces with the documented two-step invocation above: if an
## operator (or a cron job with a sloppy grep) runs the retrain target
## directly, this still stops it, because retraining on decaying data is the
## one outcome that silently makes the product worse.
rule driftGate_retrain:
    input:
        f"{AR_OUT}/drift_report.tsv"
    output:
        f"{AR_OUT}/drift_gate_ok"
    run:
        verdict = None
        degraded = ""
        for line in open(input[0]):
            if line.startswith("#verdict\t"):
                verdict = line.rstrip("\n").split("\t")[1]
            elif line.startswith("#degraded_features\t"):
                degraded = line.rstrip("\n").split("\t")[1]
        if verdict == "degradation":
            raise ValueError(
                "control shows DEGRADATION, not drift (features moved the wrong "
                "way: {}). Refusing to retrain -- a model trained on this would "
                "learn to treat the decay as normal. Investigate the library "
                "prep/sequencing run; see {} for per-feature detail."
                .format(degraded, input[0]))
        with open(output[0], "w") as fh:
            fh.write("verdict\t{}\n".format(verdict))


rule sgrepControl_retrain:
    input:
        f"{AR_OUT}/control_sorted.fastq.gz"
    output:
        f"{AR_OUT}/control_nofilter.tsv"
    shell:
        """
        gzip -dc {input} | \
        awk '{{if (NR%4==1) {{temp=$1; $1=$2; $2=temp}} }}1' | \
        awk '{{if (NR%4==1) line=line$0"\\t"; if (NR%4==2) {{print line$0; line=""}}}}' \
        > {output}
        """


## Unfiltered assembly of the control, used for QC only -- assembling before
## any filtering is what makes the QC verdict independent of the models being
## compared downstream.
rule assembleControlRaw_retrain:
    input:
        f"{AR_OUT}/control_nofilter.tsv"
    output:
        f"{AR_OUT}/qc/denovo/final_contigs_0.fa"
    benchmark:
        "Benchmarks/auto_retrain.assembleControlRaw.txt"
    params:
        d = f"{AR_OUT}/qc"
    shell:
        "mkdir -p {params.d} && cd {params.d} && "
        "{AR_PYTHON} {AR_SRC}/denovo_seed_olc.py --sequence_type se "
        "--num_processes {AR_NPROC} --r2 $OLDPWD/{input}"


rule qcPrimary_retrain:
    input:
        f"{AR_OUT}/qc/denovo/final_contigs_0.fa"
    output:
        f"{AR_OUT}/qc/primary.fa"
    shell:
        _PRIMARY_AWK + " {input} > {output}"


rule qcHits_retrain:
    input:
        query = f"{AR_OUT}/qc/primary.fa",
        ref = AR_CONTROL_REF
    output:
        f"{AR_OUT}/qc/hits.tsv"
    shell:
        "{AR_VSEARCH} --usearch_global {input.query} --db {input.ref} "
        "--id 0.75 --strand both --maxaccepts 1 --top_hits_only "
        "--userout {output} --userfields query+target+id "
        "--threads {AR_NPROC} --log {output}.log"


## Gate 1. Fails loudly (nonzero exit) so the workflow stops before training
## on a control that should not be trained on; the message names the failing
## check so the reply to the customer is about their control, not about ML.
rule controlQC_retrain:
    input:
        contigs = f"{AR_OUT}/qc/primary.fa",
        hits = f"{AR_OUT}/qc/hits.tsv"
    output:
        f"{AR_OUT}/control_qc.tsv"
    run:
        n_contigs = sum(1 for line in open(input.contigs) if line.startswith(">"))
        matched = len({line.split("\t")[0] for line in open(input.hits) if "\t" in line})
        rate = matched / n_contigs if n_contigs else 0.0
        with open(output[0], "w") as fh:
            fh.write("check\tvalue\tthreshold\tpass\n")
            fh.write("assembled_barcodes\t{}\t{}\t{}\n".format(
                n_contigs, AR_MIN_BC, n_contigs >= AR_MIN_BC))
            fh.write("reference_match_rate\t{:.4f}\t{}\t{}\n".format(
                rate, AR_MIN_16S, rate >= AR_MIN_16S))
        problems = []
        if n_contigs < AR_MIN_BC:
            problems.append(
                "only {} barcodes assembled (need >= {})".format(n_contigs, AR_MIN_BC))
        if rate < AR_MIN_16S:
            problems.append(
                "only {:.1%} of contigs match the control reference (need >= {:.0%}) "
                "-- the control looks contaminated or the library underperformed"
                .format(rate, AR_MIN_16S))
        if problems:
            raise ValueError(
                "control failed QC, refusing to retrain on it; keeping the "
                "incumbent model. " + "; ".join(problems))


## Labels: per-read identity against the control's known reference. --id is
## permissive on purpose (see label_id in the config notes above).
rule labelControl_retrain:
    input:
        qc = f"{AR_OUT}/control_qc.tsv",
        reads = f"{AR_OUT}/control_nofilter.tsv",
        ref = AR_CONTROL_REF
    output:
        fasta = temp(f"{AR_OUT}/control_reads.fa"),
        hits = f"{AR_OUT}/control_labels.tsv"
    benchmark:
        "Benchmarks/auto_retrain.labelControl.txt"
    shell:
        """
        awk -F'\\t' '{{print ">" substr($1, 23) "\\n" $2}}' {input.reads} > {output.fasta}
        {AR_VSEARCH} --usearch_global {output.fasta} --db {input.ref} \
            --id {AR_LABEL_ID} --strand both --maxaccepts 4 \
            --userout {output.hits} --userfields query+target+id \
            --threads {AR_NPROC} --log {output.hits}.log
        """


## Depends on drift_gate_ok, not merely on drift_report: reaching training at
## all requires the trigger to have cleared, so going straight to the retrain
## target cannot bypass the degradation check.
rule trainCandidate_retrain:
    input:
        features = f"{AR_OUT}/control_features.tsv",
        labels = f"{AR_OUT}/control_labels.tsv",
        drift_ok = f"{AR_OUT}/drift_gate_ok"
    output:
        model = f"{AR_OUT}/candidate.lgb",
        metrics = f"{AR_OUT}/candidate_metrics.tsv"
    benchmark:
        "Benchmarks/auto_retrain.trainCandidate.txt"
    shell:
        "{AR_PYTHON} {AR_SRC}/denovo_train_model.py "
        "--features {input.features} --labels-glob '{input.labels}' "
        "--out-model {output.model} --out-metrics {output.metrics}"


def _model_for(wildcards):
    return (AR_INCUMBENT if wildcards.model == "incumbent"
            else f"{AR_OUT}/candidate.lgb")


## Both arms filter the SAME control reads, differing only in which model
## breaks conflict-graph ties -- the comparison the promotion gate needs.
rule filterArm_retrain:
    input:
        reads = f"{AR_OUT}/control_nofilter.tsv",
        features = f"{AR_OUT}/control_features.tsv",
        model = _model_for
    output:
        f"{AR_OUT}/{{model}}_filtered.tsv"
    benchmark:
        "Benchmarks/auto_retrain.filterArm.{model}.txt"
    shell:
        "{AR_PYTHON} {AR_SRC}/denovo_read_filter.py --r2 {input.reads} "
        "--out {output} --num_processes {AR_NPROC} "
        "--ml-model {input.model} --ml-features {input.features}"


rule assembleArm_retrain:
    input:
        f"{AR_OUT}/{{model}}_filtered.tsv"
    output:
        f"{AR_OUT}/{{model}}/denovo/final_contigs_0.fa"
    benchmark:
        "Benchmarks/auto_retrain.assembleArm.{model}.txt"
    params:
        d = lambda wc: f"{AR_OUT}/{wc.model}"
    shell:
        "mkdir -p {params.d} && cd {params.d} && "
        "{AR_PYTHON} {AR_SRC}/denovo_seed_olc.py --sequence_type se "
        "--num_processes {AR_NPROC} --r2 $OLDPWD/{input}"


rule armPrimary_retrain:
    input:
        f"{AR_OUT}/{{model}}/denovo/final_contigs_0.fa"
    output:
        f"{AR_OUT}/{{model}}/primary.fa"
    shell:
        _PRIMARY_AWK + " {input} > {output}"


## The two arms are read back against the SAME unfiltered control reads.
## Using each arm's filtered read set here would let a model hide a bad
## extension by deleting the reads that disagree with it.
rule armReadback_retrain:
    input:
        reads = f"{AR_OUT}/control_nofilter.tsv",
        contigs = f"{AR_OUT}/{{model}}/primary.fa"
    output:
        f"{AR_OUT}/{{model}}_readback.tsv"
    benchmark:
        "Benchmarks/auto_retrain.armReadback.{model}.txt"
    shell:
        "{AR_PYTHON} {AR_SRC}/test/readback_qc_single.py "
        "--r2 {input.reads} --contigs {input.contigs} --num_processes {AR_NPROC} "
        "--out {output}"


## Scored against the control's real reference -- the whole point of retraining
## on a mock rather than a field sample.
rule armHits_retrain:
    input:
        query = f"{AR_OUT}/{{model}}/primary.fa",
        ref = AR_CONTROL_REF
    output:
        f"{AR_OUT}/{{model}}_hits.tsv"
    shell:
        "{AR_VSEARCH} --usearch_global {input.query} --db {input.ref} "
        "--id 0.5 --strand both --maxaccepts 0 --maxrejects 0 "
        "--userout {output} --userfields query+target+id "
        "--threads {AR_NPROC} --log {output}.log"


## Gate 2.
rule promotionGate_retrain:
    input:
        incumbent = f"{AR_OUT}/incumbent_hits.tsv",
        candidate = f"{AR_OUT}/candidate_hits.tsv",
        incumbent_contigs = f"{AR_OUT}/incumbent/primary.fa",
        candidate_contigs = f"{AR_OUT}/candidate/primary.fa",
        incumbent_readback = f"{AR_OUT}/incumbent_readback.tsv",
        candidate_readback = f"{AR_OUT}/candidate_readback.tsv",
        candidate_model = f"{AR_OUT}/candidate.lgb"
    output:
        decision = f"{AR_OUT}/retrain_decision.tsv",
        in_use = f"{AR_OUT}/models/in_use.lgb"
    shell:
        "mkdir -p $(dirname {output.in_use}) && "
        "{AR_PYTHON} {AR_SRC}/denovo_promotion_gate.py "
        "--incumbent-hits {input.incumbent} --candidate-hits {input.candidate} "
        "--incumbent-contigs {input.incumbent_contigs} --candidate-contigs {input.candidate_contigs} "
        "--incumbent-readback {input.incumbent_readback} --candidate-readback {input.candidate_readback} "
        "--incumbent-model {AR_INCUMBENT} --candidate-model {input.candidate_model} "
        "--min-improvement {AR_MIN_IMPROVEMENT} --severe-loss-points {AR_SEVERE_LOSS_POINTS} "
        "--max-severe-loss-rate {AR_MAX_SEVERE_LOSS_RATE} --min-primary-len {AR_MIN_PRIMARY_LEN} "
        "--min-yield-ratio {AR_MIN_YIELD_RATIO} --min-readback-bp-ratio {AR_MIN_READBACK_BP_RATIO} "
        "--out-decision {output.decision} --out-promoted-model {output.in_use}"
