# Standalone gated-switch selection and high-confidence QC for an existing
# final_contigs_0.fa.  It is intentionally not included in the production
# DAG: callers provide the exact read pools belonging to that assembly.
#
# Usage:
#   snakemake -s standalone_gated_switch.smk --configfile <config.yaml> -j <N>
#
# Config (under `standalone_gated_switch:`):
#   contigs             required raw final_contigs_0.fa with all candidates
#   raw_r2              required data_R2_sorted.tsv (or noisy-preprocessed TSV)
#   filtered_r2         optional read-filtered TSV; defaults to raw_r2
#   outdir              default directory containing contigs
#   highconf            default <outdir>/final_contigs_highconf.fa
#   python              default current Python executable
#   vsearch             default "vsearch"
#   num_processes       default 4
#   max_span_ratio      default 0.25
#   min_local_span_ratio default 0.0 (disabled)
#   min_placed_reads    default 5 (sec 117/123/124)
#   uchime              default false; opt-in reference-free UCHIME QC

import sys
from pathlib import Path


SG = config.get("standalone_gated_switch", {})


def sg_required(name):
    value = SG.get(name)
    if not value:
        raise ValueError("standalone_gated_switch.{} is required".format(name))
    return value


SG_CONTIGS = sg_required("contigs")
SG_RAW_R2 = sg_required("raw_r2")
# A standalone rerun may predate the optional read-filter stage; in that case
# run read-back QC against the same sorted pool used by the candidate QC.
SG_FILTERED_R2 = SG.get("filtered_r2") or SG_RAW_R2
SG_OUT = SG.get("outdir", str(Path(SG_CONTIGS).parent))
SG_HIGHCONF = SG.get("highconf", f"{SG_OUT}/final_contigs_highconf.fa")
SG_PYTHON = SG.get("python", sys.executable)
SG_VSEARCH = SG.get("vsearch", "vsearch")
SG_NPROC = SG.get("num_processes", 4)
SG_MAX_SPAN = SG.get("max_span_ratio", 0.25)
SG_MIN_LOCAL_SPAN = SG.get("min_local_span_ratio", 0.0)
SG_MIN_PLACED = SG.get("min_placed_reads", 5)
SG_UCHIME = SG.get("uchime", False)
SG_SRC = str(Path(workflow.basedir).resolve())


rule standalone_gated_switch_all:
    input:
        SG_HIGHCONF,
        f"{SG_OUT}/qc_report.tsv",
        f"{SG_OUT}/candidate_select_report.tsv"


rule junctionQCAllCandidates_standalone:
    input:
        contigs=SG_CONTIGS,
        r2=SG_RAW_R2
    output:
        f"{SG_OUT}/junction_qc_candidates.tsv"
    shell:
        "mkdir -p {SG_OUT} && "
        "{SG_PYTHON} {SG_SRC}/denovo_junction_qc.py "
        "--r2 {input.r2} --contigs {input.contigs} --all-candidates "
        "--max-span-ratio {SG_MAX_SPAN} "
        "--min-local-span-ratio {SG_MIN_LOCAL_SPAN} "
        "--num_processes {SG_NPROC} --out {output}"


rule candidateSelect_standalone:
    input:
        contigs=SG_CONTIGS,
        candidate_qc=f"{SG_OUT}/junction_qc_candidates.tsv"
    output:
        fasta=f"{SG_OUT}/gated_switch.fasta",
        report=f"{SG_OUT}/candidate_select_report.tsv",
        decision=f"{SG_OUT}/candidate_select_decision.tsv"
    shell:
        "{SG_PYTHON} {SG_SRC}/denovo_candidate_select.py "
        "--contigs {input.contigs} --candidate-qc {input.candidate_qc} "
        "--mode gated_switch --max-span-ratio {SG_MAX_SPAN} "
        "--min-placed-reads {SG_MIN_PLACED} --out-fasta {output.fasta} "
        "--out-report {output.report} --out-decision {output.decision}"


rule junctionQCSelected_standalone:
    input:
        contigs=f"{SG_OUT}/gated_switch.fasta",
        r2=SG_RAW_R2
    output:
        f"{SG_OUT}/junction_qc.tsv"
    shell:
        "{SG_PYTHON} {SG_SRC}/denovo_junction_qc.py "
        "--r2 {input.r2} --contigs {input.contigs} "
        "--max-span-ratio {SG_MAX_SPAN} "
        "--min-local-span-ratio {SG_MIN_LOCAL_SPAN} "
        "--num_processes {SG_NPROC} --out {output}"


rule readbackQC_standalone:
    input:
        contigs=f"{SG_OUT}/gated_switch.fasta",
        r2=SG_FILTERED_R2
    output:
        f"{SG_OUT}/readback_qc.tsv"
    shell:
        "{SG_PYTHON} {SG_SRC}/test/readback_qc_single.py "
        "--r2 {input.r2} --contigs {input.contigs} --n 100000000 "
        "--num_processes {SG_NPROC} --out {output}"


rule uchimeQC_standalone:
    input:
        f"{SG_OUT}/gated_switch.fasta"
    output:
        report=f"{SG_OUT}/uchime_denovo_report.tsv",
        chimeras=f"{SG_OUT}/chimeras.fasta",
        nonchimeras=f"{SG_OUT}/nonchimeras.fasta"
    run:
        if not SG_UCHIME:
            shell("touch {output.report} {output.chimeras} {output.nonchimeras}")
            return
        shell(
            "{SG_VSEARCH} --uchime_denovo {input} --uchimeout {output.report} "
            "--chimeras {output.chimeras} --nonchimeras {output.nonchimeras}"
        )


rule combineQC_standalone:
    input:
        contigs=f"{SG_OUT}/gated_switch.fasta",
        readback=f"{SG_OUT}/readback_qc.tsv",
        uchime=f"{SG_OUT}/uchime_denovo_report.tsv",
        junction=f"{SG_OUT}/junction_qc.tsv"
    output:
        report=f"{SG_OUT}/qc_report.tsv",
        highconf=SG_HIGHCONF,
        flagged=f"{SG_OUT}/gated_switch.flagged.fasta"
    shell:
        "{SG_PYTHON} {SG_SRC}/denovo_qc_combine.py "
        "--contigs {input.contigs} --readback {input.readback} "
        "--uchimeout {input.uchime} --junction {input.junction} "
        "--out-report {output.report} --out-highconf-fasta {output.highconf} "
        "--out-flagged-fasta {output.flagged}"
