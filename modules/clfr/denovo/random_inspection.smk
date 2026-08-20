# Standalone offline QA workflow: does NOT run as part of a production sample
# (nothing in denovo_olc.smk/denovo_preprocess.smk depends on this file, and
# it does not touch qc_setting/resolve_qc -- those are live per-sample
# runtime decisions, this is a benchmark). Purpose: given a raw barcode-
# tagged FASTQ, assemble it three ways -- no read filtering, the current
# default filter, and the ML-scored tie-break -- and compare the default and
# ML arms against a "local ground truth" reference assigned from the
# unfiltered arm via vsearch against a broad 16S database. Same method used
# for every denovo.md sec 51/59/62/67 result; this is that method turned into
# a rerunnable rule set instead of one-off scripts, so it can be pointed at
# a new sample without re-deriving it from scratch (denovo.md sec 76).
#
# Usage: snakemake -s random_inspection.smk --configfile <config.yaml> -j <N> \
#            random_inspection/report.md
#
# Config (under `random_inspection:`):
#   input_fastq      barcode-tagged FASTQ(.gz), e.g. data/split_read_2_trimmed.fastq.gz
#   outdir           default "random_inspection"
#   max_umis         default 3000; barcode-sorted input cap shared by all arms.
#                    Set "all" only for an intentional full-run benchmark.
#   ml_model         default modules/clfr/denovo/models/model_identity.lgb
#   greengenes_db    default modules/clfr/denovo/db/gg.fna (NOT shipped -- see db/README.md)
#   vsearch, samtools  binary paths, default to PATH lookup
#   num_processes    default 4
#   vsearch_id       reference-assignment identity floor, default 0.75 (denovo.md sec 59)
#   chimera_reference optional small, known-composition control FASTA (e.g.
#                     Zymo references). When supplied, report strict
#                     quarter-split chimera rates; it must not be Greengenes.
#   chimera_label_separator source-label delimiter in that FASTA's IDs,
#                     default "_16S" (Zymo: Bacillus_subtilis_16S_1)
# Also reuses config['params']['general_python'] and config['params']['src_dir'].

import os

RI = config.get('random_inspection', {})
RI_OUT = RI.get('outdir', 'random_inspection')
RI_SRC = config['params']['src_dir'] + "/modules/clfr/denovo"
RI_PYTHON = config['params']['general_python']
RI_ML_MODEL = RI.get('ml_model', config['params']['src_dir'] +
                      "/modules/clfr/denovo/models/model_identity.lgb")
RI_GG_DB = RI.get('greengenes_db', config['params']['src_dir'] +
                   "/modules/clfr/denovo/db/gg.fna")
RI_VSEARCH = RI.get('vsearch', 'vsearch')
RI_SAMTOOLS = RI.get('samtools', 'samtools')
RI_NPROC = RI.get('num_processes', 4)
RI_VSEARCH_ID = RI.get('vsearch_id', 0.75)
RI_CHIMERA_REFERENCE = RI.get('chimera_reference', '')
RI_CHIMERA_LABEL_SEPARATOR = RI.get('chimera_label_separator', '_16S')
RI_MAX_UMIS = RI.get('max_umis', 3000)
if RI_MAX_UMIS not in (None, "", "all"):
    RI_MAX_UMIS = int(RI_MAX_UMIS)
    if RI_MAX_UMIS < 1:
        raise ValueError("random_inspection.max_umis must be a positive integer or 'all'")

wildcard_constraints:
    arm = "nofilter|plain|ml"

# .strip() matters: this string is concatenated with " {input} > {output}",
# and a trailing newline would end the awk command before its file argument,
# leaving awk reading stdin and the filename on its own line as a command to
# execute. Braces are doubled because Snakemake format()s shell strings.
_EXTRACT_PRIMARY_AWK = r"""
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


rule random_inspection_all:
    input:
        f"{RI_OUT}/report.md"


rule sortFastq_inspect:
    input:
        RI.get('input_fastq', 'data/split_read_2_trimmed.fastq.gz')
    output:
        f"{RI_OUT}/sorted.fastq.gz"
    benchmark:
        "Benchmarks/random_inspection.sortFastq.txt"
    shell:
        "bash {RI_SRC}/denovo_sort_fastq.sh {input} {output} {RI_NPROC} {RI_SAMTOOLS}"


## Cap after barcode sort, before feature extraction and all three arms.  This
## keeps the denominator and the exact UMI cohort identical across nofilter,
## plain, and ML, while avoiding a full-data assembly/reference comparison.
rule capUmis_inspect:
    input:
        f"{RI_OUT}/sorted.fastq.gz"
    output:
        f"{RI_OUT}/subset.fastq.gz"
    params:
        max_umis=RI_MAX_UMIS
    shell:
        """
        if [ "{params.max_umis}" = "all" ] || [ -z "{params.max_umis}" ]; then
            cp {input} {output}
        else
            gzip -dc {input} | awk -v max_umis={params.max_umis} '
                NR % 4 == 1 {{
                    barcode = $2
                    if (barcode != previous) {{
                        previous = barcode
                        seen++
                    }}
                    if (seen > max_umis) exit
                }}
                {{ print }}
            ' | gzip -c > {output}
        fi
        """


rule extractFeatures_inspect:
    input:
        f"{RI_OUT}/subset.fastq.gz"
    output:
        f"{RI_OUT}/features.tsv"
    benchmark:
        "Benchmarks/random_inspection.extractFeatures.txt"
    shell:
        "{RI_PYTHON} {RI_SRC}/denovo_read_features.py "
        "--fastq {input} --out {output} --num_processes {RI_NPROC}"


## Same transform reformat_fasta2 (denovo_preprocess.smk) applies -- this arm
## IS the unfiltered read set, used only to assign each barcode's reference,
## never compared directly (denovo.md sec 59's "neutral w.r.t. the arms
## being compared" requirement).
rule buildSgrep_inspect:
    input:
        f"{RI_OUT}/subset.fastq.gz"
    output:
        f"{RI_OUT}/nofilter.tsv"
    shell:
        """
        gzip -dc {input} | \
        awk '{{if (NR%4==1) {{temp=$1; $1=$2; $2=temp}} }}1' | \
        awk '{{if (NR%4==1) line=line$0"\\t"; if (NR%4==2) {{print line$0; line=""}}}}' \
        > {output}
        """


rule readFilterPlain_inspect:
    input:
        f"{RI_OUT}/nofilter.tsv"
    output:
        f"{RI_OUT}/plain.tsv"
    benchmark:
        "Benchmarks/random_inspection.readFilterPlain.txt"
    shell:
        "{RI_PYTHON} {RI_SRC}/denovo_read_filter.py --r2 {input} --out {output} "
        "--num_processes {RI_NPROC}"


rule readFilterML_inspect:
    input:
        reads = f"{RI_OUT}/nofilter.tsv",
        features = f"{RI_OUT}/features.tsv"
    output:
        f"{RI_OUT}/ml.tsv"
    benchmark:
        "Benchmarks/random_inspection.readFilterML.txt"
    shell:
        "{RI_PYTHON} {RI_SRC}/denovo_read_filter.py --r2 {input.reads} --out {output} "
        "--num_processes {RI_NPROC} --ml-model {RI_ML_MODEL} --ml-features {input.features}"


## denovo_seed_olc.py writes to a CWD-relative denovo/final_contigs_0.fa (no
## --out flag), so each arm gets its own subdirectory to run inside of.
rule assembleArm_inspect:
    input:
        f"{RI_OUT}/{{arm}}.tsv"
    output:
        f"{RI_OUT}/{{arm}}/denovo/final_contigs_0.fa"
    benchmark:
        "Benchmarks/random_inspection.assembleArm.{arm}.txt"
    params:
        armdir = lambda wc: f"{RI_OUT}/{wc.arm}"
    shell:
        "mkdir -p {params.armdir} && "
        "cd {params.armdir} && "
        "{RI_PYTHON} {RI_SRC}/denovo_seed_olc.py --sequence_type se "
        "--num_processes {RI_NPROC} --r2 $OLDPWD/{input}"


## One representative contig per barcode (the k41_0 / "primary" convention
## denovo_olc.smk's filterOLC_longest already uses), renamed to the more
## legible denovo.fasta alongside it.
rule extractPrimary_inspect:
    input:
        f"{RI_OUT}/{{arm}}/denovo/final_contigs_0.fa"
    output:
        primary = f"{RI_OUT}/{{arm}}/primary.fa",
        aliased = f"{RI_OUT}/{{arm}}/denovo.fasta"
    shell:
        _EXTRACT_PRIMARY_AWK + " {input} > {output.primary} && cp {output.primary} {output.aliased}"


## The one step that touches the full reference database (~1.78GB; a prior
## 200k-query run against it peaked at 3.7GB resident -- denovo.md sec 76).
## Only the nofilter arm goes through this; plain/ml are scored against the
## small subset it selects (scoreIdentity_inspect below), never against the
## full db directly.
rule assignReference_inspect:
    input:
        query = f"{RI_OUT}/nofilter/primary.fa",
        db = RI_GG_DB
    output:
        f"{RI_OUT}/nofilter_hits.tsv"
    benchmark:
        "Benchmarks/random_inspection.assignReference.txt"
    shell:
        "{RI_VSEARCH} --usearch_global {input.query} --db {input.db} "
        "--id {RI_VSEARCH_ID} --strand both --maxaccepts 1 --maxrejects 32 --top_hits_only "
        "--userout {output} --userfields query+target+id "
        "--threads {RI_NPROC} --log {output}.log"


rule extractAssignedRefs_inspect:
    input:
        hits = f"{RI_OUT}/nofilter_hits.tsv",
        db = RI_GG_DB
    output:
        f"{RI_OUT}/assigned_refs.fa"
    shell:
        """
        {RI_SAMTOOLS} faidx {input.db}
        cut -f2 {input.hits} | sort -u > {output}.ids
        {RI_SAMTOOLS} faidx {input.db} -r {output}.ids > {output}
        """


## Small db (one sequence per distinct assigned reference) -- safe, fast,
## no relation to the big-db memory risk above. --maxaccepts/--maxrejects 0
## (unlimited) because the row this needs is the one matching each barcode's
## PRE-assigned target, which need not be this arm's own top hit.
rule scoreIdentity_inspect:
    input:
        query = f"{RI_OUT}/{{arm}}/primary.fa",
        db = f"{RI_OUT}/assigned_refs.fa"
    output:
        f"{RI_OUT}/{{arm}}_hits.tsv"
    wildcard_constraints:
        arm = "plain|ml"
    shell:
        "{RI_VSEARCH} --usearch_global {input.query} --db {input.db} "
        "--id 0.5 --strand both --maxaccepts 0 --maxrejects 0 "
        "--userout {output} --userfields query+target+id "
        "--threads {RI_NPROC} --log {output}.log"


rule pairedCompare_inspect:
    input:
        ref = f"{RI_OUT}/nofilter_hits.tsv",
        plain = f"{RI_OUT}/plain_hits.tsv",
        ml = f"{RI_OUT}/ml_hits.tsv"
    output:
        f"{RI_OUT}/report.tsv"
    shell:
        "{RI_PYTHON} {RI_SRC}/random_inspection_compare.py "
        "--reference-hits {input.ref} --arm-a-hits {input.plain} --arm-b-hits {input.ml} "
        "--arm-a-name plain --arm-b-name ml --out {output}"


## A strict chimera rate is meaningful only for a known-composition control
## with a small, explicitly labelled reference. It is deliberately separate
## from the broad Greengenes database used for local-reference identity.
rule splitPrimaryQuarters_inspect:
    input:
        f"{RI_OUT}/{{arm}}/primary.fa"
    output:
        f"{RI_OUT}/{{arm}}/primary_quarters.fa"
    shell:
        "{RI_PYTHON} {RI_SRC}/random_inspection_quarters.py --input {input} --out {output}"


rule scoreQuarterChimera_inspect:
    input:
        query=f"{RI_OUT}/{{arm}}/primary_quarters.fa",
        db=RI_CHIMERA_REFERENCE
    output:
        f"{RI_OUT}/{{arm}}_quarter_hits.tsv"
    shell:
        "{RI_VSEARCH} --usearch_global {input.query} --db {input.db} --id 0.75 "
        "--strand both --maxaccepts 0 --maxrejects 0 "
        "--userout {output} --userfields query+target+id --threads {RI_NPROC} "
        "--log {output}.log"


def ri_optional_chimera_inputs(wildcards):
    if not RI_CHIMERA_REFERENCE:
        return []
    return [f"{RI_OUT}/{arm}_quarter_hits.tsv" for arm in ("nofilter", "plain", "ml")]


def ri_optional_chimera_args():
    if not RI_CHIMERA_REFERENCE:
        return ""
    return (
        f" --quarter-hits-nofilter {RI_OUT}/nofilter_quarter_hits.tsv"
        f" --quarter-hits-plain {RI_OUT}/plain_quarter_hits.tsv"
        f" --quarter-hits-ml {RI_OUT}/ml_quarter_hits.tsv"
        f" --chimera-label-separator {RI_CHIMERA_LABEL_SEPARATOR}"
    )


rule randomInspectionReport_inspect:
    input:
        paired=f"{RI_OUT}/report.tsv",
        nofilter_reads=f"{RI_OUT}/nofilter.tsv",
        plain_reads=f"{RI_OUT}/plain.tsv",
        ml_reads=f"{RI_OUT}/ml.tsv",
        nofilter_contigs=f"{RI_OUT}/nofilter/primary.fa",
        plain_contigs=f"{RI_OUT}/plain/primary.fa",
        ml_contigs=f"{RI_OUT}/ml/primary.fa",
        chimera=ri_optional_chimera_inputs
    output:
        report=f"{RI_OUT}/report.md",
        metrics=f"{RI_OUT}/summary_metrics.tsv",
        identity_hist=f"{RI_OUT}/identity_histogram.tsv",
        contig_len_hist=f"{RI_OUT}/contig_length_histogram.tsv",
        reads_hist=f"{RI_OUT}/reads_per_umi_histogram.tsv"
    params:
        chimera_args=ri_optional_chimera_args()
    shell:
        "{RI_PYTHON} {RI_SRC}/random_inspection_report.py "
        "--nofilter-reads {input.nofilter_reads} --plain-reads {input.plain_reads} --ml-reads {input.ml_reads} "
        "--nofilter-contigs {input.nofilter_contigs} --plain-contigs {input.plain_contigs} --ml-contigs {input.ml_contigs} "
        "--paired-identity {input.paired} --out-md {output.report} --out-metrics {output.metrics} "
        "--out-identity-hist {output.identity_hist} --out-contig-len-hist {output.contig_len_hist} "
        "--out-reads-hist {output.reads_hist}{params.chimera_args}"
