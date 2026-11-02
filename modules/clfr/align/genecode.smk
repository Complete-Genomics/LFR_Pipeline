# clfr, se600: library composition against GENCODE (protein_coding / lncRNA / rRNA / globin / ...).
# Independent module: enabled by modules.genecode (workflows/clfr.smk adds the composition TSV/PNG to run_all).
#
# Design: subsample (a ratio needs ~1-5M reads, not the whole library) -> ONE competitive
# minimap2 pass against transcripts + rRNA + genome decoy -> per-read classification.
# The combined reference (transcripts + rRNA + genome decoy) is built once, outside this workflow, by
# build_genecode_ref.sh and shared by all samples (params.genecode_ref, default: next to genecode_transcripts).
# Input, picked when the workflow is parsed (first match wins):
#   1. params.genecode_trimmed_fq (explicit path)  \
#   2. data/split_read_2_trimmed.fastq.gz exists   /  -> already trimmed: sample it as is, no trim
#   3. otherwise data/split_read.2.fq.gz (SE600)       -> sample first, then trim only the sample (trim_sub_genecode, same
#                                                         params.adapter_ref*/trim_* settings as trim_reads), so the whole
#                                                         library is never trimmed
#   select_umis_genecode      -> genecode/sampled_umis.txt     (umi mode only)
#   subsample_reads_genecode  -> genecode/{id}.sub.fq.gz          (the sampled reads; already trimmed in case 1/2)
#   trim_sub_genecode         -> genecode/{id}.sub.trimmed.fq.gz  (case 3 only; reads dropped by minlen/maxns/maq make it a
#                                                                  bit smaller than genecode_subsample_reads)
#   map_reads_minimap_genecode reads the trimmed sample: sub.trimmed.fq.gz (case 3) or sub.fq.gz (case 1/2)
# Subsample modes (params.genecode_subsample_mode): 'head' (default) = first N reads (fast, file-order
# dependent); 'random' = seqtk read sampling, unbiased for the read-level %; 'umi' = sample whole UMIs (set explicitly), so
# every sampled UMI keeps ALL its reads -- needed only for UMI-level %, since read sampling cuts reads/UMI.
#   map_reads_minimap_genecode-> genecode/{id}.genecode.bam    (read-grouped, unsorted, unmapped kept)
#   classify_reads_genecode   -> genecode/{id}.genecode_composition.tsv (+ _ambiguous.tsv, _clip.tsv)
#   plot_composition_genecode -> genecode/{id}.genecode_composition.png (barplot, reads % and UMI %)
# Reads and UMIs are both reported: rRNA is often deep per molecule, so the two can differ a lot.
import os
import sys

GENECODE_ID = config['samples']['id']
GENECODE_THREADS = config['threads'].get('genecode_map', config['threads'].get('minimap_map', 24))
# The combined reference is sample-independent: build_genecode_ref.sh writes it next to the GENCODE files (same dir
# as params.genecode_transcripts) and every sample reuses it. params.genecode_ref overrides the location.
_gc_tx = config['params'].get('genecode_transcripts', '')
if not _gc_tx and not config['params'].get('genecode_ref'):
    raise ValueError("modules.genecode needs params.genecode_transcripts (and genecode_gtf) in this run's config.yaml, "
                     "or a prebuilt params.genecode_ref; see config/config.yaml for the genecode_* block")
GENECODE_REF = config['params'].get('genecode_ref') or (
    _gc_tx.replace('.transcripts.fa.gz', '.transcripts_rRNA_decoy.fa.gz') if _gc_tx.endswith('.transcripts.fa.gz')
    else _gc_tx + '.rRNA_decoy.fa.gz')
# The reference is NOT built by this workflow (it is a ~1 GB one-off): run build_genecode_ref.sh once first.
if not os.path.exists(GENECODE_REF):
    raise ValueError("GENCODE combined reference not found: {}\nBuild it once with "
                     "modules/clfr/align/build_genecode_ref.sh (transcripts + rRNA + genome decoy), or point "
                     "params.genecode_ref at an existing file".format(GENECODE_REF))


# See the header for the three cases. A given genecode_trimmed_fq must exist.
GENECODE_TRIMMED_GIVEN = bool(config['params'].get('genecode_trimmed_fq'))
if GENECODE_TRIMMED_GIVEN and not os.path.exists(config['params']['genecode_trimmed_fq']):
    raise ValueError("params.genecode_trimmed_fq must point to an existing split + trimmed fastq(.gz): "
                     "'{}'".format(config['params']['genecode_trimmed_fq']))
GENECODE_PIPELINE_TRIMMED = "data/split_read_2_trimmed.fastq.gz"
GENECODE_USE_TRIMMED = GENECODE_TRIMMED_GIVEN or os.path.exists(GENECODE_PIPELINE_TRIMMED)
if GENECODE_USE_TRIMMED:
    GENECODE_FQ = config['params'].get('genecode_trimmed_fq') or GENECODE_PIPELINE_TRIMMED
else:
    if config['params']['sequence_type'].lower() != 'se':
        raise ValueError("modules.genecode trims the sample itself for SE600 only (params.sequence_type: se)")
    GENECODE_FQ = "data/split_read.2.fq.gz"
print("genecode input: {} ({})".format(GENECODE_FQ, "already trimmed, no trim" if GENECODE_USE_TRIMMED
                                      else "untrimmed: sample, then trim_sub_genecode"), file=sys.stderr)
GENECODE_SUB = "genecode/{}.sub.fq.gz".format(GENECODE_ID)                  # sampled reads
GENECODE_SUB_TRIMMED = GENECODE_SUB if GENECODE_USE_TRIMMED else "genecode/{}.sub.trimmed.fq.gz".format(GENECODE_ID)  # what is mapped

GENECODE_SUB_MODE = config['params'].get('genecode_subsample_mode', 'head')
if GENECODE_SUB_MODE not in ('umi', 'head', 'random'):
    raise ValueError("params.genecode_subsample_mode must be umi, head or random")


rule select_umis_genecode:
    # Same input as denovo_preprocess.select_denovo_barcodes (split_stat_read1.log: idx, reads, UMI),
    # then a seeded reservoir sample of K UMIs -> reproducible UMI-level subsample.
    input:
        "split_stat_read1.log"
    output:
        "genecode/sampled_umis.txt"
    params:
        min_reads = config['params']['reads_per_BC'],
        n_umis = config['params'].get('genecode_umi_n', 50000),
        seed = config['params'].get('genecode_subsample_seed', 100)
    shell:
        """
        mkdir -p genecode
        awk -F '\\t' -v cutoff={params.min_reads} \
            'NR > 4 && NF >= 3 && $2 + 0 >= cutoff && $3 ~ /^[ACGT]+$/ {{print $3}}' {input} \
        | awk -v k={params.n_umis} -v s={params.seed} \
            'BEGIN {{srand(s)}} {{n++; if (n <= k) r[n] = $0; else {{j = int(rand() * n) + 1; if (j <= k) r[j] = $0}}}} \
             END {{m = (n < k ? n : k); for (i = 1; i <= m; i++) print "BX:Z:" r[i]; \
                  print "eligible UMIs: " n ", sampled: " m > "/dev/stderr"}}' \
        > {output}
        [ -s {output} ]
        """


rule subsample_reads_genecode:
    input:
        fq = GENECODE_FQ,
        umis = "genecode/sampled_umis.txt" if GENECODE_SUB_MODE == 'umi' else []
    output:
        fq = GENECODE_SUB
    params:
        n = config['params'].get('genecode_subsample_reads', 400000),
        seed = config['params'].get('genecode_subsample_seed', 100),
        mode = GENECODE_SUB_MODE,
        seqtk = config['params'].get('seqtk', 'seqtk'),
        bgzip = config['params'].get('bgzip', 'gzip')
    shell:
        """
        set -euo pipefail
        N={params.n}
        if [ "{params.mode}" = "umi" ]; then
            # One full pass: reads of a UMI are scattered through the file, so no early exit is possible.
            # Header is "<id>#<bc>#/2<TAB>BX:Z:<bc>" (split_barcode_LFR.py); same filter as denovo filter_reads2.
            gzip -dc {input.fq} \
            | awk -F '\\t' 'NR == FNR {{keep[$1]; next}} FNR % 4 == 1 {{ok = ($2 in keep)}} ok' {input.umis} - \
            | {params.bgzip} -c > {output.fq}
        elif [ "{params.mode}" = "random" ]; then
            {params.seqtk} sample -s {params.seed} {input.fq} $N | gzip > {output.fq}
        else
            # head closes the pipe early, so gzip -dc gets SIGPIPE: relax pipefail for this line only
            ( set +o pipefail; gzip -dc {input.fq} | head -n $((4*N)) | gzip > {output.fq} )
        fi
        n_out=$(( $(gzip -dc {output.fq} | wc -l) / 4 ))
        echo "subsampled reads: $n_out (mode {params.mode})" >&2
        [ "$n_out" -gt 0 ]
        """


if not GENECODE_USE_TRIMMED:
    # Same bbduk recipe and the same config switches as splitreads.smk trim_reads (SE only), applied to the sample:
    # up to four bbduk stages chained through pipes (5' contaminant -> 3' linker -> 3' adapter + qtrim -> polyA); only the
    # last writes to disk. Stage order matters: polyA is only a tail once the 3' adapter behind it is removed.
    rule trim_sub_genecode:
        input:
            GENECODE_SUB
        output:
            GENECODE_SUB_TRIMMED
        params:
            bbduk = config['params']['bbduk'],
            adapter_ref = config['params'].get('adapter_ref', ''),
            adapter_ref_k = int(config['params'].get('adapter_ref_k', 23)),
            trim_5p_contam = str(config['params'].get('adapter_ref_5p', False)).strip().lower() in ("true", "1", "yes"),
            adapter_ref_5p_file = config['params'].get('adapter_ref_5p_file', ''),
            trim_3p_linker = str(config['params'].get('adapter_ref_3p', False)).strip().lower() in ("true", "1", "yes"),
            adapter_ref_3p_file = config['params'].get('adapter_ref_3p_file', ''),
            trim_polya = int(config['params'].get('trim_polya', 0)),
            trim_q = int(config['params'].get('trim_q', 6)),
            trim_maq = int(config['params'].get('trim_maq', 0)),
            trim_mgi = str(config['params'].get('trim_mgi_adapter', True)).strip().lower() in ("true", "1", "yes"),
            trim_qtrim = str(config['params'].get('trim_qtrim', True)).strip().lower() in ("true", "1", "yes"),
            trim_5p_restrict_left = int(config['params'].get('trim_5p_restrict_left', 30)),
            trim_minlen = int(config['params'].get('trim_minlen', 10)),
            trim_maxns = int(config['params'].get('trim_maxns', 0))
        benchmark:
            "Benchmarks/main.trim_sub_genecode.txt"
        shell:
            """
            ADAPTER_ARGS=""
            if [[ -n "{params.adapter_ref}" && "{params.trim_mgi}" == "True" ]]; then
                ADAPTER_ARGS="ktrim=r ref={params.adapter_ref} k={params.adapter_ref_k} hdist=1"
            fi
            ORIG_IN="in={input}"; PIPE_IN="in=stdin.fq int=f"; MID_OUT="out=stdout.fq"; FINAL_OUT="out={output}"
            MAIN_ARGS=""
            if [[ "{params.trim_qtrim}" == "True" ]]; then MAIN_ARGS="qtrim=rl trimq={params.trim_q}"; fi
            if [[ {params.trim_maxns} -ge 0 ]]; then MAIN_ARGS="$MAIN_ARGS maxns={params.trim_maxns}"; fi
            MAIN_ARGS="$MAIN_ARGS minlen={params.trim_minlen}"
            if [[ {params.trim_maq} -gt 0 ]]; then MAIN_ARGS="$MAIN_ARGS maq={params.trim_maq}"; fi

            CMD=""
            MAIN_IN="$ORIG_IN"
            if [[ "{params.trim_5p_contam}" == "True" ]]; then
                CMD="{params.bbduk} $ORIG_IN $MID_OUT ktrim=l k=17 mink=11 hdist=1 restrictleft={params.trim_5p_restrict_left} ref={params.adapter_ref_5p_file} | "
                MAIN_IN="$PIPE_IN"
            fi
            if [[ "{params.trim_3p_linker}" == "True" ]]; then
                CMD="${{CMD}}{params.bbduk} $MAIN_IN $MID_OUT ktrim=r k=12 hdist=0 ref={params.adapter_ref_3p_file} | "
                MAIN_IN="$PIPE_IN"
            fi
            if [[ {params.trim_polya} -gt 0 ]]; then
                CMD="${{CMD}}{params.bbduk} $MAIN_IN $MID_OUT $MAIN_ARGS $ADAPTER_ARGS | {params.bbduk} $PIPE_IN $FINAL_OUT trimpolya={params.trim_polya}"
            else
                CMD="${{CMD}}{params.bbduk} $MAIN_IN $FINAL_OUT $MAIN_ARGS $ADAPTER_ARGS"
            fi
            eval "$CMD"
            """


rule map_reads_minimap_genecode:
    input:
        ref = GENECODE_REF,
        fq = GENECODE_SUB_TRIMMED
    output:
        bam = "genecode/{}.genecode.bam".format(GENECODE_ID)
    threads:
        GENECODE_THREADS
    params:
        MINIMAP = config['params']['minimap2'],
        max_hits = config['params'].get('genecode_max_hits', 20)
    log:
        "genecode/minimap_genecode.log"
    benchmark:
        "Benchmarks/main.map_reads_genecode.txt"
    shell:
        """
        set -euo pipefail
        # No sort: minimap2 emits each read's records together, which the classifier needs.
        # No --sam-hit-only: unmapped reads stay in the BAM so they count in the denominator.
        # -N/-p keep near-best secondary hits so reads hitting several classes are flagged ambiguous.
        {params.MINIMAP} -ax map-ont -k 15 -w 10 \
            -N {params.max_hits} -p 0.9 --secondary=yes \
            -t {threads} \
            {input.ref} {input.fq} 2>> {log} \
        | samtools view -b -@ 4 -o {output.bam} - 2>> {log}
        """


rule classify_reads_genecode:
    input:
        bam = "genecode/{}.genecode.bam".format(GENECODE_ID)
    output:
        tsv = "genecode/{}.genecode_composition.tsv".format(GENECODE_ID),
        amb = "genecode/{}.genecode_ambiguous.tsv".format(GENECODE_ID),
        clip = "genecode/{}.genecode_clip.tsv".format(GENECODE_ID)
    params:
        script = config['params']['src_dir'] + "modules/clfr/align/genecode_classify.py",
        gtf = config['params'].get('genecode_gtf', ''),
        min_cov = config['params'].get('genecode_min_cov', 0.5)
    shell:
        """
        set -euo pipefail
        gtf_arg=""
        if [ -n "{params.gtf}" ]; then gtf_arg="--gtf {params.gtf}"; fi
        python3 {params.script} --bam {input.bam} --out {output.tsv} --out_ambiguous {output.amb} \
            --out_clip {output.clip} --min_cov {params.min_cov} $gtf_arg
        cat {output.tsv} >&2
        """


rule plot_composition_genecode:
    input:
        tsv = "genecode/{}.genecode_composition.tsv".format(GENECODE_ID)
    output:
        png = "genecode/{}.genecode_composition.png".format(GENECODE_ID)
    params:
        script = config['params']['src_dir'] + "modules/clfr/align/genecode_plot.py",
        title = "GENCODE composition: {}".format(config['samples'].get('name', GENECODE_ID))
    shell:
        "python3 {params.script} --tsv {input.tsv} --out {output.png} --title '{params.title}'"
