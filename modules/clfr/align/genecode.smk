# clfr, se600: library composition against GENCODE (protein_coding / lncRNA / rRNA / globin / ...).
# Selected by params.mrna_mapper == 'minimap_genecode' (workflows/clfr.smk; run_all = composition TSV).
#
# Design: subsample (a ratio needs ~1-5M reads, not the whole library) -> ONE competitive
# minimap2 pass against transcripts + rRNA + genome decoy -> per-read classification.
# The combined reference (transcripts + rRNA + genome decoy) is built once, outside this workflow, by
# build_genecode_ref.sh and shared by all samples (params.genecode_ref, default: next to genecode_transcripts).
# Input is params.genecode_trimmed_fq: an already split + trimmed fastq (a sample of the library is enough); this
# workflow never starts from raw reads and does not run split_reads / trim_reads.
#   select_umis_genecode      -> genecode/sampled_umis.txt     (umi mode only)
#   subsample_reads_genecode  -> genecode/{id}.sub.fq.gz
# Subsample modes (params.genecode_subsample_mode): 'head' (default) = first N reads (fast, file-order
# dependent); 'random' = seqtk read sampling, unbiased for the read-level %; 'umi' = sample whole UMIs (set explicitly), so
# every sampled UMI keeps ALL its reads -- needed only for UMI-level %, since read sampling cuts reads/UMI.
#   map_reads_minimap_genecode-> genecode/{id}.genecode.bam    (read-grouped, unsorted, unmapped kept)
#   classify_reads_genecode   -> genecode/{id}.genecode_composition.tsv (+ _ambiguous.tsv)
#   plot_composition_genecode -> genecode/{id}.genecode_composition.png (barplot, reads % and UMI %)
# Reads and UMIs are both reported: rRNA is often deep per molecule, so the two can differ a lot.
import os

GENECODE_ID = config['samples']['id']
GENECODE_THREADS = config['threads'].get('genecode_map', config['threads'].get('minimap_map', 24))
# The combined reference is sample-independent: build_genecode_ref.sh writes it next to the GENCODE files (same dir
# as params.genecode_transcripts) and every sample reuses it. params.genecode_ref overrides the location.
_gc_tx = config['params'].get('genecode_transcripts', '')
if not _gc_tx and not config['params'].get('genecode_ref'):
    raise ValueError("minimap_genecode needs params.genecode_transcripts (and genecode_gtf) in this run's config.yaml, "
                     "or a prebuilt params.genecode_ref; see config/config.yaml for the genecode_* block")
GENECODE_REF = config['params'].get('genecode_ref') or (
    _gc_tx.replace('.transcripts.fa.gz', '.transcripts_rRNA_decoy.fa.gz') if _gc_tx.endswith('.transcripts.fa.gz')
    else _gc_tx + '.rRNA_decoy.fa.gz')
# The reference is NOT built by this workflow (it is a ~1 GB one-off): run build_genecode_ref.sh once first.
if not os.path.exists(GENECODE_REF):
    raise ValueError("GENCODE combined reference not found: {}\nBuild it once with "
                     "modules/clfr/align/build_genecode_ref.sh (transcripts + rRNA + genome decoy), or point "
                     "params.genecode_ref at an existing file".format(GENECODE_REF))


GENECODE_FQ = config['params'].get('genecode_trimmed_fq', '')
if not GENECODE_FQ or not os.path.exists(GENECODE_FQ):
    raise ValueError("params.genecode_trimmed_fq must point to an existing split + trimmed fastq(.gz): "
                     "'{}'".format(GENECODE_FQ))

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
        min_reads = config['params'].get('genecode_umi_min_reads', 25),
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
        fq = "genecode/{}.sub.fq.gz".format(GENECODE_ID)
    params:
        n = config['params'].get('genecode_subsample_reads', 2000000),
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


rule map_reads_minimap_genecode:
    input:
        ref = GENECODE_REF,
        fq = "genecode/{}.sub.fq.gz".format(GENECODE_ID)
    output:
        bam = "genecode/{}.genecode.bam".format(GENECODE_ID)
    threads:
        GENECODE_THREADS
    params:
        MINIMAP = config['params']['minimap2'],
        max_hits = config['params'].get('genecode_max_hits', 50)
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
        amb = "genecode/{}.genecode_ambiguous.tsv".format(GENECODE_ID)
    params:
        script = config['params']['src_dir'] + "modules/clfr/align/genecode_classify.py",
        gtf = config['params'].get('genecode_gtf', ''),
        min_cov = config['params'].get('genecode_min_cov', 0.8)
    shell:
        """
        set -euo pipefail
        gtf_arg=""
        if [ -n "{params.gtf}" ]; then gtf_arg="--gtf {params.gtf}"; fi
        python3 {params.script} --bam {input.bam} --out {output.tsv} --out_ambiguous {output.amb} \
            --min_cov {params.min_cov} $gtf_arg
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
