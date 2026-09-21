# clfr, t1 se600, align with minimap
import os

MINIMAP_TOTAL_THREADS = config['threads'].get('minimap_total', 32)
MINIMAP_MAP_THREADS = config['threads'].get('minimap_map', 24)
MINIMAP_SORT_THREADS = config['threads'].get('minimap_sort', 8)
if MINIMAP_MAP_THREADS + MINIMAP_SORT_THREADS > MINIMAP_TOTAL_THREADS:
    raise ValueError(
        "threads.minimap_map + threads.minimap_sort must not exceed "
        "threads.minimap_total"
    )

rule map_reads_minimap:
    input:
        ref = REF,
        # fqs = expand("{sample}", sample=config['samples']['fastq'])
        # fq1 = 'data/split_read.1.fq.gz',
        fq2 = 'data/split_read_2_trimmed.fastq.gz'
    output:
        bam = "keep/Align/{}.sort.bam".format(config['samples']['id']), # may want this to be a temp file
        bai = "keep/Align/{}.sort.bam.bai".format(config['samples']['id']),
    threads:
        MINIMAP_TOTAL_THREADS
    params:
        MINIMAP = config['params']['minimap2'],
        anno_bed = config['params']['minimap2_anno_bed'],
        STAR = config['params']['star'],
        star_index = config['params']['star_index'],
        hisat2 = config['params']['hisat2'],
        hisat2_index = config['params']['hisat2_index'],
        hisat2_splicesites = config['params']['hisat2_splicesites'],
        bwa_mem = config['params']['bwa_mem'],
        gatk_install = config['params']['gatk_install'],
        reference = config['params'].get('minimap2_index') or REF,
        map_threads = MINIMAP_MAP_THREADS,
        sort_threads = MINIMAP_SORT_THREADS,
        sort_mem = config['params'].get('minimap_sort_mem', '2G'),
        readgroup = r'@RG\tID:{0}\tSM:{0}\tPL:{1}'.format(config['samples']['id'],
                                                          config['params']['platform'])
    benchmark:
        "Benchmarks/main.map_reads.txt"
    shell:
        """
        set -euo pipefail

        # ============ 1. set env ============
        #log=Align/minimap2.log
        SORT_TMP=Align/tmp/minimap_tmp_$$
        echo "minimap2 samtools sort tmp: $SORT_TMP" >&2
        mkdir -p $SORT_TMP
        trap 'rm -rf "$SORT_TMP"' EXIT


        # ============ 2. minimap2 mapping ============
        {params.MINIMAP} -ax splice:sr \
            -t {params.map_threads} \
            --secondary=no \
            --sam-hit-only \
            --junc-bed {params.anno_bed} \
            {params.reference} {input.fq2} \
            2>> minimap.log \
        | samtools sort -@ {params.sort_threads} -m {params.sort_mem} \
            -T $SORT_TMP/sort \
            -o {output.bam} - \
            2>> minimap.log

        # ============ 3. index ============
        samtools index -@ {params.sort_threads} {output.bam} {output.bai} 2>> minimap.log

        """
