# clfr, t1 se600, align with minimap -- chunked variant of align.minimap.smk.
# Benchmarks on main.map_reads.txt showed near-linear scaling with thread
# count (43.8h @ t=15 vs 27.7h @ t=24, predicted 44.3h), meaning the step is
# genuinely compute-bound rather than I/O-bound. Splitting the input fastq
# into independent chunks lets those threads be spread across multiple
# Slurm jobs/nodes instead of being capped by one node's core count, while
# producing the same keep/Align/{id}.sort.bam(.bai) outputs as the
# unsplit rule so downstream rules are unaffected.
import os

MINIMAP_N_CHUNKS = config['threads'].get('minimap_chunks', 4)
MINIMAP_CHUNK_TOTAL_THREADS = config['threads'].get('minimap_chunk_total', 8)
MINIMAP_CHUNK_MAP_THREADS = config['threads'].get('minimap_chunk_map', 6)
MINIMAP_CHUNK_SORT_THREADS = config['threads'].get('minimap_chunk_sort', 2)
if MINIMAP_CHUNK_MAP_THREADS + MINIMAP_CHUNK_SORT_THREADS > MINIMAP_CHUNK_TOTAL_THREADS:
    raise ValueError(
        "threads.minimap_chunk_map + threads.minimap_chunk_sort must not exceed "
        "threads.minimap_chunk_total"
    )

wildcard_constraints:
    chunk = r"\d+"

rule split_fastq_for_mapping:
    input:
        fq2 = 'data/split_read_2_trimmed.fastq.gz'
    output:
        expand("Align/tmp/minimap_chunks/chunk_{chunk}.fastq.gz", chunk=range(MINIMAP_N_CHUNKS))
    params:
        n = MINIMAP_N_CHUNKS,
        prefix = "Align/tmp/minimap_chunks/chunk_"
    threads:
        config['threads'].get('minimap_split', 4)
    benchmark:
        "Benchmarks/main.split_fastq_for_mapping.txt"
    shell:
        """
        set -euo pipefail
        mkdir -p Align/tmp/minimap_chunks
        gzip -dc {input.fq2} | awk -v n={params.n} -v prefix={params.prefix} '
        {{
            rec = rec $0 "\\n"
            if (NR % 4 == 0) {{
                idx = int((NR - 1) / 4) % n
                printf "%s", rec | ("gzip -c > " prefix idx ".fastq.gz")
                rec = ""
            }}
        }}
        END {{
            for (i = 0; i < n; i++) close("gzip -c > " prefix i ".fastq.gz")
        }}'
        """

rule map_minimap_chunk:
    input:
        ref = REF,
        fq = "Align/tmp/minimap_chunks/chunk_{chunk}.fastq.gz"
    output:
        bam = "Align/tmp/minimap_chunks/chunk_{chunk}.sort.bam"
    threads:
        MINIMAP_CHUNK_TOTAL_THREADS
    params:
        MINIMAP = config['params']['minimap2'],
        anno_bed = config['params']['minimap2_anno_bed'],
        reference = config['params'].get('minimap2_index') or REF,
        map_threads = MINIMAP_CHUNK_MAP_THREADS,
        sort_threads = MINIMAP_CHUNK_SORT_THREADS,
        sort_mem = config['params'].get('minimap_sort_mem', '2G'),
    benchmark:
        "Benchmarks/main.map_reads_chunk.{chunk}.txt"
    shell:
        """
        set -euo pipefail

        SORT_TMP=Align/tmp/minimap_tmp_{wildcards.chunk}_$$
        mkdir -p $SORT_TMP
        trap 'rm -rf "$SORT_TMP"' EXIT

        {params.MINIMAP} -ax splice:sr \
            -t {params.map_threads} \
            --secondary=no \
            --sam-hit-only \
            --junc-bed {params.anno_bed} \
            {params.reference} {input.fq} \
            2>> minimap.chunk{wildcards.chunk}.log \
        | samtools sort -@ {params.sort_threads} -m {params.sort_mem} \
            -T $SORT_TMP/sort \
            -o {output.bam} - \
            2>> minimap.chunk{wildcards.chunk}.log
        """

rule merge_minimap_chunks:
    input:
        bams = expand("Align/tmp/minimap_chunks/chunk_{chunk}.sort.bam", chunk=range(MINIMAP_N_CHUNKS))
    output:
        bam = "keep/Align/{}.sort.bam".format(config['samples']['id']),
        bai = "keep/Align/{}.sort.bam.bai".format(config['samples']['id']),
    threads:
        config['threads'].get('minimap_merge', 8)
    benchmark:
        "Benchmarks/main.merge_minimap_chunks.txt"
    shell:
        """
        set -euo pipefail
        samtools merge -@ {threads} -f {output.bam} {input.bams}
        samtools index -@ {threads} {output.bam} {output.bai}
        """
