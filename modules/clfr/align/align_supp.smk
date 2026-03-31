import os

rule skip_mark_duplicates:
    input:
        "keep/Align/{id}.sort.bam",
    output:
        tmp = "keep/Align/{id}.sort.markdup.bam",
        bam="Align/{id}.sort.removedup_rm000.bam",
        bai="Align/{id}.sort.removedup_rm000.bam.bai"
    # benchmark:
    #     "Benchmarks/calc_frag_len.remove_duplicates.{id}.txt"
    threads:
        config['threads']['bwa']
    run:
        # if config['params']['bc_condition'] == 'random_bc' or config['params']['bc_condition'] == 'random_bc_umi_rc' or config['params']['bc_condition'] == 'pcrfree':
        shell("cd keep/Align && ln -s data.sort.bam data.sort.markdup.bam && cd ../../Align && ln -s ../keep/Align/data.sort.bam data.sort.removedup_rm000.bam && cd .. && samtools index {output.bam}")        # elif 'standard' in config['params']['bc_condition']:
        #     shell("samtools view -b -h -F 0x400 {input}  > {output.bam} && samtools index {output.bam}")
        

# idxstats: Primary Alignment
rule mapped_uniq_bc_bases_count:
    input:
        bam="Align/{}.sort.removedup_rm000.bam".format(config['samples']['id']),
        bai="Align/{}.sort.removedup_rm000.bam.bai".format(config['samples']['id'])
    output:
        "Align/idxstats_removedup_rm000.bam.txt",
    shell:
        "samtools idxstats {input.bam} > {output} "

rule count_idxstats:
    input:
        "Align/idxstats_removedup_rm000.bam.txt",
        "Align/picard_align_metrics.txt"
        # "Align/picard_align_metrics_BCreads_only.txt"
    output:
        "Align/mapped_uniq_bc_bases_count.txt"
    params:
        src=config['params']['src_dir'],
        my_python= config['params']['general_python'],
    shell:
        "{params.my_python} {params.src}/modules/stlfr/align/lfr_supp.py --module idxstats && cat {output} >> Align/picard_align_metrics.txt"