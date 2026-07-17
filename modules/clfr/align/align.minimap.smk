# clfr, t1 se600, align with minimap
import os

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
        config['threads']['bwa']
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
        readgroup = r'@RG\tID:{0}\tSM:{0}\tPL:{1}'.format(config['samples']['id'],
                                                          config['params']['platform'])
    benchmark:
        "Benchmarks/main.map_reads.txt"
    shell:
        """
        # ============ 1. set env ============
        #log=Align/minimap2.log
        #tmp_prefix=/dev/shm/minimap_tmp
        
        SORT_TMP=/dev/shm/minimap_tmp_$$
        echo "minimap2 samtools sort tmp: $SORT_TMP" >&2
        mkdir -p $SORT_TMP


        # ============ 2. minimap2 mapping ============
        {params.MINIMAP} -ax splice:sr \
            -t 15 \
            --secondary=no \
            --sam-hit-only \
            --junc-bed {params.anno_bed} \
            {input.ref} {input.fq2} \
            2>> minimap.log \
        | samtools view -@ 4 -b - \
        | samtools sort -@ 10 -m 2G \
            -T $SORT_TMP/sort \
            -o {output.bam} - \
            2>> minimap.log

        # ============ 3. index ============
        samtools index -@ 8 {output.bam} {output.bai} 2>> minimap.log

        # ============ 4. clean up ============
        rm -rf $SORT_TMP

        """


