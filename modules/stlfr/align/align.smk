# Rule for mapping reads, ref=trascripton for mrna data
rule map_reads:
    input:
        ref = REF,
        fq1 = 'data/split_read.1.fq.gz',
        fq2 = 'data/split_read.2.fq.gz'
    output:
        bam = "keep/Align/{}.sort.bam".format(config['samples']['id']),
    threads:
        config['threads']['bwa']
    params:
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
    run:
        # map with readgroup and comments appended, this creates a BX tag with the barcode
        # Also includes sorting
        if config['params']['library_type']=='gdna':
            if SEQUENCE_TYPE=='pe':
                shell("/opt/cgi-tools/bin/bwa {params.bwa_mem} -M -R '{params.readgroup}' -C "
                    "-t {threads} {input.ref} {input.fq1} {input.fq2} 2>Align/aln.err | tee Align/data.aln_mem.sam | "
                "samtools sort -o {output.bam} -@ {threads} -O bam -")
            elif SEQUENCE_TYPE=='se':
                shell("/opt/cgi-tools/bin/bwa {params.bwa_mem} -M -R '{params.readgroup}' -C "
                    "-t {threads} {input.ref} data/split_read.2.fq.gz 2>Align/aln.err | tee Align/data.aln_mem.sam | "
                "samtools sort -o {output.bam} -@ {threads} -O bam -")
        elif config['params']['library_type']=='mrna':
            if MRNA_MAPPER=='hisat2':
                if config['params']['sequence_type']=='pe':
                    shell("{params.hisat2} -x {params.hisat2_index} -1 {input.fq1} -2 {input.fq2} --known-splicesite-infile {params.hisat2_splicesites} --dta-cufflinks --score-min L,0,-0.2 --max-intronlen 1000000 --pen-noncansplice 30 --add-chrname --no-softclip --novel-splicesite-outfile Align/novel-splicesite-outfile.txt -p 20 | samtools sort -o {output.bam} -@ 20 -O bam -")
                elif config['params']['sequence_type']=='se':
                    shell("{params.hisat2} -x {params.hisat2_index} -U {input.fq2} --known-splicesite-infile {params.hisat2_splicesites} --dta-cufflinks --score-min L,0,-0.2 --max-intronlen 1000000 --pen-noncansplice 30 --add-chrname --no-softclip --novel-splicesite-outfile Align/novel-splicesite-outfile.txt -p 20 | samtools sort -o {output.bam} -@ 20 -O bam -")
            elif MRNA_MAPPER=='star':
                if config['params']['sequence_type']=='pe':
                    shell("{params.STAR} "
                    "--runThreadN {threads} "
                    "--genomeDir {params.star_index} "
                    "--outSAMtype BAM SortedByCoordinate "
                    "--outFileNamePrefix Align/ "
                    "--outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 0 "
                    "--readFilesCommand zcat --readFilesIn data/split_read.1.fq.gz data/split_read.2.fq.gz && mv Align/Aligned.sortedByCoord.out.bam {output.bam} "
                    )
                elif config['params']['sequence_type']=='se':
                    shell("{params.STAR} "
                    "--runThreadN {threads} "
                    "--genomeDir {params.star_index} "
                    "--outSAMtype BAM SortedByCoordinate "
                    "--outFileNamePrefix Align/ "
                    "--outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 0 "
                    "--readFilesCommand zcat --readFilesIn data/split_read.2.fq.gz && mv Align/Aligned.sortedByCoord.out.bam {output.bam} "
                    )


# Perform the deduplication step and generate metrics
rule mark_dups:
    input:
        bam = "keep/Align/{id}.sort.bam",
    output:
        bam = "keep/Align/{id}.sort.markdup.bam",
        metrics = "Align/{id}_dedup_metrics.txt"
    params:
        gatk_install = config['params']['gatk_install']
    benchmark:
        "Benchmarks/main.mark_dups.{id}.txt"
    run:
        if 'random_bc' in config['params']['bc_condition']:
            shell("{params.gatk_install} MarkDuplicates -I {input.bam} "
            "-O {output.bam} "
            "-BARCODE_TAG BX "
            "-M {output.metrics}")
        elif 'standard' in config['params']['bc_condition']:
            shell("{params.gatk_install} MarkDuplicates -I {input.bam} "
            "-O {output.bam} "
            "-BARCODE_TAG BC "
            "-M {output.metrics}")
        else:
            shell("{params.gatk_install} MarkDuplicates -I {input.bam} "
            "-O {output.bam} "
            "-M {output.metrics}")
        



# Index markdups bam
rule index_mark_dups:
    input:
        bam = "keep/Align/{id}.sort.markdup.bam"
    output:
        bai = "keep/Align/{id}.sort.markdup.bam.bai"
    threads:
        config['threads']['bwa']
    shell:
        "samtools index -@ {threads} {input}"

# the mark duplicates rule just addds the flag, the below removes duplicates, secondary, supplementary
# hapcut uses this sam file
# We also remove barcode uninformative reads (0_0_0)
rule remove_duplicates:
    input:
        "keep/Align/{id}.sort.markdup.bam"
    output:
        bam="Align/{id}.sort.removedup_rm000.bam",
        bai="Align/{id}.sort.removedup_rm000.bam.bai"
    threads:
        config['threads']['bwa']
    run:
        if config['params']['bc_condition'] == 'random_bc' or config['params']['bc_condition'] == 'random_bc_umi_rc' or config['params']['bc_condition'] == 'pcrfree':
            shell("cd Align && ln -s ../keep/Align/data.sort.markdup.bam data.sort.removedup_rm000.bam && cd .. && samtools index {output.bam}")
        elif 'standard' in config['params']['bc_condition']:
            shell("samtools view -b -h -F 0x400 {input}  > {output.bam} && samtools index {output.bam}")


# This step parses the duplicate metrics and creates a more readable summary
rule mark_dups_txt:
    input:
        "Align/{id}_dedup_metrics.txt"
    output:
        "Align/{id}_dedup_metrics2.txt"
    params:
        toolsdir = config['params']['toolsdir']
    benchmark:
        "Benchmarks/main.mark_dups_txt.{id}.txt"
    shell:
        "perl {params.toolsdir}/tools/mark_dups_txt.pl {input} {output}"

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
    output:
        "Align/mapped_uniq_bc_bases_count.txt"
    params:
        src=config['params']['src_dir'],
        my_python= config['params']['general_python'],
    shell:
        "{params.my_python} {params.src}/modules/stlfr/align/lfr_supp.py --module idxstats && cat {output} >> Align/picard_align_metrics.txt"
