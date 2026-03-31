# Rule for mapping reads, ref=trascripton for mrna data
rule map_reads:
    input:
        ref = REF,
        # fqs = expand("{sample}", sample=config['samples']['fastq'])
        fq1 = 'data/split_read.1.fq.gz',
        fq2 = 'data/split_read.2.fq.gz'
    output:
        bam = "keep/Align/{}.sort.bam".format(config['samples']['id']), # may want this to be a temp file
        # sam = "Align/{}.aln_mem.sam".format(config['samples']['id'])
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
                    # "--outSAMattrRGline ID:{params.readgroup} "
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


