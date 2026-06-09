## TODO: as it's full length of reads, need to refactor picard.py to measure mapping rate

rule split_reads_BCgDNA:
    input:
        expand("data/read_{i}.fq.gz", i=range(1,3))
    output:
        out1 =expand("data/split_read_{bc_type}.2.fq.gz", bc_type=["nobc1", "nobc2bc3", "noall3", "nobc3"])
        # "split_stat_read1.log"
    params:
        r2_len = config['params']['read_len'],
        r1_len = config['params']['read_len_r1'], 
        toolsdir = config['params']['toolsdir'],
        src_dir = config['params']['src_dir'],
        randomBC_dir = config['params']['randomBC_dir'],
        general_python = config['params']['general_python'],
        bc_len = config['params']['sbc_len'],
        bc_len_redundant = config['params']['cbc_len'],
        cbc_len = config['params']['cbc_len'],
        bc_start = config['params']['bc_start'],
        gdna_start = config['params']['gdna_start'],
        additional_bc_start = config['params']['additional_bc_start'],
        additional_bc_len = config['params']['additional_bc_len'],
        gdna_start_r1 = config['params']['gdna_start_r1'], 
        additional_bc_len_r1 = config['params']['additional_bc_len_r1'],
        adapter_len = config['params']['adapter_len'],
    threads:
        config['threads']['split_group']
    run:
        if config['params']['bc_condition'] == 'BCgDNA':
            params.barcode = determine_barcode_list(400000)
            shell("{params.general_python} {params.src_dir}/modules/shared/splitreads/split_barcode_stLFR.py "
                    "--barcode {params.barcode} "
                    "--r1 {input[0]} --r2 {input[1]} --read_len {params.r2_len} --output data/split_read "
                    "--bc_start {params.bc_start} --gdna_start {params.gdna_start} "
                    "--additional_bc_start {params.additional_bc_start} --additional_bc_len {params.additional_bc_len} "
                    "--gdna_start_r1 {params.gdna_start_r1} --read_len_r1 {params.r1_len} "
                    "--adapter_len {params.adapter_len} --swap none --output_mode stratified --threads {threads} "
                    "2> data/split_stat_read.err")

rule map_gDNA:
    input:
        fq= "data/split_read_{bc_type}.2.fq.gz"
    output:
        bam="Align/{bc_type}.sort.bam",
        sam="Align/{bc_type}.sort.sam"
    params:
        ref = REF,
        bwa_mem = config['params']['bwa_mem'],
        gatk_install = config['params']['gatk_install'],
        readgroup = r'@RG\tID:{0}\tSM:{0}\tPL:{1}'.format(config['samples']['id'],
                                                          config['params']['platform'])
    run:
        if SEQUENCE_TYPE=='pe':
            shell("bwa {params.bwa_mem} -M -R '{params.readgroup}' -C "
                "-t {threads} {params.ref} data/split_read_{wildcards.bc_type}.1.fq.gz data/split_read_{wildcards.bc_type}.2.fq.gz 2>Align/aln.err | tee {output.sam} | "
            "samtools sort -o {output.bam} -@ {threads} -O bam -")
        elif SEQUENCE_TYPE=='se':
            shell("bwa {params.bwa_mem} -M -R '{params.readgroup}' -C "
                "-t {threads} {params.ref} data/split_read_{wildcards.bc_type}.2.fq.gz 2>Align/aln.err | tee {output.sam} | "
            "samtools sort -o {output.bam} -@ {threads} -O bam -") 


rule filter_r2:
    input:
        "Align/{bc_type}.sort.bam"
    output:
        "Align/{bc_type}.sort.r2.bam.txt"
    run:
        shell("samtools view -f 128 {input} | head -200 > {output} --quiet || true ")

# rule filter_cigar:
#     input:
#         "Align/{bc_type}.sort.r2.bam"
#     output:
#         "Align/{bc_type}.sort.r2.bam.txt"
#     run:
#         shell("samtools view {input} ")
                                                  

rule filter_pos_tag:
    input:
        "Align/{bc_type}.sort.r2.bam.txt"
    output:
        "Align/{bc_type}.sort.r2.bam.pos.txt"
    params:
        python = config['params']['general_python'],
        src_dir = config['params']['src_dir']
    run:
        shell("{params.python} {params.src_dir}/modules/stlfr/bcgdna/BCgDNA_filterCIGAR.py --bam_r2_type {wildcards.bc_type} > {output} ")
