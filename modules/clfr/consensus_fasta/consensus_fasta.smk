'''
part of clfr workflow, get (transcriptome) fasta
'''
NUM_SPLITS_CONSENSUS=5
rule get_chr_dict:
    input:
        bam="Align/{}.sort.removedup_rm000.bam".format(config['samples']['id']),
        bai="Align/{}.sort.removedup_rm000.bam.bai".format(config['samples']['id'])
    output:
        "Align/samtools_idx.txt"
    shell:
        "samtools idxstats {input.bam} > {output}"
        
rule reformat_readid:
    input:
        "Make_Vcf/step3_hapcut/step1_modify_bam/{id}_sort.markdup_{chr}.bam"
    output:
        "Align/tmp/{id}_{chr}.name.bam"
    params:
        python = config['params']['general_python'],
        src_dir = config['params']['src_dir'],
        seq_type = SEQUENCE_TYPE, 
        umi_len = config['params']['cbc_len'],
    run:
        command = ["{params.python}",
                    "{params.src_dir}/modules/clfr/consensus_fasta/consensus_fasta_supp.py",
                    "--module reformat_readid",
                    "--seq_type {params.seq_type} ",
                    "--umi_len {params.umi_len} ",
                    "--chr_name {wildcards.chr}"] 
        shell(" ".join(command))

        # "python {params.src_dir}/modules/clfr/consensus_fasta/consensus_fasta_supp.py --module reformat_readid --chr_name {wildcards.chr}"

rule sort_reformated_bam:
    input:
        "Align/tmp/{id}_{chr}.name.bam"
    output:
        "Align/tmp/{id}_{chr}.name.sort.bam"
    shell:
        "samtools sort -@ 2 -n -o {output} {input} "

# samtools sort -n -o Align/${chr}.name.sort.bam Align/${chr}.name.bam


rule get_consensus_fasta:
    input:
        bam="Align/tmp/{id}_{chr}.name.sort.bam", 
        chr_dict = "Align/samtools_idx.txt"     
    output:
        "Align/tmp/{chr}/{id}_{chr}_{split_idx}.fasta"
    threads: 1
    params:
        python = config['params']['general_python'],
        src_dir = config['params']['src_dir'], 
        ref = config['params']['ref_fa'],
        min_reads = config['frag_de_novo']['reads_per_BC']
    run:
        command = ["{params.python}",
                    "{params.src_dir}/modules/clfr/consensus_fasta/consensus_fasta.py",
                    "--bam {input.bam}",
                    "--ref_fasta {params.ref}",
                    "--output_fasta {output}",
                    "--chrom {wildcards.chr}",
                    "--dict_file {input.chr_dict}",
                    "--split_index {wildcards.split_idx}",
                    "--min_reads {params.min_reads}"
                    ] 
        shell(" ".join(command))


split_cnt = list(range(NUM_SPLITS_CONSENSUS))  # [0, 1, 2, ..., 19]

rule merge_consensus_fasta:
    input:
        fa = expand("Align/tmp/{chr}/{id}_{chr}_{i}.fasta", id=config['samples']['id'],chr=CHROMS, i=split_cnt),
    output:
        "Align/consensus/consensus.fasta"
    shell:
        "cat {input.fa} > {output} "

rule fasta_frag_len_distribution_consensus:
    input:
        "Align/consensus/consensus.fasta"
    output:
        "Align/consensus/consensus_frag_length_distribution.pdf",
        # "Align/frag_length_distribution.txt"
    params:
        python = config['params']['general_python'],
        src_dir = config['params']['src_dir'], 
        minreads_fasta = config['calc_frag']['minreads_fasta'], 
    run:
        command = ["{params.python}",
                    "{params.src_dir}/modules/clfr/consensus_fasta/exon2fasta.py",
                    "--fasta {input}",
                    "--outdir Align/consensus/",
                    "--name consensus",
                    "--minreads_fasta {params.minreads_fasta}",
                    "--module metrics_basic"] 
        shell(" ".join(command))

rule map_fasta_consensus:
    input:
        "Align/consensus/consensus.fasta"
    output:
        "Align/consensus/consensus.paf"
    params:
        minimap = config['frag_de_novo']['minimap'],
        refgenome = config['params']['ref_fa_mrna']
    run:
        command = ["{params.minimap} -x asm20 -t 20 ",
                   "{params.refgenome} ",
                   "{input} > {output} "
                   ] 
        shell(" ".join(command))

rule correc_direction_consensus:
    input:
        consensus_fa = "Align/consensus/consensus.fasta",
        consensus_paf = "Align/consensus/consensus.paf"
    output:
        "Align/consensus/consensus.fixRC.fasta"
    params:
        python = config['params']['general_python'],
        flank_end = config['frag_de_novo']['flank_end'],
        adapter_seq = config['frag_de_novo']['adapter_seq'],
        src_dir = config['params']['src_dir']
    run:
        if config['params']['correct_rc']==True:
            command = ["{params.python} ",
                    "{params.src_dir}/modules/clfr/denovo/denovo_supp.py ",
                    "--adapter_seq {params.adapter_seq} ",
                    "--flank_end {params.flank_end} ",
                    "--fasta {input.consensus_fa}",
                    "--module fix_fa_rc ",
                    "--outdir Align/ "
                    ] 
            shell(" ".join(command))
        else:
            shell("cd Align/consensus && ln -s consensus.fasta consensus.fixRC.fasta ")

rule eval_Mandalorion:
    input:
        "Align/consensus/consensus.fixRC.fasta"
    output:
        "Align/tmp/mando/Isoforms.filtered.clean.gtf"
    params:
        gtf = config['params']['mandalorion_gtf'],
        ref = config['params']['mandalorion_ref'],
        src = config['params']['mandalorion'],
        python_exec = config['params']['mandalorion_python'],
        cpus = config['threads']['metrics']
    shell:
        """
        {params.python_exec} {params.src} -g {params.gtf} -G {params.ref} -f {input} -R 1 -t {params.cpus} -p Align/tmp/mando
        """




rule eval_SQANTI3:
    input:
        "Align/tmp/mando/Isoforms.filtered.clean.gtf"
    output:
        "Align/tmp/SQANTI3_QC_output/data_SQANTI3_report.pdf"
    params:
        python_exec = config['params']['sqanti3_python'],
        sqanti3_script = config['params']['sqanti3_script'],
        ref_gtf = config['params']['mandalorion_gtf'],
        ref_fasta = config['params']['mandalorion_ref'],
        cage_peak_bed = config['params']['sqanti3_cage_peak'],
        polyA_motif_list = config['params']['sqanti3_polyA_motif'],
        cpus = config['threads']['metrics']
    shell:
        """
        {params.python_exec} {params.sqanti3_script} \
            --isoforms {input} \
            --refGTF {params.ref_gtf} \
            --refFasta {params.ref_fasta} \
            --polyA_motif_list {params.polyA_motif_list} \
            --CAGE_peak {params.cage_peak_bed} \
            --force_id_ignore \
            --skipORF \
            -o data \
            --cpus {params.cpus} \
            -d Align/tmp/SQANTI3_QC_output \
            --report both
        """

rule report_consensus:
    input:
        "Align/tmp/SQANTI3_QC_output/data_SQANTI3_report.pdf"
    output:
        "Align/consensus/consensus.fixRC_SQANTI3_report.pdf"
    run:
        shell("cp {input} {output} ")

# rule split_hapcut_bam:
#     input:
#         bam="Align/{id}.sort.removedup_rm000.bam",
#         bai="Align/{id}.sort.removedup_rm000.bam.bai"
#     output:
#         "Make_Vcf/step3_hapcut/step1_modify_bam/{id}_sort.markdup_{chr}.bam"
#     params:
#         id = config['samples']['id']
#     shell:
#         "samtools view -bh {input.bam} {wildcards.chr} > "
#             "{output} && "
#             "samtools index {output}"

# rule index_bam_consensus:
#     input:
#         "Align/{id}.sort.removedup_rm000.bam"
#     output:
#         "Align/{id}.sort.removedup_rm000.bam.bai"
#     run:
#         shell("samtools index {input}")


