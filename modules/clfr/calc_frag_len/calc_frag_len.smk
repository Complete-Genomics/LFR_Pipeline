rule slice_chr:
    input:
        bam = "Align/data.sort.markdup.bam",
        bai = "Align/data.sort.markdup.bam.bai"
    output:
        bam = "Align/chr19.dedup.bam",
        bai = "Align/chr19.dedup.bam.bai"
    shell:
        "samtools view -b {input.bam} chr19 > {output.bam} && samtools index {output.bam}"


# this formats the chroms to be used as a flag for calc_frag_len.py
def get_chroms(wildcards):
    if len(CHROMS) > 1:
        return ",".join(CHROMS)
    else:
        return CHROMS

# calculate fragment lengths
rule calc_frag_len:
    input:
        bam = "Align/{}.sort.removedup_rm000.bam".format(config['samples']['id']),
    output:
        "Calc_Frag_Length_{split}/frag_length_distribution.pdf",
        "Calc_Frag_Length_{split}/frag_and_bc_summary.txt",
        "Calc_Frag_Length_{split}/frag_minreads"+str(MINREADS)+"_merged.bed",
        "Calc_Frag_Length_{split}/frag_and_bc_dataframe.tsv",
        "Calc_Frag_Length_{split}/barcode_collection",
        "Calc_Frag_Length_{split}/frag_summary_minreads"+str(MINREADS)+".txt"
    params:
        min_frag = MIN_FRAG,
        read_len = config['params']['read_len'],
        chroms = get_chroms,
        toolsdir = config['params']['toolsdir'],
        src_dir = config['params']['src_dir'],
        include_dups = config['calc_frag']['include_dups'],
        minreads = MINREADS,
        n_tolerance = config['calc_frag']['n_tolerance'],
        cbc_len = config['params']['cbc_len'],
        sequence_type = SEQUENCE_TYPE,
        writeouttsvs = config['calc_frag']['writeouttsvs'],
        mapping_quality = MAPPING_QUALITY,
        BC_type = config['params']['bc_condition'],
        python = config['params']['general_python'],

    threads:
        config['threads']['calc_frag']
    benchmark:
        "Benchmarks/calc_frag_len.calc_frag_len_{split}.txt"
    run:
        # fix error 'too long input for bash' for cDNA ref
        command = ["{params.python}",
                "{params.src_dir}/src/calc_frag_len_optimized.py",
                "--minfrag {params.min_frag}",
                "--splitdist {wildcards.split}",
                "--readlen {params.read_len}",
                "--lfr_type clfr",
                "--threads {threads}",
                "--minreads {params.minreads}",
                "--n_tolerance {params.n_tolerance}", 
                "--sequence_type {params.sequence_type}",
                "--cbc_len {params.cbc_len}",
                "--mapping_quality {params.mapping_quality}",
                "--BC_type {params.BC_type} ",
                "--outdir Calc_Frag_Length_{wildcards.split}"]
        if len(params.chroms)>0:
            command.append("--chroms {params.chroms}")

        if params.include_dups:
            command.append("--includedups")
        if params.writeouttsvs:
            command.append("--writeouttsvs")
        command.append("--bampath {input.bam}")

        shell(" ".join(command))

rule coverage_bias:
    params:
        num_bin_in_frag = config['calc_frag']['num_bin_in_frag'],
        minreads = MINREADS,
        read_len = config['params']['read_len'],
        python = config['params']['general_python'],
        src_dir = config['params']['src_dir'], 
        library_type = config['params']['library_type']
    input:
        "Calc_Frag_Length_{split}/frag_and_bc_dataframe.tsv",
        "Calc_Frag_Length_{split}/frag_summary_minreads"+f'{MINREADS}'+".txt"
    output:
        "Calc_Frag_Length_{split}/bin_num10/mean_cov_bin_"+f'{MINREADS}'+"_done",
        "Calc_Frag_Length_{split}/bin_size100/mean_cov_bin_"+f'{MINREADS}'+"_done"

    run:
        command = ["{params.python}",
                   "{params.src_dir}/modules/clfr/calc_frag_len/frag_mean_cov.py",
                   "--read_len {params.read_len}",
                   "--minreads {params.minreads}",
                   "--num_bin_in_frag {params.num_bin_in_frag}",
                   "--library_type {params.library_type}",
                   "--dirname Calc_Frag_Length_{wildcards.split}/"] 
        shell(" ".join(command))

rule coverage_bias_boxplot:
    params:
        name_sample = config['samples']['lanes'],
        num_bin_in_frag = config['calc_frag']['num_bin_in_frag'],
        read_len = config['params']['read_len'],
        python = config['params']['general_python'],
        src_dir = config['params']['src_dir'],
        ylim = config['calc_frag']['ylim'],
        minreads = MINREADS,
        library_type = config['params']['library_type']
    input:
        "Calc_Frag_Length_{split}/bin_num10/mean_cov_bin_"+f'{MINREADS}'+"_done",
        "Calc_Frag_Length_{split}/bin_size100/mean_cov_bin_"+f'{MINREADS}'+"_done"
    output:
        "Calc_Frag_Length_{split}/bin_num10/mean_cov_bin_boxplot_3000.png",
        "Calc_Frag_Length_{split}/mean_cov_bin_boxplot_"+f'{MINREADS}'+"_done"

    run:
        command = ["{params.python}",
                   "{params.src_dir}/modules/clfr/calc_frag_len/frag_mean_cov_boxplot.py",
                   "--name_sample {params.name_sample}",
                   "--read_len {params.read_len}",
                   "--num_bin_in_frag {params.num_bin_in_frag}",
                   "--ylim {params.ylim}",
                   "--minreads {params.minreads}",
                   "--library_type {params.library_type}",
                   "--dirname Calc_Frag_Length_{wildcards.split}/"] 
        shell(" ".join(command))

rule reads_per_bc_distribution:
    params:
        src_dir = config['params']['src_dir'],
        python = config['params']['general_python']
    input:
        "Calc_Frag_Length_{split}/frags_reads_per_bc.tsv"
    output:
        "Calc_Frag_Length_{split}/reads_per_bc_distribution_regardless_mapping.tsv"
    run:
        command = ["{params.python}",
                   "{params.src_dir}/modules/clfr/calc_frag_len/reads_per_bc_distribution.py",
                   "--dirname ./"] 
        shell(" ".join(command))

rule ercc_count:
    input:
        "Align/data.sort.removedup_rm000.bam"
    output:
        "Align/ercc_count.txt"
    params:
        ercc_gtf = config['params']['ercc_gtf']
    run:
        if SEQUENCE_TYPE=='pe':
            shell("featureCounts -T 64 -p -a {params.ercc_gtf} -g gene_id -o Align/ercc_count.txt Align/data.sort.markdup.bam && sed -i '1,1d' Align/ercc_count.txt")
        elif SEQUENCE_TYPE=='se':
            shell("featureCounts -T 64 -a {params.ercc_gtf} -g gene_id -o Align/ercc_count.txt Align/data.sort.markdup.bam && sed -i '1,1d' Align/ercc_count.txt")


rule ercc_plot:
    input:
        cnt = "Align/ercc_count.txt",
    output:
        "Align/ercc_count_ratio.txt"
    params:
        src_dir = config['params']['src_dir'],
        python = config['params']['general_python'],
        ercc_ref_gc = config['params']['ercc_ref_gc']
    shell:
        "{params.python} {params.src_dir}/modules/clfr/calc_frag_len/ercc_count.py --ercc_count {input.cnt} --ercc_ref {params.ercc_ref_gc} --output {output}"
