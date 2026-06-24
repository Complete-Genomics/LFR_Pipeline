# this formats the chroms to be used as a flag for calc_frag_len.py
def get_chroms(wildcards):
    if len(CHROMS) > 1:
        return ",".join(CHROMS)
    else:
        assert len(CHROMS) == 1
        return CHROMS


# calculate fragment lengths, intial calc_frag_len.py remove reads by same position to dedup, now remove dup in BC-dependt way ahead and no longer deup in calc_frag_len.py
rule calc_frag_len:
    input:
        bam = "Align/{}.sort.removedup_rm000.bam".format(config['samples']['id']),
    output:
        "Calc_Frag_Length_{split}/frag_length_distribution.pdf",
        "Calc_Frag_Length_{split}/n_read_distribution.pdf",
        "Calc_Frag_Length_{split}/frag_and_bc_summary.txt",
        "Calc_Frag_Length_{split}/frag_and_bc_dataframe.tsv",
        "Calc_Frag_Length_{split}/frag_summary_minreads"+str(MINREADS)+".txt"

    params:
        min_frag = MIN_FRAG,
        read_len = config['params']['read_len'],
        chroms = get_chroms,
        toolsdir = config['params']['toolsdir'],
        include_dups = config['calc_frag']['include_dups'],
        minreads = MINREADS,
        ray_mem = config['calc_frag']['ray_mem'],
        src_dir = config['params']['src_dir'],
        writeouttsvs = config['calc_frag']['writeouttsvs'],
        mapping_quality = MAPPING_QUALITY,
        python = config['params']['calc_frag_python']
    threads:
        config['threads']['calc_frag']
    benchmark:
        "Benchmarks/calc_frag_len.calc_frag_len_{split}.txt"
    run:
        command = ["{params.python}",
                   "{params.src_dir}/src/calc_frag_len_optimized.py",
                   "--minfrag {params.min_frag}",
                   "--splitdist {wildcards.split}",
                   "--readlen {params.read_len}",
                   "--lfr_type stlfr",
                   "--threads {threads}",
                   "--minreads {params.minreads}",
                   "--ray_mem {params.ray_mem} ",
                   "--mapping_quality {params.mapping_quality} ",
                   "--outdir Calc_Frag_Length_{wildcards.split}/"] 
        if len(params.chroms)>0:
            command.append("--chroms {params.chroms}")                               
        if params.include_dups:
            command.append("--includedups")
        if params.writeouttsvs:
            command.append("--writeouttsvs")
        command.append("--bampath {input.bam}")

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
