# Run HaplotypeCaller to generate a vcf without filtering
rule run_variant_calling:
    input:
        bam = "keep/Align/{id}.sort.markdup.bam",
        bai = "keep/Align/{id}.sort.markdup.bam.bai",
        ref = config['params']['ref_fa']
    output:
        "keep/Make_Vcf/step1_haplotyper/{id}_gatk.vcf"
    params:
        gatk_install = config['params']['gatk_install'],
        dbsnp = config['params']['dbsnp_path']
    benchmark:
        "Benchmarks/make_vcf.run_hpcaller.{id}.txt"
    run:
        command = ["{params.gatk_install}", "HaplotypeCaller",
                   "-R", "{input.ref}", "-I", "{input.bam}",
                   "-O", "{output}"]
        # if dbsnp isn't supplied, omit the -d flag
        if params.dbsnp:
            command.extend(["-D", "{params.dbsnp}"])
        shell(" ".join(command))

# Filter to keep pass vars
# This is necessary for LongHap and HapCut since they don't  check the filter column
# rule keep_pass_vars:
#     input:
#         "Make_Vcf/step1_haplotyper/{id}_gatk.vcf"
#     output:
#         "Make_Vcf/step1_haplotyper/{id}_gatk_pass_vars.vcf"
#     benchmark:
#         "Benchmarks/make_vcf.keep_pass_vars.{id}.txt"
#     shell:
#         """
#         awk '($1~/^#/ || $7=="PASS" || $7=="."){{print}}' {input} > {output}
#         """
rule keep_pass_vars:
    input:
        "keep/Make_Vcf/step1_haplotyper/{id}_gatk.vcf"
    params: 
        g=11,
        G=61,
        m=0.11,
        M=0.265,
        x=5.5,
        X=3.95,
        id = config['samples']['id'],
        toolsdir = config['params']['toolsdir'],
        general_python = config['params']['general_python']
    output:
        "Make_Vcf/step1_haplotyper/{id}_gatk_pass_vars.vcf"
    benchmark:
        "Benchmarks/make_vcf.keep_pass_vars.{id}.txt"
    shell:
        """
        {params.general_python} {params.toolsdir}/tools/vcffilter.py  -g {params.g} -G {params.G} -m {params.m} -M {params.M} -x {params.x} -X {params.X} -infile {input} -sample {params.id}
        """

# Filter to keep SNPs
# This is used for evaluating against a benchmark
rule select_snps:
    input:
        "keep/Make_Vcf/step1_haplotyper/{id}_gatk_pass_vars.vcf"
    output:
        "Make_Vcf/step2_benchmarking/{id}.snp.vcf.gz"
    benchmark:
        "Benchmarks/make_vcf.select_snps.{id}.txt"
    params:
        bcftools = config['params']['bcftools']
    shell:
        "{params.bcftools} view -O z --type snps {input} > {output}"


# Filter to keep InDels
# This is used for evaluating against a benchmark
rule select_indels:
    input:
        "Make_Vcf/step1_haplotyper/{id}_gatk_pass_vars.vcf"
    output:
        "Make_Vcf/step2_benchmarking/{id}.indel.vcf.gz"
    benchmark:
        "Benchmarks/make_vcf.select_indels.{id}.txt"
    params:
        bcftools = config['params']['bcftools']
    shell:
        "{params.bcftools} view -O z --type indels {input} > {output}"


# index the snps
rule index_snps:
    input:
        "Make_Vcf/step2_benchmarking/{id}.snp.vcf.gz"
    output:
        "Make_Vcf/step2_benchmarking/{id}.snp.vcf.gz.tbi"
    benchmark:
        "Benchmarks/make_vcf.index_snps.{id}.txt"
    params:
        tabix = config['params']['tabix']
    shell:
        "{params.tabix} -p vcf -f {input}"


# index the indels
rule index_indels:
    input:
        "Make_Vcf/step2_benchmarking/{id}.indel.vcf.gz"
    output:
        "Make_Vcf/step2_benchmarking/{id}.indel.vcf.gz.tbi"
    benchmark:
        "Benchmarks/make_vcf.index_indels.{id}.txt"
    params:
        tabix = config['params']['tabix']
    shell:
        "{params.tabix} -p vcf -f {input}"


# evaluates snps against the benchmark
rule eval_snps:
    input:
        snp = "Make_Vcf/step2_benchmarking/{}.snp.vcf.gz".format(config['samples']['id']),
        index = "Make_Vcf/step2_benchmarking/{}.snp.vcf.gz.tbi".format(config['samples']['id'])
    output:
        "Make_Vcf/step2_benchmarking/snp_compare/summary.txt"
    params:
        benchmark_snp = config['benchmark']['benchmark_snp'],
        ref_sdf = config['benchmark']['ref_sdf'],
        bedfile = config['benchmark']['bedfile'],
        direc = "Make_Vcf/step2_benchmarking/snp_compare",
        rtg_install = config['params']['rtg_install']
    benchmark:
        "Benchmarks/make_vcf.eval_snps.txt"
    shell:
        # snakemake generates expected directories, but vcfeval won't run if the directory is already generated
        # so we have to remove it
        "rm -r {params.direc};"
        "{params.rtg_install} vcfeval "
            "-b {params.benchmark_snp} "
            "-c {input.snp} "
            "-e {params.bedfile} "
            "-t {params.ref_sdf} "
            "-o {params.direc}"


# evaluate indels against benchmark
rule eval_indels:
    input:
        indel = "Make_Vcf/step2_benchmarking/{}.indel.vcf.gz".format(config['samples']['id']),
        index = "Make_Vcf/step2_benchmarking/{}.indel.vcf.gz.tbi".format(config['samples']['id'])
    output:
        "Make_Vcf/step2_benchmarking/indel_compare/summary.txt"
    params:
        benchmark_indel = config['benchmark']['benchmark_indel'],
        ref_sdf = config['benchmark']['ref_sdf'],
        bedfile = config['benchmark']['bedfile'],
        direc = "Make_Vcf/step2_benchmarking/indel_compare",
        rtg_install = config['params']['rtg_install']
    benchmark:
        "Benchmarks/make_vcf.eval_indels.txt"
    shell:
        # snakemake generates expected directories, but vcfeval won't run if the directory is already generated
        # so we have to remove it
        "rm -r {params.direc};"
        "{params.rtg_install} vcfeval "
            "-b {params.benchmark_indel} "
            "-c {input.indel} "
            "-e {params.bedfile} "
            "-t {params.ref_sdf} "
            "-o {params.direc} "
