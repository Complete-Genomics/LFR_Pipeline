# generate bam for hapcut

# split the bam by chromosome
rule split_hapcut_bam:
    input:
        "keep/Align/{id}.sort.removedup_rm000.bam"
    output:
        "Make_Vcf/step3_hapcut/step1_modify_bam/{id}_sort.markdup_{chr}.bam"
    params:
        id = config['samples']['id']
    shell:
        "samtools view -bh {input} {wildcards.chr} > "
            "{output} && "
            "samtools index {output}"

# gzip the pass vars vcf
rule gzip_hapcut_vcf:
    input:
        "Make_Vcf/step1_haplotyper/{id}_gatk_pass_vars.vcf"
    output:
        gz_vcf = "keep/Make_Vcf/step3_hapcut/step2_split_vcf/{id}_gatk_pass_vars.vcf.gz",
        index = "keep/Make_Vcf/step3_hapcut/step2_split_vcf/{id}_gatk_pass_vars.vcf.gz.tbi"
    params:
        bcftools = config['params']['bcftools'],
        tabix = config['params']['tabix']
    shell:
        "{params.bcftools} view -O z {input} > "
            "{output.gz_vcf} && "
            "{params.tabix} -p vcf -f {output.gz_vcf}"


# split the gzipped vcf by chromosome
rule split_hapcut_vcf:
    input:
        vcf = "keep/Make_Vcf/step3_hapcut/step2_split_vcf/{id}_gatk_pass_vars.vcf.gz",
        index = "keep/Make_Vcf/step3_hapcut/step2_split_vcf/{id}_gatk_pass_vars.vcf.gz.tbi"
    output:
        "Make_Vcf/step3_hapcut/step2_split_vcf/{id}_gatk_pass_vars_{chr}.vcf"
    params:
        bcftools = config['params']['bcftools'],
        id = config['samples']['id']
    shell:
        "{params.bcftools} view -O v "
            "-r {wildcards.chr} "
            "{input.vcf} > "
            "{output}"


# run extract hairs to get fragments for hapcut
rule get_hapcut_fragments:
    input:
        bam = "Make_Vcf/step3_hapcut/step1_modify_bam/{id}_sort.markdup_{chr}.bam",
        vcf = "Make_Vcf/step3_hapcut/step2_split_vcf/{id}_gatk_pass_vars_{chr}.vcf"
    output:
        "Make_Vcf/step3_hapcut/step3_run_hapcut2_10xpipeline/s1_unlinked_frag/{id}_unlinked_fragment_{chr}"
    params:
        toolsdir = config['params']['toolsdir'],
        HapCUT2 = config['params']['hapcut2']
    shell:
        "{params.HapCUT2}/extractHAIRS --10X 1 --indels 1 "
            "--bam {input.bam} "
            "--VCF {input.vcf} "
            "--out {output}"


# run LinkFragments.py
# This preps files for hapcut
rule link_hapcut_fragments:
    input:
        bam = "Make_Vcf/step3_hapcut/step1_modify_bam/{id}_sort.markdup_{chr}.bam",
        vcf = "Make_Vcf/step3_hapcut/step2_split_vcf/{id}_gatk_pass_vars_{chr}.vcf",
        frag = "Make_Vcf/step3_hapcut/step3_run_hapcut2_10xpipeline/s1_unlinked_frag/{id}_unlinked_fragment_{chr}"
    output:
        "Make_Vcf/step3_hapcut/step3_run_hapcut2_10xpipeline/s2_link_frag_files/{id}_linked_frag_{chr}"
    params:
        toolsdir = config['params']['toolsdir'],
        linkdist = config['params']['hapcut_link_dist'],
        HapCUT2 = config['params']['hapcut2']
    shell:
        "/usr/bin/python3 {params.HapCUT2}/LinkFragments.py --bam {input.bam} "
            "--VCF {input.vcf} "
            "--fragments {input.frag} "
            "--out {output} "
            "-d {params.linkdist}"


# Run hapcut for each chromosome
rule run_hapcut2:
    input:
        frag = "Make_Vcf/step3_hapcut/step3_run_hapcut2_10xpipeline/s2_link_frag_files/{id}_linked_frag_{chr}",
        vcf = "Make_Vcf/step3_hapcut/step2_split_vcf/{id}_gatk_pass_vars_{chr}.vcf"
    output:
        blocks = "Make_Vcf/step3_hapcut/step3_run_hapcut2_10xpipeline/s3_hapcut_output/{id}_hapblock_{chr}",
        vcf = "Make_Vcf/step3_hapcut/step3_run_hapcut2_10xpipeline/s3_hapcut_output/{id}_hapblock_{chr}.phased.VCF"
    params:
        toolsdir = config['params']['toolsdir'],
        HapCUT2 = config['params']['hapcut2'],
        # prune = config['params']['prune'],
    shell:
        "{params.HapCUT2}/HAPCUT2 --nf 1 "
            "--fragments {input.frag} "
            "--VCF {input.vcf} "
            "--output {output.blocks}"
            "--threshold {params.prune} "


# evaluate hapcut
rule evaluate_hapcut2:
    input:
        hapblock = expand("Make_Vcf/step3_hapcut/step3_run_hapcut2_10xpipeline/s3_hapcut_output/{id}_hapblock_{chr}", id=config['samples']['id'], chr=CHROMS),
        frag = expand("Make_Vcf/step3_hapcut/step3_run_hapcut2_10xpipeline/s2_link_frag_files/{id}_linked_frag_{chr}", id=config['samples']['id'], chr=CHROMS),
        vcf = expand("Make_Vcf/step3_hapcut/step2_split_vcf/{id}_gatk_pass_vars_{chr}.vcf", id=config['samples']['id'], chr=CHROMS)
    output:
        "Make_Vcf/step3_hapcut/step4_compare_with_refphasing/hapcut_eval.txt"
    params:
        toolsdir = config['params']['toolsdir'],
        src_dir = config['params']['src_dir'],
        truth_chroms = config['benchmark']['truth_vcf_dir'],
        ref_index = REF + '.fai',
        ref_id = config['params']['ref_id']
    run:
        import sys
        from pathlib import Path

        # if there's a truth_chroms parameter, we should be evaluating a benchmarked sample
        # the truth_chroms_{chr}.vcf should contain the phased truthset VCF for that chromosome
        if params.ref_id in ['hg002', 'hg001', 'hg005']:
            truth_chroms = expand(params.truth_chroms + "/truth_chroms_{chr}.vcf", chr=CHROMS)
            shell("""python3 {params.src_dir}/modules/shared/phasing/calculate_haplotype_statistics.py --haplotype_blocks {input.hapblock} --vcf {input.vcf} --fragments {input.frag} --phased_vcf {truth_chroms} --contig_size_file {params.ref_index} > {output}""")
        # Otherwise we import the hapblock_vcf_error_rate_multiple module from hapcut
        # and then use our vcf as the truth VCF
        # This means we won't get any switch statistics, but without knowing the truth
        # we can't figure that out anyway
        # This does allow us to get the number of phased SNPs, N50 and NA50
        else:
            sys.path.append(str(Path(params.toolsdir) / 'tools'))
            from calculate_haplotype_statistics_no_truth import hapblock_vcf_error_rate_multiple

            err = hapblock_vcf_error_rate_multiple(input.hapblock, input.vcf, input.vcf, False)
            with open(output[0], "w") as outf:
                print(err, file=outf)

rule PSblock2bed:
    input:
        "Make_Vcf/step3_hapcut/step4_compare_with_refphasing/hapcut_eval.txt"
    params:
        python = config['params']['general_python'],
        src_dir = config['params']['src_dir']
    output:
        "Make_Vcf/step3_hapcut/step3_run_hapcut2_10xpipeline/s3_hapcut_output/data_hapblock.bed"
    run:
        command = ["{params.python}",
                   "{params.src_dir}/modules/shared/phasing/hapcutPS2bed.py --input_dir Make_Vcf/step3_hapcut/step3_run_hapcut2_10xpipeline/s3_hapcut_output "] 
        shell(" ".join(command))

# Aggregate all of the hapcut VCFs as a full phased VCF
rule aggregate_hapcut_vcf:
    input:
        vcfs = expand("Make_Vcf/step3_hapcut/step3_run_hapcut2_10xpipeline/s3_hapcut_output/{id}_hapblock_{chr}.phased.VCF", id=config['samples']['id'], chr=CHROMS)
    output:
        agg_vcf = "Make_Vcf/step3_hapcut/step4_compare_with_refphasing/{}_hapcut.phased.vcf".format(config['samples']['id'])
    run:
        import sys, re, glob
        from pathlib import Path

        contigs = []
        # open up the first vcf to get the header and output it to the aggregte vcf
        with open(input.vcfs[0], "r") as vcf_header, open(output.agg_vcf, "w") as outfile:
            for line in vcf_header:
                if line.startswith("#"):
                    print(line.strip(), file=outfile)
                    # scrape the header for contigs present
                    if "contig=" in line:
                        # capture contig info
                        id_name = re.split("<|>", line.strip())[1]
                        # split into a dict of ID: contig name, length: <len>, assembly: ref
                        keys = dict(parts.split("=") for parts in id_name.split(","))
                        contigs.append(keys['ID'])
                else:
                    break


        file_dict = {}
        # iterate through input VCFs
        # to create a dictionary of contigs and file paths
        for file in input.vcfs:
            name = Path(file).name
            tail = name.split("_hapblock_")[-1]
            contig = tail.split(".phased.VCF")[0]
            file_dict[contig] = file

        # iterate through contigs, this means the results will be in order
        for contig in contigs:
            # some contigs won't be present as files if they didn't have SNPs
            # this is more relevant for alignments to de novo assemblies
            try:
                file = file_dict[contig]
            except KeyError:
                print(f"No file for contig {contig}, skipping...", file=sys.stderr)
                continue
            # if we find the file, open it and output results to the aggregate file
            with open(file, "r") as vcf, open(output.agg_vcf, "a") as outfile:
                for line in vcf:
                    if line.startswith("#"):
                        continue
                    else:
                        print(line.strip(), file=outfile)
