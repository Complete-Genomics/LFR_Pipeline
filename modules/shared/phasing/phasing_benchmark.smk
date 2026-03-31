# phasing use benchmark vcf

vcf_benchmark_dir = config['benchmark']['vcf_benchmark'] 

rule get_hapcut_fragments_benchmark:
    input:
        bam = "Make_Vcf/step3_hapcut/step1_modify_bam/{id}_sort.markdup_{chr}.bam",
        vcf = lambda wildcards: f"{vcf_benchmark_dir}/truth_chroms_{wildcards.chr}.vcf",

    output:
        "Benchmarks/step3_hapcut/step3_run_hapcut2_10xpipeline/s1_unlinked_frag/{id}_unlinked_fragment_{chr}"
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
rule link_hapcut_fragments_benchmark:
    input:
        bam = "Make_Vcf/step3_hapcut/step1_modify_bam/{id}_sort.markdup_{chr}.bam",
        vcf = lambda wildcards: f"{vcf_benchmark_dir}/truth_chroms_{wildcards.chr}.vcf",
        frag = "Benchmarks/step3_hapcut/step3_run_hapcut2_10xpipeline/s1_unlinked_frag/{id}_unlinked_fragment_{chr}"
    output:
        "Benchmarks/step3_hapcut/step3_run_hapcut2_10xpipeline/s2_link_frag_files/{id}_linked_frag_{chr}"
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
rule run_hapcut2_benchmark:
    input:
        frag = "Benchmarks/step3_hapcut/step3_run_hapcut2_10xpipeline/s2_link_frag_files/{id}_linked_frag_{chr}",
        vcf = lambda wildcards: f"{vcf_benchmark_dir}/truth_chroms_{wildcards.chr}.vcf",
    output:
        blocks = "Benchmarks/step3_hapcut/step3_run_hapcut2_10xpipeline/s3_hapcut_output/{id}_hapblock_{chr}",
        vcf = "Benchmarks/step3_hapcut/step3_run_hapcut2_10xpipeline/s3_hapcut_output/{id}_hapblock_{chr}.phased.VCF"
    params:
        toolsdir = config['params']['toolsdir'],
        HapCUT2 = config['params']['hapcut2']
    shell:
        "{params.HapCUT2}/HAPCUT2 --nf 1 "
            "--fragments {input.frag} "
            "--VCF {input.vcf} "
            "--output {output.blocks}"


# evaluate hapcut
rule evaluate_hapcut2_benchmark:
    input:
        hapblock = expand("Benchmarks/step3_hapcut/step3_run_hapcut2_10xpipeline/s3_hapcut_output/{id}_hapblock_{chr}", id=config['samples']['id'], chr=CHROMS),
        frag = expand("Benchmarks/step3_hapcut/step3_run_hapcut2_10xpipeline/s2_link_frag_files/{id}_linked_frag_{chr}", id=config['samples']['id'], chr=CHROMS),
        vcf = expand("{vcf_benchmark_dir}/truth_chroms_{chr}.vcf", vcf_benchmark_dir=vcf_benchmark_dir, chr=CHROMS)
    output:
        "Benchmarks/step3_hapcut/step4_compare_with_refphasing/hapcut_eval.txt"
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
        # if params.ref_id in ['hg002', 'hg001', 'hg005']:
        # We import the hapblock_vcf_error_rate_multiple module from hapcut
        # and then use our vcf as the truth VCF
        # This means we won't get any switch statistics, but without knowing the truth
        # we can't figure that out anyway
        # This does allow us to get the number of phased SNPs, N50 and NA50
        sys.path.append(str(Path(params.toolsdir) / 'tools'))
        from calculate_haplotype_statistics_no_truth import hapblock_vcf_error_rate_multiple

        err = hapblock_vcf_error_rate_multiple(input.hapblock, input.vcf, input.vcf, False)
        with open(output[0], "w") as outf:
            print(err, file=outf)

# rule PSblock2bed_benchmark:
#     input:
#         "Benchmarks/step3_hapcut/step4_compare_with_refphasing/hapcut_eval.txt"
#     params:
#         python = config['params']['general_python'],
#         src_dir = config['params']['src_dir']
#     output:
#         "Benchmarks/step3_hapcut/step3_run_hapcut2_10xpipeline/s3_hapcut_output/data_hapblock.bed"
#     run:
#         command = ["{params.python}",
#                    "{params.src_dir}/src/hapcutPS2bed.py"] 
#         shell(" ".join(command))

# Aggregate all of the hapcut VCFs as a full phased VCF
rule aggregate_hapcut_vcf_benchmark:
    input:
        vcfs = expand("Benchmarks/step3_hapcut/step3_run_hapcut2_10xpipeline/s3_hapcut_output/{id}_hapblock_{chr}.phased.VCF", id=config['samples']['id'], chr=CHROMS)
    output:
        agg_vcf = "Benchmarks/step3_hapcut/step4_compare_with_refphasing/{}_hapcut.phased.vcf".format(config['samples']['id'])
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


