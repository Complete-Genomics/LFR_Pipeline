rule merge2bam:
    input:
        r1="Align/r1.data.sort.bam", 
        r2="Align/r2.data.sort.bam"
    output:
        b1="Align/data.sort.bam",
        b2="Align/data.sort.markdup.bam",
    params: 
        samtools = config['params']['samtools']
    shell:
        "{params.samtools} merge -@ 20 -o {output.b1} {input.r1} {input.r2} cd && Align && ln -s {output.b1} {output.b2} && cd .. "        


# Index markdups bam
rule index_mark_dups:
    input:
        bam = "Align/{id}.sort.markdup.bam"
    output:
        bai = "Align/{id}.sort.markdup.bam.bai"
    threads:
        config['threads']['bwa']
    shell:
        "samtools index -@ {threads} {input}"

# the mark duplicates rule just addds the flag, the below removes duplicates
# hapcut uses this sam file
# We also remove barcode uninformative reads (0_0_0)
rule remove_duplicates:
    input:
        "Align/{id}.sort.markdup.bam"
    output:
        "Align/{id}.sort.removedup_rm000.sam"
    benchmark:
        "Benchmarks/calc_frag_len.remove_duplicates.{id}.txt"
    threads:
        config['threads']['bwa']
    run:
        if config['params']['bc_condition'] == 'random_bc' or config['params']['bc_condition'] == 'random_bc_umi_rc' or config['params']['bc_condition'] == 'pcrfree':
            shell("samtools view {input} > {output}")
        elif config['params']['bc_condition'] == 'standard':
            shell("samtools view -h -F 0x400 {input} | "
            "awk -F $'\t' '($1!~/#0_0_0$/){{print}}' > {output}")
        

rule sam2bam_removedup:
    input:
        bam = "Align/{id}.sort.markdup.bam"
    output:
        "Align/{id}.sort.removedup_rm000.bam"
    run:
        if 'random_bc' in config['params']['bc_condition']:
            shell("cd Align && ln -s data.sort.markdup.bam data.sort.removedup_rm000.bam && cd .. && samtools index {output}")
        elif config['params']['bc_condition'] == 'standard':
            shell("samtools view -bh {input.sam} > {output} && samtools index {output}")


rule remove_000:
    input:
        "Align/{id}.sort.markdup.bam"
    output:
        "Align/{id}.sort.markdup_rm000.sam"
    run:
        if 'random_bc' in config['params']['bc_condition']:
            shell("touch  {output}")
        elif config['params']['bc_condition'] == 'standard':
            shell("samtools view -h {input} | "
            "awk -F $'\t' '($1!~/#0_0_0$/){{print}}' > {output}")

rule sam2bam_markdup:
    input:
        "Align/{id}.sort.markdup_rm000.sam"
    output:
        "Align/{id}.sort.markdup_rm000.bam"
    run:
        if 'random_bc' in config['params']['bc_condition']:
                shell("touch {output}")
        elif config['params']['bc_condition'] == 'standard':
            shell("samtools view -bh {input} > {output} && samtools index {output}")

# def summary_report_input(wildcards):
#     summary_report_files = []
#     # if config['params']['bc_condition'] == 'random_bc':
#     #     summary_report_files.append("Align/picard_align_metrics_bc.txt")
#     # if config['params']['bc_condition'] == 'standard':
#     #     summary_report_files.append("Align/picard_align_metrics_BCreads_only.txt")
#     # Add stLFR specific metrics that the summary report uses
#     if config['modules']['LFR']:
#         calc_frag_file = ["Calc_Frag_Length/frag_length_distribution.pdf",
#                           "Calc_Frag_Length/n_read_distribution.pdf",
#                           "Calc_Frag_Length/frag_and_bc_summary.txt",
#                           "Calc_Frag_Length/frags_per_bc.pdf",
#                           "Calc_Frag_Length/frag_summary_minreads"+str(MINREADS)+".txt"
#                           ]


#         for split_dist in config['calc_frag']['split_dist']:
#             for outfile in calc_frag_file:
#                 parts = outfile.split("/")
#                 summary_report_files.append(parts[0] + "_" + str(split_dist) + "/" + parts[1])

#     # Add phasing specific metrics that the summary report uses
#     if config['modules']['phasing']:
#         # summary_report_files.append("Make_Vcf/step4_longhap/longhap_results.txt")
#         summary_report_files.append("Make_Vcf/step3_hapcut/step4_compare_with_refphasing/hapcut_eval.txt")

#     return summary_report_files


# rule to generate the summary report
# This report has pretty much everything that goes into the online spreadsheet
# It will search for all the potential files
# If it can't find a file it'll just emit a message to stderr and keep going
# rule generate_summary_report:
#     input:
#         summary_report_input
#     output:
#         "summary_report.txt"
#     params:
#         toolsdir = config['params']['toolsdir'],
#         samp = config['samples']['id'],
#         read_length = config['params']['read_len'],
#         min_frag = config['calc_frag']['min_frag']
#     benchmark:
#         "Benchmarks/metrics.generate_summary_report.txt"
#     shell:
#         "python3 {params.toolsdir}/tools/summary_report_v5.py | "
#         "tee > {output}"