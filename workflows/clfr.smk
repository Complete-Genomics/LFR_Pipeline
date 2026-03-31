# include the config file
configfile: "config.yaml"
import pysam

SEQUENCE_TYPE =config['params']['sequence_type'].lower()
MRNA_MAPPER = config['params']['mrna_mapper'].lower()
# Define the ref based on the config file
# Sort of acts like a global variable so you don't need to always type the whole thing
if config['params']['library_species']=='human':
    if config['params']['library_type']=='gdna':
        REF = config['params']['ref_fa']
    elif config['params']['library_type']=='mrna':
        REF = config['params']['ref_fa']
else:
    REF = config['params']['ref_fa_other']

if config['params']['bc_condition']=='standard' or config['params']['bc_condition']=='combined':
    if config['calc_frag']['minreads']=="default" and config['calc_frag']['mapping_quality']=="default":
        MINREADS = 4
        MAPPING_QUALITY = 30
    else:
        MINREADS = config['calc_frag']['minreads']
        MAPPING_QUALITY = config['calc_frag']['mapping_quality']
        

elif 'random_bc' in config['params']['bc_condition']:
    if config['calc_frag']['minreads']=="default" and config['calc_frag']['mapping_quality']=="default":
        MINREADS = 100
        MAPPING_QUALITY = 0
        
    else:
        MINREADS = config['calc_frag']['minreads']
        MAPPING_QUALITY = config['calc_frag']['mapping_quality']
 
MIN_FRAG = config['calc_frag']['min_frag']
# This defines the logic for which chromosomes to use
# Check if the chroms are defined
if not config['samples']['chroms']:
    CHROMS = []
else:
    CHROMS = config['samples']['chroms']


# define a function to return target files based on config settings
def run_all_input(wildcards):
    run_all_files = ["data/split_read.2.fq.gz"]

    if config['modules']['mapping']==True:
        run_all_files.extend(["Align/{}.sort.bam".format(config['samples']['id']),])

    if config['modules']['eval_metrics']==True:
        run_all_files.extend([
                        "Align/gatk_metrics_{}.alignment_summary_metrics".format(config['samples']['id'])])
        
        run_all_files.extend(["summary_report.txt"])
        if config['params']['library_type']=='mrna' and any('ERCC' in chrom for chrom in CHROMS):
            run_all_files.append("Align/ercc_count_ratio.txt")

    # if stLFR is set to true add fragment calculation metrics
    if config['modules']['LFR']:
        calc_frag_file = [
                          "Calc_Frag_Length/frag_and_bc_summary.txt",
                          "Calc_Frag_Length/mean_cov_bin_boxplot_"+f'{MINREADS}'+"_done",
                          "Calc_Frag_Length/frag_minreads"+str(MINREADS)+"_merged.bed",
                          "Calc_Frag_Length/frag_summary_minreads"+f'{MINREADS}'+".txt"]


        # add output directories for the various split distances
        for split_dist in config['calc_frag']['split_dist']:
            for outfile in calc_frag_file:
                parts = outfile.split("/")
                run_all_files.append(parts[0] + "_" + str(split_dist) + "/" + parts[1])

    # if variant_calling is true add VCF
    if config['modules']['variant_calling']:
        run_all_files.append("Make_Vcf/step1_haplotyper/{}_gatk.vcf".format(config['samples']['id']))

    # if benchmarking is set add benchmark summaries
    if config['modules']['benchmarking']:
        run_all_files.append("Make_Vcf/step2_benchmarking/snp_compare/summary.txt")
        run_all_files.append("Make_Vcf/step2_benchmarking/indel_compare/summary.txt")

    # if mate_pair_analysis is set add it's output
    if config['modules']['mate_pair_analysis']:
        run_all_files.append("Align/Distance_Between_Reads/mate_pairs_distance.pdf")

    # if read_overlap_analysis is set add it's output
    if config['modules']['read_overlap_analysis']:
        run_all_files.append("Align/Distance_Between_Reads/distance_between_reads.pdf")

    # if duplicate plot is set add it's output
    if config['modules']['duplicate_plot']:
        run_all_files.append("Benchmarks/metrics.duplicate_plot.txt")

    # if phasing is set add targets for hapcut and longhap analysis
    if config['modules']['phasing']:
        for CHR in CHROMS:
            run_all_files.append("Make_Vcf/step3_hapcut/step3_run_hapcut2_10xpipeline/s3_hapcut_output/{}_hapblock_{}".format(config['samples']['id'], CHR))

        run_all_files.extend(["Make_Vcf/step3_hapcut/step4_compare_with_refphasing/hapcut_eval.txt",
                              "Make_Vcf/step3_hapcut/step4_compare_with_refphasing/{}_hapcut.phased.vcf".format(config['samples']['id']),
                              "Make_Vcf/step3_hapcut/step3_run_hapcut2_10xpipeline/s3_hapcut_output/data_hapblock.bed",
                              ])

    if config['modules']['frag_de_novo']:
        run_all_files.append(["data/bc_seq_hash.json",
            "frag_de_novo/bc_list.txt",
            "frag_de_novo/done.fq"])

    if config['modules']['exon2fasta'] == True:
        run_all_files.extend(["Align/frag_coverage_done", "Align/frag_length_distribution_N100.pdf"])

    if config['modules']['consensus_fasta'] == True:
        run_all_files.extend([ 'Align/consensus/consensus.fasta', "Align/consensus/consensus_frag_length_distribution.pdf",'Align/consensus/consensus.fixRC.fasta',"Align/consensus/consensus.fixRC_SQANTI3_report.pdf"])

    RNA_16S_MODE = config['modules'].get('rna_16s', 'align_ref')
    if RNA_16S_MODE == 'meta_denovo':
        run_all_files.extend(["rna_16s/meta_denovo/contigs.fasta", "rna_16s/quast/meta_denovo/report.txt"])
    if RNA_16S_MODE == 'align_ref':
        run_all_files.extend([ "rna_16s/align_ref/abundance_align_ref.png"])
    if RNA_16S_MODE == 'frag_denovo':
        run_all_files.extend(["rna_16s/frag_denovo/all.contigs_max.fasta", "rna_16s/quast/frag_denovo/report.txt"])


    return run_all_files


# rule run all, the files above are the targets for snakemake
rule run_all:
    input:
        run_all_input
        

# Include other modules and rules to run

src_dir = config['params']['src_dir']
mrna_mapper = config['params']['mrna_mapper']
include: src_dir+"modules/clfr/calc_frag_len/calc_frag_len.smk"
include: src_dir+"modules/shared/metrics/metrics.smk"
include: src_dir+"modules/clfr/align/align_supp.smk"
if mrna_mapper =='minimap2':
    include: src_dir+"modules/clfr/align/align.minimap.smk"
else:
    include: src_dir+"modules/clfr/align/align.main.smk"
include: src_dir+"modules/shared/splitreads/splitreads.smk"
include: src_dir+"modules/shared/variant_calling/make_vcf.smk"
include: src_dir+"modules/clfr/consensus_fasta/consensus_fasta.smk"
include: src_dir+"modules/clfr/exon2fasta/exon2fasta.smk"
include: src_dir+"modules/clfr/rna_16s/rna_16s.smk"
include: src_dir+"modules/clfr/denovo/denovo_clfr.smk"
