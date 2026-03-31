# include the config file
configfile: "config.yaml"

SEQUENCE_TYPE =config['params']['sequence_type'].lower()
MRNA_MAPPER =config['params']['mrna_mapper'].lower()
# Define the ref based on the config file
# Sort of acts like a global variable so you don't need to always type the whole thing
if config['params']['library_species']=='human':
    if config['params']['library_type']=='gdna':
        REF = config['params']['ref_fa']
    elif config['params']['library_type']=='mrna':
        REF = config['params']['ref_fa_mrna']
elif config['params']['library_species']=='other':
    REF = config['params']['ref_fa_other']

MINREADS = config['calc_frag']['minreads']
MIN_FRAG = config['calc_frag']['min_frag']
MAPPING_QUALITY = config['calc_frag']['mapping_quality']
BC_SPLIT = config['modules']['BCsplit']
if config['params']['bc_condition']=='pcrfree':
    BC_SPLIT = False

# SE needs half reads
if config['params']['bc_condition']=='standard' or config['params']['bc_condition']=='pcrfree' :
    if SEQUENCE_TYPE=='pe':
        MINREADS = 4
    else:
        MINREADS = 2

# This defines the logic for which chromosomes to use
# Check if the chroms are defined
if not config['samples']['chroms']:
    CHROMS = []
    # If they aren't defined try and scrape the fasta.fai index
    try:
        with open(REF+".fai", "r") as fasta_index:
            for line in fasta_index:
                CHROMS.append(line.strip().split()[0])

    except Exception as e:
        print(f"{REF}.fai cannot be opened for parsing.")
        print(f"Exception: {e}")
        sys.exit(1)
# if the chroms are defined just use those
else:
    CHROMS = config['samples']['chroms']


# define a function to return target files based on config settings
def run_all_input(wildcards):
    # split only
    run_all_files = ["data/split_read.2.fq.gz",]

    # with alignment
    if config['modules']['mapping']==True:
        run_all_files.extend(["keep/Align/{}.sort.bam".format(config['samples']['id'])])

    if config['modules']['eval_metrics']==True:
        if BC_SPLIT==True:
            run_all_files.extend(['Align/mapped_uniq_bc_bases_count.txt'])
        run_all_files.extend([
                        "summary_report.txt",
                        "Align/{}_dedup_metrics2.txt".format(config['samples']['id']),
                        "Align/flagstat_metric.txt",
                        "Align/gatk_metrics_{}.alignment_summary_metrics".format(config['samples']['id'])])
        
    # if stLFR is set to true add fragment calculation metrics
    if config['modules']['LFR']:
        calc_frag_file = ["Calc_Frag_Length/frag_length_distribution.pdf",
                          "Calc_Frag_Length/n_read_distribution.pdf",
                          "Calc_Frag_Length/frag_and_bc_summary.txt",
                          "Calc_Frag_Length/frag_summary_minreads"+str(MINREADS)+".txt"]
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

    # phasing use vc called
    if config['modules']['phasing']==True:
        if config['modules']['variant_calling']==True:
            for CHR in CHROMS:
                run_all_files.append("Make_Vcf/step3_hapcut/step3_run_hapcut2_10xpipeline/s3_hapcut_output/{}_hapblock_{}".format(config['samples']['id'], CHR))

            run_all_files.extend(["Make_Vcf/step3_hapcut/step4_compare_with_refphasing/hapcut_eval.txt",
                                "Make_Vcf/step3_hapcut/step4_compare_with_refphasing/{}_hapcut.phased.vcf".format(config['samples']['id']),
                                "Make_Vcf/step3_hapcut/step3_run_hapcut2_10xpipeline/s3_hapcut_output/data_hapblock.bed",
                                ])
        else:
            for CHR in CHROMS:
                run_all_files.append("Benchmarks/step3_hapcut/step3_run_hapcut2_10xpipeline/s3_hapcut_output/{}_hapblock_{}".format(config['samples']['id'], CHR))
            run_all_files.extend(["Benchmarks/step3_hapcut/step4_compare_with_refphasing/hapcut_eval.txt",
                                "Benchmarks/step3_hapcut/step4_compare_with_refphasing/hapcut_eval.txt",
                                ])


    ## if BCgDNA, no run of other modules
    if config['modules']['BCgDNA'] == True:
        _txt = ["Align/{}.sort.r2.bam.pos.txt".format(i) for i in ["nobc1", "nobc2bc3", "noall3", "nobc3"]]
        run_all_files=["data/split_read.2.fq.gz"]+_txt

    return run_all_files


# rule run all, the files above are the targets for snakemake
rule run_all:
    input:
        run_all_input
        

# Include other modules and rules to run
src_dir = config['params']['src_dir']
include: src_dir+"modules/stlfr/calc_frag_len/calc_frag_len.smk"
include: src_dir+"modules/shared/variant_calling/make_vcf.smk"
include: src_dir+"modules/shared/metrics/metrics.smk"
include: src_dir+"modules/stlfr/align/align.smk"
include: src_dir+"modules/shared/splitreads/splitreads.smk"
include: src_dir+"modules/shared/phasing/phasing.smk"
include: src_dir+"modules/shared/phasing/phasing_benchmark.smk"
include: src_dir+"modules/stlfr/bcgdna/BCgDNA.smk"
