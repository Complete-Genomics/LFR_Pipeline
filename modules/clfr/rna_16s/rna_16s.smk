# rna_16s.smk - 16S rRNA Analysis Pipeline, as a module of a larger workflow
# Input: Align/{SAMPLE_ID}.sort.bam from upstream step
# Implements 3 methods: meta_denovo, align_ref, frag_denovo
# Selected by config['modules']['rna_16s']: 'meta_denovo' | 'align_ref' | 'frag_denovo' | False
import os
from pathlib import Path

# ============================================================================
# Configuration from config.yaml
# ============================================================================
RNA_16S_MODE = config['modules'].get('rna_16s', 'align_ref')
SAMPLE_ID = config['samples'].get('id', 'data')
REF_16S = config['params'].get('ref_fa_other', '')
SEQUENCE_TYPE = config['params'].get('sequence_type', 'pe').lower()

# Frag denovo params
N_FRAG = config['params'].get('n_frag', 100)
MIN_READS_BC = config['params'].get('min_reads_bc', 50)
MAX_READS_BC = config['params'].get('max_reads_bc', 1000)

# Tool paths from config
PYTHON = config['params'].get('general_python', 'python3')
SPADES = config['frag_de_novo'].get('denovo_assembler', 'SPAdes-3.14.0-Linux/bin/spades.py')
QUAST = config['frag_de_novo'].get('quast_dir', 'quast/quast.py')
MINIMAP = config['frag_de_novo'].get('minimap', 'minimap2-2.16_x64-linux/minimap2')
BWA='bwa'
SAMTOOLS='samtools'
SCRIPT_DIR = str(Path(workflow.basedir).parent / 'rna_16s')

# Thread counts
THREADS_BWA = config['threads'].get('bwa', 20)
THREADS_QUAST = config['threads'].get('quast', 4)


# ============================================================================
# Method 1: Meta Denovo (MetaSPAdes assembly)
# ============================================================================
rule meta_denovo:
    input:
        r1="data/split_read.1.fq.gz",
        r2="data/split_read.2.fq.gz"
    output:
        contigs="rna_16s/meta_denovo/contigs.fasta"
    params:
        spades=SPADES,
        out_dir="rna_16s/meta_denovo",
        threads=THREADS_BWA
    shell:
        """
        {params.spades} -k 55,77,99,127 --meta \
            -t {params.threads} \
            -1 {input.r1} -2 {input.r2} \
            -o {params.out_dir}
        """


# ============================================================================
# Method 2: Align to Reference (BWA + SAMtools), get .bam from align.smk
# ============================================================================

rule align_idxstats:
    input:
        bam=f"Align/{SAMPLE_ID}.sort.bam",
        bai=f"Align/{SAMPLE_ID}.sort.bam.bai"
    output:
        stats=f"Align/{SAMPLE_ID}.idxstats.txt"
    params:
        samtools=SAMTOOLS
    shell:
        """
        {params.samtools} idxstats {input.bam} > {output.stats}
        """


# ============================================================================
# Method 3: Frag Denovo (per-barcode assembly)
# ============================================================================

# Step 3.2: Extract FASTQ for all selected barcodes
rule frag_denovo_bc2fq:
    input:
        r1="data/split_read.1.fq.gz",
        r2="data/split_read.2.fq.gz",
        bc_list="rna_16s/frag_denovo/bc4denovo.txt"
    output:
        done=touch("rna_16s/frag_denovo/fq/bc2fq.done")
    params:
        python=PYTHON,
        script=SCRIPT_DIR + '/bc2fq.py',
        outdir="rna_16s/frag_denovo/fq"
    shell:
        """
        mkdir -p {params.outdir}
        {params.python} {params.script} --r1r2 1 --dirname ./ --bc_list {input.bc_list} --outdir {params.outdir}
        {params.python} {params.script} --r1r2 2 --dirname ./ --bc_list {input.bc_list} --outdir {params.outdir}
        """


# Step 3.3: SPAdes assembly + get max contig for each barcode
rule frag_denovo_spades:
    input:
        done="rna_16s/frag_denovo/fq/bc2fq.done",
        r1="rna_16s/frag_denovo/fq/{bc}.1.fq.gz",
        r2="rna_16s/frag_denovo/fq/{bc}.2.fq.gz"
    output:
        contigs="rna_16s/frag_denovo/spades/{bc}/contigs_max.fasta"
    params:
        spades=SPADES,
        python=PYTHON,
        get_max_fa=SCRIPT_DIR + '/get_max_fa.py',
        out_dir="rna_16s/frag_denovo/spades/{bc}",
        threads=4
    shell:
        """
        {params.spades} -k 55 -t {params.threads} \
            -1 {input.r1} -2 {input.r2} \
            -o {params.out_dir}
        {params.python} {params.get_max_fa} \
            --input {params.out_dir}/contigs.fasta \
            --output {output.contigs}
        """


# Step 3.4: Merge all barcode contigs (uses checkpoint for dynamic barcode list)
def get_bc_contigs(wildcards):
    bc_list_file = checkpoints.frag_denovo_bc_stats.get().output.bc_list
    with open(bc_list_file) as f:
        bcs = [line.strip() for line in f if line.strip()]
    return expand("rna_16s/frag_denovo/spades/{bc}/contigs_max.fasta", bc=bcs)

# Step 3.1: Select barcodes with appropriate read count (checkpoint for dynamic DAG)
rule frag_denovo_bc_stats:
    input:
        "split_stat_read1.log"
    output:
        bc_list="rna_16s/frag_denovo/bc4denovo.txt"
    params:
        min_reads=MIN_READS_BC,
        max_reads=MAX_READS_BC,
        n_frag=N_FRAG
    run:
        bcs = []
        with open(input[0]) as f:
            for line in f:
                cols = line.strip().split('\t')
                if len(cols) >= 3:
                    try:
                        n_reads = int(cols[1])
                        if params.min_reads <= n_reads <= params.max_reads:
                            bcs.append(cols[2])
                    except ValueError:
                        continue
        with open(output.bc_list, 'w') as out:
            for bc in bcs[:params.n_frag]:
                out.write(bc + '\n')

rule frag_denovo_merge:
    input:
        get_bc_contigs
    output:
        "rna_16s/frag_denovo/all.contigs_max.fasta"
    shell:
        """
        cat {input} > {output}
        """


# ============================================================================
# Quast Evaluation
# ============================================================================
rule quast_meta_denovo:
    input:
        "rna_16s/meta_denovo/contigs.fasta"
    output:
        report="rna_16s/quast/meta_denovo/report.txt"
    params:
        ref=REF_16S,
        out_dir="rna_16s/quast/meta_denovo",
        quast=QUAST,
        threads=THREADS_QUAST
    shell:
        """
        {params.quast} {input} -o {params.out_dir} -r {params.ref} \
            --min-contig 100 --no-snps --threads {params.threads} --space-efficient
        """


rule quast_frag_denovo:
    input:
        "rna_16s/frag_denovo/all.contigs_max.fasta"
    output:
        report="rna_16s/quast/frag_denovo/report.txt"
    params:
        ref=REF_16S,
        out_dir="rna_16s/quast/frag_denovo",
        quast=QUAST,
        threads=THREADS_QUAST
    shell:
        """
        {params.quast} {input} -o {params.out_dir} -r {params.ref} \
            --min-contig 100 --no-snps --threads {params.threads} --space-efficient
        """


# ============================================================================
# ZymoBIOMICS Abundance Stats (for control samples)
# Depends on quast_frag_denovo which produces the coords file as a byproduct
# ============================================================================
rule abundance_frag_denovo:
    input:
        quast_report="rna_16s/quast/frag_denovo/report.txt"
    output:
        "rna_16s/abundance_frag_denovo.txt"
    params:
        python=PYTHON,
        script=SCRIPT_DIR + '/rna_16s.py',
        coords_dir="rna_16s/quast/frag_denovo"
    shell:
        """
        {params.python} {params.script} --module metrics_ref_zymo \
            --outdir {params.coords_dir} > {output}
        """

rule abundance_meta_denovo:
    input:
        quast_report="rna_16s/quast/meta_denovo/report.txt"
    output:
        "rna_16s/abundance_meta_denovo.txt"
    params:
        python=PYTHON,
        script=SCRIPT_DIR + '/rna_16s.py',
        coords_dir="rna_16s/quast/meta_denovo"
    shell:
        """
        {params.python} {params.script} --module metrics_ref_zymo \
            --outdir {params.coords_dir} > {output}
        """

rule abundance_align_ref:
    input:
        bam=f"Align/{SAMPLE_ID}.sort.bam"
    output:
        csv="rna_16s/align_ref/abundance_align_ref.csv",
        png="rna_16s/align_ref/abundance_align_ref.png"
    params:
        python=PYTHON,
        script=SCRIPT_DIR + '/align_ref.py',
        outdir="rna_16s/align_ref",
        ref_fasta=REF_16S
    shell:
        """
        {params.python} {params.script} --bam {input.bam} \
            --outdir {params.outdir} --ref_fasta {params.ref_fasta}
        """