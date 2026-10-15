'''
Standalone patch (not included by clfr.smk): drop globin consensus fragments and re-plot the length distribution.

  consensus/consensus.fixRC.fasta + consensus/consensus.paf
    -> consensus/consensus.fixRC.NOglobin.fasta
       consensus/globin_removed_ids.tsv, consensus/globin_removed_stats.tsv
       consensus/consensus_frag_length_distribution_NOglobin.png (+ .txt)

Run from the analysis directory (the one that holds consensus/ and config.yaml):
  snakemake -s <pipeline>/modules/clfr/consensus_fasta/globin.smk -j 4
Reads from config.yaml: params.general_python, params.minimap2_mrna_index or params.ref_fa_mrna
(the transcriptome consensus.paf was mapped to), calc_frag.minreads_fasta (default 100),
optional params.globin_genes. The helper scripts are taken from this file's own directory.
'''
from pathlib import Path

configfile: "config.yaml"

SCRIPT_DIR = Path(workflow.basedir)
PYTHON = config.get('params', {}).get('general_python', 'python')
GLOBIN_GENES = config.get('params', {}).get('globin_genes', "HBA1,HBA2,HBB,HBD,HBE1,HBG1,HBG2,HBM,HBQ1,HBZ")

rule all:
    input:
        "consensus/consensus.fixRC.NOglobin.fasta",
        "consensus/consensus_frag_length_distribution_NOglobin.png"

rule remove_globin_consensus:
    input:
        fasta = "consensus/consensus.fixRC.fasta",
        paf = "consensus/consensus.paf"
    output:
        fasta = "consensus/consensus.fixRC.NOglobin.fasta",
        ids = "consensus/globin_removed_ids.tsv",
        stats = "consensus/globin_removed_stats.tsv"
    params:
        script = str(SCRIPT_DIR / "globin_filter.py"),
        # same reference map_fasta_consensus mapped to
        ref = config['params'].get('minimap2_mrna_index') or config['params']['ref_fa_mrna'],
        genes = GLOBIN_GENES
    benchmark:
        "Benchmarks/consensus_fasta.remove_globin_consensus.txt"
    shell:
        "{PYTHON} {params.script} "
        "--fasta {input.fasta} --paf {input.paf} --ref {params.ref} --genes {params.genes} "
        "--out_fasta {output.fasta} --out_ids {output.ids} --out_stats {output.stats}"

# same plot as fasta_frag_len_distribution_consensus (exon2fasta.py metrics_basic), on the globin-free fasta.
# metrics_basic writes into --outdir and appends to a shared frag_length_distribution.txt, so run it in a
# scratch dir and move the png/txt out under the NOglobin names.
rule fasta_frag_len_distribution_NOglobin:
    input:
        "consensus/consensus.fixRC.NOglobin.fasta"
    output:
        png = "consensus/consensus_frag_length_distribution_NOglobin.png",
        txt = "consensus/consensus_frag_length_distribution_NOglobin.txt"
    params:
        script = str(SCRIPT_DIR / "exon2fasta.py"),
        minreads_fasta = config.get('calc_frag', {}).get('minreads_fasta', 100),
        tmpdir = "consensus/tmp_NOglobin_plot"
    benchmark:
        "Benchmarks/consensus_fasta.fasta_frag_len_distribution_NOglobin.txt"
    shell:
        """
        rm -rf {params.tmpdir} && mkdir -p {params.tmpdir}
        {PYTHON} {params.script} \
            --fasta {input} --outdir {params.tmpdir} --name consensus_NOglobin \
            --minreads_fasta {params.minreads_fasta} --module metrics_basic
        mv {params.tmpdir}/consensus_NOglobin_frag_length_distribution.png {output.png}
        mv {params.tmpdir}/frag_length_distribution.txt {output.txt}
        rm -rf {params.tmpdir}
        """
