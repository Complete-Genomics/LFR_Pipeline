# clfr: composition of the consensus transcripts against GENCODE (protein_coding / lncRNA / rRNA / globin / ...).
# Independent module: enabled by modules.genecode_consensus (workflows/clfr.smk adds the composition TSV/PNG to run_all).
# The read-level counterpart is genecode.smk (subsampled SE600 reads vs transcripts + rRNA + genome decoy).
#
# Input is consensus_fasta.smk's consensus/consensus.fixRC.fasta (built on demand if missing).
# Reference = GENCODE transcripts + rRNA, no genome decoy: consensus come from genome-mapped reads, so the question is
# only which annotated transcript each one matches. Built once, outside this workflow, by
#   build_genecode_ref.sh -t <genecode_transcripts> -r <genecode_rrna_fa> -o <...>.transcripts_rRNA.fa.gz
# and shared by all samples (params.genecode_consensus_ref, default: next to genecode_transcripts).
#   map_consensus_genecode      -> genecode/{id}.consensus.paf   (minimap2 -x asm20, all hits)
#   classify_consensus_genecode -> genecode/{id}.consensus.{composition,per_read,genes}.tsv + .png
# Classification (genecode_consensus_anno.py): best hit per consensus = most matching bases; globin by gene name and
# rRNA are kept whatever the coverage; otherwise best-hit query coverage < genecode_consensus_min_cov -> low_cov;
# consensus absent from the PAF -> unmapped (intron / intergenic / unannotated isoform / chimera all land here).
import os

GENECODE_CONS_ID = config['samples']['id']
_gcc_tx = config['params'].get('genecode_transcripts', '')
if not _gcc_tx and not config['params'].get('genecode_consensus_ref'):
    raise ValueError("modules.genecode_consensus needs params.genecode_transcripts in this run's config.yaml, "
                     "or a prebuilt params.genecode_consensus_ref; see config/config.yaml for the genecode_* block")
GENECODE_CONS_REF = config['params'].get('genecode_consensus_ref') or (
    _gcc_tx.replace('.transcripts.fa.gz', '.transcripts_rRNA.fa.gz') if _gcc_tx.endswith('.transcripts.fa.gz')
    else _gcc_tx + '.rRNA.fa.gz')
if not os.path.exists(GENECODE_CONS_REF):
    raise ValueError("GENCODE transcripts + rRNA reference not found: {}\nBuild it once with "
                     "modules/clfr/align/build_genecode_ref.sh -t <genecode_transcripts> -r <genecode_rrna_fa> -o {} "
                     "(no -g), or point params.genecode_consensus_ref at an existing file".format(
                         GENECODE_CONS_REF, GENECODE_CONS_REF))


rule map_consensus_genecode:
    input:
        ref = GENECODE_CONS_REF,
        fasta = "consensus/consensus.fixRC.fasta"
    output:
        paf = "genecode/{}.consensus.paf".format(GENECODE_CONS_ID)
    threads:
        config['threads'].get('consensus_minimap', 20)
    params:
        MINIMAP = config['params']['minimap2']
    log:
        "genecode/minimap_consensus_genecode.log"
    benchmark:
        "Benchmarks/main.map_consensus_genecode.txt"
    shell:
        """
        mkdir -p genecode
        # same preset as map_fasta_consensus; secondary hits stay in the PAF, the classifier keeps the best per query
        {params.MINIMAP} -x asm20 -t {threads} {input.ref} {input.fasta} > {output.paf} 2>> {log}
        """


rule classify_consensus_genecode:
    input:
        paf = "genecode/{}.consensus.paf".format(GENECODE_CONS_ID),
        fasta = "consensus/consensus.fixRC.fasta"
    output:
        tsv = "genecode/{}.consensus.composition.tsv".format(GENECODE_CONS_ID),
        per_read = "genecode/{}.consensus.per_read.tsv".format(GENECODE_CONS_ID),
        genes = "genecode/{}.consensus.genes.tsv".format(GENECODE_CONS_ID),
        png = "genecode/{}.consensus.png".format(GENECODE_CONS_ID)
    params:
        python = config['params']['general_python'],
        script = config['params']['src_dir'] + "modules/clfr/align/genecode_consensus_anno.py",
        prefix = "genecode/{}.consensus".format(GENECODE_CONS_ID),
        min_cov = config['params'].get('genecode_consensus_min_cov', 0.8)
    benchmark:
        "Benchmarks/main.classify_consensus_genecode.txt"
    shell:
        """
        {params.python} {params.script} {input.paf} {params.prefix} --fasta {input.fasta} --min_cov {params.min_cov}
        """
