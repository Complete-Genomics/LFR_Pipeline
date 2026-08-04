
SEQUENCE_TYPE = config['params']['sequence_type'].lower()
MRNA_MAPPER = config['params']['mrna_mapper'].lower()

# Preprocessing (barcode selection -> read filtering -> sgrep TSV reformat) lives
# in denovo_preprocess.smk, included unconditionally by workflows/clfr.smk so both
# the megahit (denovo_clfr.smk) and seedext/OLC (this file) assembler branches
# share it -- this file only covers denovo_seed_olc.py onward.

rule run_denovoOLC_parallel:
    input:
        "denovo/data_R1_sgrep.tsv",
        "denovo/data_R2_sgrep.tsv"
    output:
        "denovo/frag_denovo_done"
    params:
        num_processes = config['frag_de_novo']['num_processes'],
        sequence_type = config['params']['sequence_type'],
        python = config['params']['general_python'],
        min_ctg_len = config['frag_de_novo']['min_ctg_len'],
        n_umi = config['frag_de_novo'].get('assembly_N_umi'),
        run_parallel = config['frag_de_novo'].get('run_parallel', False),
        max_contigs = config['frag_de_novo'].get('max_contigs', 8),
        src_dir = config['params']['src_dir']
    run:
        if not params.run_parallel:
            shell("mkdir -p denovo && touch {output}")
            return

        # chunk to run on multiple nodes, or a dummy number 300000000 to run on single node
        params.end_idx = config['frag_de_novo'].get('end_idx')
        params.start_idx = config['frag_de_novo'].get('start_idx', 0)

        # pure-Python OLC assembler, no megahit / no subprocess fork per UMI
        command = ["{params.python}",
                    "{params.src_dir}/modules/clfr/denovo/denovo_seed_olc.py",
                    "--num_processes {params.num_processes} ",
                    "--sequence_type {params.sequence_type} ",
                    "--n_line_chunk 2000000 ",
                    "--start_idx {params.start_idx} ",
                    "--min_ctg_len {params.min_ctg_len} ",
                    "--max_contigs {params.max_contigs} ",
                    "--nth_of_nodes 0"]
        if params.n_umi not in (None, "", "all"):
            command.append("--n {params.n_umi} ")

        if params.end_idx not in (None, "", "all"):
            command.append("--end_idx {params.end_idx} ")
        shell(" ".join(command))


## denovo_seed_olc.py already writes each barcode's contigs longest-first
## (k41_0 == longest), so "keep only k41_0" == "keep only the longest per UMI".
rule filterOLC_longest:
    input:
        "denovo/frag_denovo_done"
    output:
        "denovo/denovo.longest.fasta"
    params:
        run_parallel = config['frag_de_novo'].get('run_parallel', False)
    run:
        if not params.run_parallel:
            shell("touch {output}")
            return

        shell("""
            awk '
                /^>/ {{
                    if (seq != "" && keep) print header"\\n"seq
                    header=$0
                    keep = (header ~ /k41_0$/)
                    seq=""
                    next
                }}
                {{ seq = seq $0 }}
                END {{ if (seq != "" && keep) print header"\\n"seq }}
            ' denovo/final_contigs_0.fa > {output}
        """)

rule plotOLC_frag_len_distribution:
    input:
        "denovo/denovo.longest.fasta"
    output:
        "denovo/frag_length_distribution.pdf"
    params:
        python = config['params']['general_python'],
        src_dir = config['params']['src_dir'],
        run_parallel = config['frag_de_novo'].get('run_parallel', False)
    run:
        if not params.run_parallel:
            shell("touch {output}")
            return

        command = ["{params.python}",
                    "{params.src_dir}/modules/clfr/denovo/denovo_supp.py",
                    "--fasta denovo/denovo.longest.fasta",
                    "--outdir denovo/",
                    "--module metrics_olc"]
        shell(" ".join(command))

# rule filterOLC_fasta_1k:
#     input:
#         inputfile="denovo/denovo.longest.fasta"
#     output:
#         outputfile="denovo/denovo.longest.1k.fasta"
#     params:
#         run_parallel = config['frag_de_novo'].get('run_parallel', False)
#     run:
#         if not params.run_parallel:
#             shell("touch {output.outputfile}")
#             return

#         from Bio import SeqIO
#         outfile = open(output.outputfile, 'w')

#         with open(input.inputfile, "r") as input_fasta:
#             for record in SeqIO.parse(input_fasta, "fasta"):
#                 if len(record.seq) >= 1000:
#                     SeqIO.write(record, outfile, "fasta")
#         outfile.close()

rule olc_done:
    input:
        "denovo/denovo.longest.fasta",
        "denovo/frag_length_distribution.pdf"
    output:
        touch("denovo/done.fq")
