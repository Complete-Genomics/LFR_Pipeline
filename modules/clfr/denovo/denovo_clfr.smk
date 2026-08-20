
SEQUENCE_TYPE = config['params']['sequence_type'].lower()
MRNA_MAPPER = config['params']['mrna_mapper'].lower()

megahit = config['params']['megahit']
# rg = config['params'].get('rg', 'rg')
bbduk = config['params']['bbduk']
bgzip = config['params']['bgzip']

# Preprocessing (barcode selection -> read filtering -> sorted TSV reformat) now
# lives in denovo_preprocess.smk, included unconditionally by workflows/clfr.smk
# so both the megahit and seedext/OLC assembler branches share it.

def count_fq_len():
    file_path = f'denovo/data_R2_sorted.tsv'
    with open(file_path, 'r') as f:
        total_lines = sum(1 for line in f)
    return total_lines


rule run_denovo_parallel:
    input:
        r1="denovo/data_R1_sorted.tsv",
        r2=noisy_preprocess_reads,
        noisy_qc="denovo/noisy_preprocess_decision.tsv"
    output:
        "denovo/frag_denovo_done"
    params:
        num_processes = config['frag_de_novo']['num_processes'],
        sequence_type = config['params']['sequence_type'],
        python = config['params']['general_python'],
        min_ctg_len = config['frag_de_novo']['min_ctg_len'],
        k_min = 41,
        k_max = 41,
        megahit = config['params']['megahit'],
        tmp_dir = config['frag_de_novo'].get('tmp_dir', '/dev/shm'),
        run_parallel = config['frag_de_novo'].get('run_parallel', False),
        src_dir = config['params']['src_dir']
    run:
        if not params.run_parallel:
            shell("mkdir -p denovo && touch {output}")
            return

        # chunk to run on mutiple nodes, or a dummy number 300000000 to run on single node
        params.end_idx = config['frag_de_novo'].get('end_idx')
        params.start_idx = config['frag_de_novo'].get('start_idx', 0)

        command = ["{params.python}",
                   "{params.src_dir}/modules/clfr/denovo/denovo_clfr_ram.py",
                   "--num_processes {params.num_processes} ",
                   "--sequence_type {params.sequence_type} ",
                   "--r1 {input.r1} ",
                   "--r2 {input.r2} ",
                   "--n_line_chunk 2000000 ",
                   "--start_idx {params.start_idx} ",
                   "--module denovo_parallel ",
                   "--min_ctg_len {params.min_ctg_len} ",
                   "--megahit {params.megahit} ",
                   "--tmp_dir {params.tmp_dir} ",
                   "--nth_of_nodes 0"]

        if params.end_idx not in (None, "", "all"):
            command.append("--end_idx {params.end_idx} ")
        shell(" ".join(command))

rule map_denovo:
    input:
        "denovo/frag_denovo_done"
    output:
        "denovo/denovo.paf"
    params:
        minimap = config['params']['minimap2'],
        refgenome = config['params']['ref_fa_mrna'],
        run_parallel = config['frag_de_novo'].get('run_parallel', False)
    run:
        if not params.run_parallel:
            shell("touch {output}")
            return

        if config['frag_de_novo']['denovo_type'] == 'correct_rc':
            command = ["{params.minimap} -x asm20 -t 20 ",
                   "{params.refgenome} ",
                   "denovo/final_contigs_0.fa > {output} "
                   ] 
            shell(" ".join(command))
        else:
            shell("touch denovo/denovo.paf")




rule correc_direction_denovo:
    input:
        denovo_fa = "denovo/frag_denovo_done",
        denovo_paf = "denovo/denovo.paf"
    output:
        "denovo/denovo.fixRC.fasta"
    params:
        python = config['params']['general_python'],
        flank_end = config['frag_de_novo']['flank_end'],
        adapter_seq = config['frag_de_novo']['adapter_seq'],
        src_dir = config['params']['src_dir'],
        run_parallel = config['frag_de_novo'].get('run_parallel', False)
    run:
        if not params.run_parallel:
            shell("touch {output}")
            return

        if config['frag_de_novo']['denovo_type'] == 'correct_rc':
            command = ["{params.python} ",
                    "{params.src_dir}/modules/clfr/denovo/denovo_supp.py ",
                        "--adapter_seq {params.adapter_seq} ",
                    "--flank_end {params.flank_end} ",
                    "--module fix_fa_rc "
                    ] 
            shell(" ".join(command))
        else:
            shell("cd denovo && ln -s final_contigs_0.fa denovo.fixRC.fasta ")

## output 1bc1frag.fasta, 
rule plot_denovo_frag_len_distribution:
    input:
        "denovo/frag_denovo_done",
        "denovo/denovo.fixRC.fasta"
    output:
        "denovo/frag_length_distribution.pdf", 
        "denovo/denovo.fixRC.1bc1frag.fasta"
    params:
        python = config['params']['general_python'],
        src_dir = config['params']['src_dir'],
        run_parallel = config['frag_de_novo'].get('run_parallel', False)
    run:
        if not params.run_parallel:
            shell("touch {output[0]} {output[1]}")
            return

        command = ["{params.python}",
                    "{params.src_dir}/modules/clfr/denovo/denovo_supp.py",
                    "--fasta denovo/final_contigs_0.fa",
                    "--outdir denovo/",
                    "--module metrics_basic"] 
        shell(" ".join(command))

rule filter_fasta_1k:
    input:
        inputfile="denovo/denovo.fixRC.1bc1frag.fasta"
    output:
        outputfile="denovo/denovo.fixRC.1bc1frag.1k.fasta"
    params:
        run_parallel = config['frag_de_novo'].get('run_parallel', False)
    run:
        if not params.run_parallel:
            shell("touch {output.outputfile}")
            return

        from Bio import SeqIO
        outfile=open(output.outputfile, 'w')

        with open(input.inputfile, "r") as input_fasta:
            for record in SeqIO.parse(input_fasta, "fasta"):
                if len(record.seq) >= 1000:
                    SeqIO.write(record, outfile, "fasta")
        outfile.close()

rule denovo_done:
    input:
        "denovo/denovo.fixRC.1bc1frag.1k.fasta"
    output:
        touch("denovo/done.fq")
