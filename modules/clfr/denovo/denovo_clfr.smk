
SEQUENCE_TYPE = config['params']['sequence_type'].lower()
MRNA_MAPPER = config['params']['mrna_mapper'].lower()

megahit = config['params']['megahit']
# rg = config['params'].get('rg', 'rg')
bbduk = config['params']['bbduk']
bgzip = config['params']['bgzip']

def count_lines():
    ## count filtered_barcode_freq.txt to determine starts of fq header line in filter_reads1 function
    file_path = "denovo/filtered_barcode_freq.txt"
    with open(file_path, 'r') as f:
        line_count = sum(1 for _ in f)
    return (((line_count%4)+1)%4)

rule select_denovo_barcodes:
    input:
        "split_stat_read1.log"
    output:
        "denovo/filtered_barcode_freq.txt"
    params:
        reads_per_BC = config['frag_de_novo']['reads_per_BC']
    shell:
        """
        mkdir -p denovo
        awk -F '\\t' -v cutoff={params.reads_per_BC} \
            'NR > 4 && NF >= 3 && $2 + 0 >= cutoff {{print "BX:Z:" $3}}' \
            {input} > {output}
        """

rule filter_reads1:
    input:
        barcode_freq="denovo/filtered_barcode_freq.txt",
        read="data/split_read.1.fq.gz"
    output:
        "denovo/data_R1_filtered.fastq.gz"
    params:
        bgzip = bgzip,
    run:
        if config['params']['sequence_type'].lower()=='pe':
            params.nth_line = count_lines()
            shell("awk -F '\t' -v OFS='\t' 'FNR == NR{{a[$1]++; next}} {{if (NR % 4 == {params.nth_line} ) {{ok=0; if($2 in a) ok = 1}}; if (ok == 1) print $0}}' {input.barcode_freq} <(zcat {input.read}) | {params.bgzip} -c > {output} ")
        elif config['params']['sequence_type'].lower()=='se':
            shell("touch denovo/data_R1_filtered.fastq.gz")


rule filter_reads2:
    input:
        barcode_freq="denovo/filtered_barcode_freq.txt",
        read="data/split_read.2.fq.gz"
    output:
        "denovo/data_R2_filtered.fastq.gz"
    params:
        bgzip = bgzip,
    run:
        params.nth_line = count_lines()
        shell("awk -F '\t' -v OFS='\t' 'FNR == NR{{a[$1]++; next}} {{if (NR % 4 == {params.nth_line} ) {{ok=0; if($2 in a) ok = 1}}; if (ok == 1) print $0}}' {input.barcode_freq} <(zcat {input.read}) | {params.bgzip} -c > {output} ")



# rule trim_reads:
#     input:
#         f1="denovo/data_R1_filtered.fastq.gz",
#         f2="denovo/data_R2_filtered.fastq.gz"
#     output:
#         t1= "denovo/data_R1_filtered_trimmed.fastq.gz",
#         t2= "denovo/data_R2_filtered_trimmed.fastq.gz"
#     params:
#         bbduk = bbduk,
#         sequence_type= config['params']['sequence_type'].lower()
#     shell:
#         """
#         if [[ "{params.sequence_type}" == "pe" ]]; then
#             {params.bbduk} in1={input.f1} in2={input.f2} out1={output.t1} out2={output.t2} qtrim=rl 
#         elif [[ "{params.sequence_type}" == "se" ]]; then
#             {params.bbduk} in={input.f2} out={output.t2} qtrim=rl && touch denovo/data_R1_filtered_trimmed.fastq.gz
#         else
#             echo "Unknown type {params.sequence_type}" >&2;
#             exit 1;
#         fi
#         """


rule reformat_fasta1:
    input:
        "denovo/data_R1_filtered.fastq.gz"
    output:
        "denovo/data_R1_sgrep.tsv"
    params:
        sequence_type = config['params']['sequence_type'].lower()
    shell:
        """
        if [[ "{params.sequence_type}" == "pe" ]]; then
            zcat {input} | \
            awk '{{if (NR%4==1) {{temp=$1; $1=$2; $2=temp}} }}1' | \
            awk '{{if (NR%4==1) line=line$0"\\t"; if (NR%4==2) {{print line$0; line=""}}}}' | \
            sort -T /tmp/ -S 60% > {output}
        elif [[ "{params.sequence_type}" == "se" ]]; then
            touch denovo/data_R1_sgrep.tsv
        else
            echo "Unknown type {params.sequence_type}" >&2;
            exit 1;
        fi
        """

rule reformat_fasta2:
    input:
        "denovo/data_R2_filtered.fastq.gz"
    output:
        "denovo/data_R2_sgrep.tsv"
    shell:
        """
        zcat {input} | \
        awk '{{if (NR%4==1) {{temp=$1; $1=$2; $2=temp}} }}1' | \
        awk '{{if (NR%4==1) line=line$0"\\t"; if (NR%4==2) {{print line$0; line=""}}}}' | \
        sort -T /tmp/ -S 60% > {output}
        """

def count_fq_len():
    file_path = f'denovo/data_R2_sgrep.tsv'
    with open(file_path, 'r') as f:
        total_lines = sum(1 for line in f)
    return total_lines


rule run_denovo_parallel:
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
        k_min = 41,
        k_max = 41,
        megahit = config['params']['megahit'],
        tmp_dir = config['frag_de_novo'].get('tmp_dir', '/dev/shm'),
        src_dir = config['params']['src_dir']
    run:
        # chunk to run on mutiple nodes, or a dummy number 300000000 to run on single node
        params.end_idx = config['frag_de_novo'].get('end_idx')
        params.start_idx = config['frag_de_novo'].get('start_idx', 0)
        command = ["{params.python}",
                   "{params.src_dir}/modules/clfr/denovo/denovo_clfr_ram.py",
                   "--num_processes {params.num_processes} ",
                   "--sequence_type {params.sequence_type} ",
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
        refgenome = config['params']['ref_fa_mrna']
    run:
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
        src_dir = config['params']['src_dir']
    run:
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
    run:
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
    run:
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
