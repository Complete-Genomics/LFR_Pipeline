
# Assembler-agnostic preprocessing: barcode selection -> read filtering ->
# sgrep TSV reformat. Shared by both denovo_clfr.smk (megahit) and
# denovo_olc.smk (seedext/OLC), which branch on config['frag_de_novo']['assembler']
# only from run_denovo_parallel onward -- included unconditionally so neither
# assembler-specific module needs its own copy of this chain.

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
        read="data/split_read_1_trimmed.fastq.gz"
    output:
        "denovo/data_R1_filtered.fastq.gz"
    params:
        bgzip = config['params']['bgzip'],
    run:
        if config['params']['sequence_type'].lower()=='pe':
            params.nth_line = count_lines()
            shell("awk -F '\t' -v OFS='\t' 'FNR == NR{{a[$1]++; next}} {{if (NR % 4 == {params.nth_line} ) {{ok=0; if($2 in a) ok = 1}}; if (ok == 1) print $0}}' {input.barcode_freq} <(zcat {input.read}) | {params.bgzip} -c > {output} ")
        elif config['params']['sequence_type'].lower()=='se':
            shell("touch denovo/data_R1_filtered.fastq.gz")


rule filter_reads2:
    input:
        barcode_freq="denovo/filtered_barcode_freq.txt",
        read="data/split_read_2_trimmed.fastq.gz"
    output:
        "denovo/data_R2_filtered.fastq.gz"
    params:
        bgzip = config['params']['bgzip'],
    run:
        params.nth_line = count_lines()
        shell("awk -F '\t' -v OFS='\t' 'FNR == NR{{a[$1]++; next}} {{if (NR % 4 == {params.nth_line} ) {{ok=0; if($2 in a) ok = 1}}; if (ok == 1) print $0}}' {input.barcode_freq} <(zcat {input.read}) | {params.bgzip} -c > {output} ")


rule reformat_fasta1:
    input:
        "denovo/data_R1_filtered.fastq.gz"
    output:
        "denovo/data_R1_sgrep.tsv"
    params:
        sequence_type = config['params']['sequence_type'].lower(),
        sort_mem = config['frag_de_novo'].get('sort_mem', '4G'),
        sort_tmp_dir = config['frag_de_novo'].get('sort_tmp_dir', 'denovo/tmp_sort')
    shell:
        """
        mkdir -p {params.sort_tmp_dir}
        if [[ "{params.sequence_type}" == "pe" ]]; then
            zcat {input} | \
            awk '{{if (NR%4==1) {{temp=$1; $1=$2; $2=temp}} }}1' | \
            awk '{{if (NR%4==1) line=line$0"\\t"; if (NR%4==2) {{print line$0; line=""}}}}' | \
            LC_ALL=C sort -T {params.sort_tmp_dir} -S {params.sort_mem} > {output}
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
    params:
        sort_mem = config['frag_de_novo'].get('sort_mem', '4G'),
        sort_tmp_dir = config['frag_de_novo'].get('sort_tmp_dir', 'denovo/tmp_sort')
    shell:
        """
        mkdir -p {params.sort_tmp_dir}
        zcat {input} | \
        awk '{{if (NR%4==1) {{temp=$1; $1=$2; $2=temp}} }}1' | \
        awk '{{if (NR%4==1) line=line$0"\\t"; if (NR%4==2) {{print line$0; line=""}}}}' | \
        LC_ALL=C sort -T {params.sort_tmp_dir} -S {params.sort_mem} > {output}
        """