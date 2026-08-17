
# Assembler-agnostic preprocessing: barcode selection -> read filtering ->
# sorted TSV reformat. Shared by both denovo_clfr.smk (megahit) and
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
    # Also rejects barcodes containing an ambiguous base call ('N') here, at
    # the one place both assemblers' barcode admission funnels through. An
    # 'N' surviving into the final assigned barcode should be unreachable
    # for one that passed whitelist correction, and signals a base-calling
    # artifact, not a real UMI (denovo.md sec 29, sec 60).
    #
    # An earlier version of this check instead rejected any barcode >=80%
    # one base (matching denovo_qc_probe.py's is_low_complexity(), fine for
    # that function's own job of nudging a diagnostic probe's sample away
    # from artifacts, since over-excluding there costs nothing). Reused as a
    # hard rejection gate here it was a false-positive generator: real
    # barcodes in this UMI design are legitimately single-base-dominated
    # (a plain "AAAAAAAAAAAAAAA" with normal reads and a normal ~1.3kb
    # contig is common), and checking that threshold against a real
    # 3000-barcode sample found 339 (11.3%) that would have been silently
    # dropped, all producing ordinary contigs -- see denovo.md sec 60 for
    # the full false-positive investigation, including why no fraction
    # threshold can cleanly separate real incidents from ordinary
    # low-diversity barcodes here. The two known incident barcodes with no
    # 'N' are instead bounded by denovo_seed_olc.py's _extend_one_contig
    # max_contig_len backstop, which doesn't need to guess the cause.
    input:
        "split_stat_read1.log"
    output:
        "denovo/filtered_barcode_freq.txt"
    benchmark:
        "Benchmarks/denovo.select_denovo_barcodes.txt"
    params:
        reads_per_BC = config['frag_de_novo']['reads_per_BC']
    shell:
        """
        mkdir -p denovo
        awk -F '\\t' -v cutoff={params.reads_per_BC} \
            'NR > 4 && NF >= 3 && $2 + 0 >= cutoff && $3 !~ /N/ {{print "BX:Z:" $3}}' \
            {input} > {output}
        """

rule filter_reads1:
    input:
        barcode_freq="denovo/filtered_barcode_freq.txt",
        read="data/split_read_1_trimmed.fastq.gz"
    output:
        "denovo/data_R1_filtered.fastq.gz"
    benchmark:
        "Benchmarks/denovo.filter_reads1.txt"
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
    benchmark:
        "Benchmarks/denovo.filter_reads2.txt"
    params:
        bgzip = config['params']['bgzip'],
    run:
        params.nth_line = count_lines()
        shell("awk -F '\t' -v OFS='\t' 'FNR == NR{{a[$1]++; next}} {{if (NR % 4 == {params.nth_line} ) {{ok=0; if($2 in a) ok = 1}}; if (ok == 1) print $0}}' {input.barcode_freq} <(zcat {input.read}) | {params.bgzip} -c > {output} ")


rule reformat_fasta1:
    input:
        "denovo/data_R1_filtered.fastq.gz"
    output:
        "denovo/data_R1_sorted.tsv"
    benchmark:
        "Benchmarks/denovo.reformat_fasta1.txt"
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
            touch denovo/data_R1_sorted.tsv
        else
            echo "Unknown type {params.sequence_type}" >&2;
            exit 1;
        fi
        """

rule reformat_fasta2:
    input:
        "denovo/data_R2_filtered.fastq.gz"
    output:
        "denovo/data_R2_sorted.tsv"
    benchmark:
        "Benchmarks/denovo.reformat_fasta2.txt"
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


## Zymo-relative sample gate for the hs1-like "depth inflated by short reads"
## failure mode (denovo.md sec 78).  Both depth and read-length distributions
## must have large adverse drift, and projected retained depth must be safe,
## before `auto` applies max_reads_per_umi=300 + read length >=300.  Automatic
## filtering is SE-only; PE remains report-only until pair-aware salvage exists.
checkpoint noisyPreprocessQC:
    input:
        reads="denovo/data_R2_sorted.tsv",
        depth_baseline=(config['params']['src_dir'] +
                        "/modules/clfr/denovo/models/" +
                        "zymo_reads_per_umi_distribution.tsv"),
        length_baseline=(config['params']['src_dir'] +
                         "/modules/clfr/denovo/models/" +
                         "zymo_read_length_distribution.tsv")
    output:
        "denovo/noisy_preprocess_decision.tsv"
    benchmark:
        "Benchmarks/denovo.noisyPreprocessQC.txt"
    params:
        python = config['params']['general_python'],
        src_dir = config['params']['src_dir'],
        mode = config['frag_de_novo'].get('noisy_preprocess', 'auto'),
        sequence_type = config['params']['sequence_type'].lower(),
        min_read_length = config['frag_de_novo'].get(
            'noisy_min_read_length', 300),
        max_reads_per_umi = config['frag_de_novo'].get(
            'noisy_max_reads_per_umi', 300)
    run:
        mode = str(params.mode).strip().lower()
        if mode not in ('off', 'report', 'auto'):
            raise ValueError(
                "frag_de_novo.noisy_preprocess={!r}; expected 'off', "
                "'report', or 'auto'".format(params.mode))
        if mode == 'auto' and params.sequence_type != 'se':
            print("noisy preprocess: PE auto salvage is not pair-aware; "
                  "running report-only and preserving reads")
            mode = 'report'
        params.effective_mode = mode
        shell(
            "{params.python} "
            "{params.src_dir}/modules/clfr/denovo/denovo_noisy_preprocess_qc.py "
            "--r2 {input.reads} --r2-format tsv "
            "--depth-baseline {input.depth_baseline} "
            "--read-length-baseline {input.length_baseline} "
            "--mode {params.effective_mode} "
            "--min-read-length {params.min_read_length} "
            "--max-reads-per-umi {params.max_reads_per_umi} "
            "--out {output}"
        )


def noisy_preprocess_reads(wildcards):
    decision_path = checkpoints.noisyPreprocessQC.get().output[0]
    with open(decision_path) as report:
        decision = dict(line.rstrip("\n").split("\t", 1)
                        for line in report if "\t" in line)
    if decision.get('action') == 'salvage':
        return "denovo/data_R2_salvaged.tsv"
    return "denovo/data_R2_sorted.tsv"


## This rule enters the DAG only when the checkpoint selected salvage. Healthy
## samples return data_R2_sorted.tsv directly, with no copy and no extra I/O.
rule applyNoisyPreprocess:
    input:
        reads="denovo/data_R2_sorted.tsv",
        decision="denovo/noisy_preprocess_decision.tsv"
    output:
        "denovo/data_R2_salvaged.tsv"
    params:
        python = config['params']['general_python'],
        src_dir = config['params']['src_dir']
    shell:
        "{params.python} "
        "{params.src_dir}/modules/clfr/denovo/denovo_noisy_preprocess_qc.py "
        "--r2 {input.reads} --r2-format tsv "
        "--apply-decision {input.decision} --reads-out {output}"
