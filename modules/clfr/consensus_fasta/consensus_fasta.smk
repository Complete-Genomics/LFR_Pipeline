'''
part of clfr workflow, get (transcriptome) fasta
'''

from pathlib import Path

NUM_SPLITS_CONSENSUS=5
wildcard_constraints:
    split_idx = r"\d+",
    i = r"\d+"

# ---- drop globin + globin pseudogene reads from the bam consensus starts from ----
# Mapping stays on the standard reference. This only saves consensus time (reformat / sort / per-UMI assembly);
# the unfiltered Align/{id}.sort.removedup_rm000.bam is untouched and still used by everything else.
#   params.consensus_exclude_globin: true|false   default false (off)
#   params.consensus_exclude_bed: optional BED (reference contig names); replaces the packaged one; only used when consensus_exclude_globin is true
# Default BED = globin_regions.bed next to this file (GRCh38, chr-prefixed; how it was built: globin_regions.README.md).
# It is a checked-in file, not rebuilt per run. For a reference with other contig names, give your own BED.
# get_chr_dict (samtools_idx.txt: per-chr read counts that set the 5 split windows) reads the SAME filtered bam,
# so the counts match what get_consensus_fasta sees.
CONSENSUS_EXCLUDE_BED_USER = config['params'].get('consensus_exclude_bed') or ''
CONSENSUS_EXCLUDE_GLOBIN = bool(config['params'].get('consensus_exclude_globin', False))
CONSENSUS_EXCLUDE_BED = CONSENSUS_EXCLUDE_BED_USER or str(Path(config['params']['src_dir']) / "modules" / "clfr" / "consensus_fasta" / "globin_regions.bed")
CONSENSUS_BAM_PAT = "Align/{id}.sort.removedup_rm000" + (".noglobin" if CONSENSUS_EXCLUDE_GLOBIN else "") + ".bam"
CONSENSUS_BAM = CONSENSUS_BAM_PAT.replace("{id}", config['samples']['id'])

# ---- optional salvage of low-depth UMI reads (params.consensus_salvage, default false) ----
# Moves reads of UMIs with < reads_per_BC reads into a valid UMI on the same chr when the UMI differs by <= salvage_max_mm
# substitutions and the read lies inside that UMI's fragment (see salvage_umi.py). get_consensus_fasta then reads the
# salvaged name-sorted bam instead of Align/tmp/{id}_{chr}.name.sort.bam.
CONSENSUS_SALVAGE = bool(config['params'].get('consensus_salvage', False))
CONSENSUS_NAMESORT_BAM = ("Align/tmp/{id}_{chr}.name.salvage.sort.bam" if CONSENSUS_SALVAGE
                          else "Align/tmp/{id}_{chr}.name.sort.bam")

if CONSENSUS_EXCLUDE_GLOBIN:
    rule filter_globin_bam:
        input:
            bam = "Align/{id}.sort.removedup_rm000.bam",
            bai = "Align/{id}.sort.removedup_rm000.bam.bai",
            bed = CONSENSUS_EXCLUDE_BED
        output:
            bam = "Align/{id}.sort.removedup_rm000.noglobin.bam",
            bai = "Align/{id}.sort.removedup_rm000.noglobin.bam.bai",
            stats = "Align/{id}.globin_filter_stats.txt"
        threads:
            config['threads'].get('consensus_name_sort', 2)
        params:
            samtools = config['params'].get('samtools', 'samtools')
        benchmark:
            "Benchmarks/consensus_fasta.filter_globin_bam.{id}.txt"
        shell:
            """
            set -euo pipefail
            # -L keeps reads overlapping the BED; -U receives everything else, which is the output we want (do not add -M: it
            # switches -L to index-region iteration and -U would then miss the reads outside the BED)
            {params.samtools} view -@ {threads} -b -L {input.bed} -U {output.bam} -o /dev/null {input.bam}
            {params.samtools} index -@ {threads} {output.bam}
            before=$({params.samtools} idxstats {input.bam} | awk '{{n+=$3}} END{{print n}}')
            after=$({params.samtools} idxstats {output.bam} | awk '{{n+=$3}} END{{print n}}')
            printf "mapped_before\\t%s\\nmapped_after\\t%s\\nremoved\\t%s\\n" $before $after $((before-after)) > {output.stats}
            """

rule get_chr_dict:
    input:
        bam=CONSENSUS_BAM,
        bai=CONSENSUS_BAM + ".bai"
    output:
        "Align/samtools_idx.txt"
    params:
        samtools = config['params'].get('samtools', 'samtools')
    benchmark:
        "Benchmarks/consensus_fasta.get_chr_dict.txt"
    shell:
        "{params.samtools} idxstats {input.bam} > {output}"

rule split_consensus_bam:
    input:
        bam = CONSENSUS_BAM_PAT,
        bai = CONSENSUS_BAM_PAT + ".bai"
    output:
        bam = "Make_Vcf/step3_hapcut/step1_modify_bam/{id}_sort.markdup_{chr}.bam"
    params:
        samtools = config['params'].get('samtools', 'samtools')
    benchmark:
        "Benchmarks/consensus_fasta.split_consensus_bam.{id}.{chr}.txt"
    shell:
        "{params.samtools} view -bh {input.bam} {wildcards.chr} > {output.bam}"
        
rule reformat_readid:
    input:
        bam = "Make_Vcf/step3_hapcut/step1_modify_bam/{id}_sort.markdup_{chr}.bam"
    output:
        "Align/tmp/{id}_{chr}.name.bam"
    params:
        python = config['params']['general_python'],
        src_dir = config['params']['src_dir'],
        seq_type = SEQUENCE_TYPE, 
        umi_len = config['params']['cbc_len'],
    benchmark:
        "Benchmarks/consensus_fasta.reformat_readid.{id}.{chr}.txt"
    run:
        command = ["{params.python}",
                    "{params.src_dir}/modules/clfr/consensus_fasta/consensus_fasta_supp.py",
                    "--module reformat_readid",
                    "--seq_type {params.seq_type} ",
                    "--umi_len {params.umi_len} ",
                    "--chr_name {wildcards.chr}",
                    "--input_bam {input.bam}",
                    "--output_bam {output}"] 
        shell(" ".join(command))

        # "python {params.src_dir}/modules/clfr/consensus_fasta/consensus_fasta_supp.py --module reformat_readid --chr_name {wildcards.chr}"

rule sort_reformated_bam:
    input:
        "Align/tmp/{id}_{chr}.name.bam"
    output:
        "Align/tmp/{id}_{chr}.name.sort.bam"
    params:
        samtools = config['params'].get('samtools', 'samtools')
    threads:
        config['threads'].get('consensus_name_sort', 2)
    benchmark:
        "Benchmarks/consensus_fasta.sort_reformated_bam.{id}.{chr}.txt"
    shell:
        "{params.samtools} sort -@ {threads} -n -o {output} {input} "

# samtools sort -n -o Align/${chr}.name.sort.bam Align/${chr}.name.bam

rule salvage_umi_map:
    input:
        bam = "Align/tmp/{id}_{chr}.name.sort.bam",
        split_log = "split_stat_read1.log"
    output:
        reads = "consensus/salvage/{id}_{chr}.salvage_reads.tsv",
        stats = "consensus/salvage/{id}_{chr}.salvage_stats.tsv"
    threads: 1
    resources:
        mem_mb = config['params'].get('salvage_mem_mb', 8000)
    params:
        python = config['params']['general_python'],
        src_dir = config['params']['src_dir'],
        samtools = config['params'].get('samtools', 'samtools'),
        min_reads = config['params']['reads_per_BC'],
        max_mm = config['params'].get('salvage_max_mm', 2),
        pad = config['params'].get('salvage_pad', 0)
    benchmark:
        "Benchmarks/consensus_fasta.salvage_umi_map.{id}.{chr}.txt"
    shell:
        "{params.python} {params.src_dir}/modules/clfr/consensus_fasta/salvage_umi.py --module map "
        "--input_bam {input.bam} --chrom {wildcards.chr} --split_log {input.split_log} "
        "--min_reads {params.min_reads} --max_mm {params.max_mm} --pad {params.pad} "
        "--samtools {params.samtools} --output_reads {output.reads} --output_stats {output.stats}"

# renaming changes the name order, so salvage_umi_apply re-sorts by name (consensus_fasta.py groups consecutive reads by UMI)
rule salvage_umi_apply:
    input:
        bam = "Align/tmp/{id}_{chr}.name.sort.bam",
        reads = "consensus/salvage/{id}_{chr}.salvage_reads.tsv"
    output:
        bam = "Align/tmp/{id}_{chr}.name.salvage.sort.bam"
    params:
        python = config['params']['general_python'],
        src_dir = config['params']['src_dir'],
        samtools = config['params'].get('samtools', 'samtools')
    threads:
        config['threads'].get('consensus_name_sort', 2)
    benchmark:
        "Benchmarks/consensus_fasta.salvage_umi_apply.{id}.{chr}.txt"
    shell:
        "{params.python} {params.src_dir}/modules/clfr/consensus_fasta/salvage_umi.py --module apply "
        "--input_bam {input.bam} --reads_table {input.reads} --output_bam {output.bam}.unsorted.bam && "
        "{params.samtools} sort -@ {threads} -n -o {output.bam} {output.bam}.unsorted.bam && "
        "rm -f {output.bam}.unsorted.bam"

rule salvage_summary:
    input:
        expand("consensus/salvage/{id}_{chr}.salvage_stats.tsv", id=config['samples']['id'], chr=CHROMS)
    output:
        "consensus/salvage_summary.tsv"
    benchmark:
        "Benchmarks/consensus_fasta.salvage_summary.txt"
    run:
        header, rows = None, []
        for f in input:
            with open(f) as fh:
                header = next(fh).rstrip("\n").split("\t")
                rows.append(next(fh).rstrip("\n").split("\t"))
        with open(output[0], "w") as out:
            out.write("\t".join(header) + "\n")
            for r in rows:
                out.write("\t".join(r) + "\n")
            totals = [str(sum(int(r[i]) for r in rows)) for i in range(1, len(header))]
            out.write("\t".join(["all"] + totals) + "\n")


rule get_consensus_fasta:
    input:
        bam=CONSENSUS_NAMESORT_BAM,
        chr_dict = "Align/samtools_idx.txt"     
    output:
        fasta = "consensus/tmp/{chr}/{id}_{chr}_{split_idx}.fasta",
        stats = "consensus/tmp/{chr}/{id}_{chr}_{split_idx}.assembly_counts.tsv"
    threads: 1
    resources:
        mem_mb = 4000
    params:
        python = config['params']['general_python'],
        src_dir = config['params']['src_dir'], 
        ref = config['params']['ref_fa'],
        samtools = config['params'].get('samtools', 'samtools'),
        min_reads = config['params']['reads_per_BC'],
        downsample_ratio = config['params'].get('downsample_ratio', 1.0),
        use_samtools_reference = config['params'].get('use_samtools_reference', False),
        drop_intronic_mono_exon = config['params'].get('consensus_drop_intronic_mono_exon', False),
        consensus_gtf = config['params'].get('consensus_gtf', ''),
        intronic_frac = config['params'].get('consensus_intronic_frac', 0.5)
    benchmark:
        "Benchmarks/consensus_fasta.get_consensus_fasta.{id}.{chr}.{split_idx}.txt"
    run:
        command = ["{params.python}",
                    "{params.src_dir}/modules/clfr/consensus_fasta/consensus_fasta.py",
                    "--bam {input.bam}",
                    "--ref_fasta {params.ref}",
                    "--output_fasta {output.fasta}",
                    "--chrom {wildcards.chr}",
                    "--dict_file {input.chr_dict}",
                    "--split_index {wildcards.split_idx}",
                    "--min_reads {params.min_reads}",
                    "--downsample_ratio {params.downsample_ratio}",
                    "--stats_file {output.stats}",
                    "--samtools {params.samtools}"
                    ] 
        if params.use_samtools_reference:
            command.append("--use_samtools_reference")
        if params.drop_intronic_mono_exon:
            command.append("--drop_intronic_mono_exon --intronic_frac {params.intronic_frac}")
            if params.consensus_gtf:
                command.append("--gtf {params.consensus_gtf}")
        shell(" ".join(command))


split_cnt = list(range(NUM_SPLITS_CONSENSUS))  # [0, 1, 2, ..., 19]
MIN_FRAG_LEN = 400

rule fix_consensus_format:
    input:
        "consensus/tmp/{chr}/{id}_{chr}_{i}.fasta"
    output:
        "consensus/tmp/{chr}/{id}_{chr}_{i}.noN.fix.fasta"
    params:
        python = config['params']['general_python'],
        src_dir = config['params']['src_dir'],
        min_frag_len = MIN_FRAG_LEN
    benchmark:
        "Benchmarks/consensus_fasta.fix_consensus_format.{id}.{chr}.{i}.txt"
    shell:
        "{params.python} {params.src_dir}/modules/clfr/consensus_fasta/consensus_fasta_supp.py "
        "--module fix_fasta --input_fasta {input} --output_fasta {output}.tmp && "
        "awk -v min_frag_len={params.min_frag_len} "
        "'function emit() {{if (header != \"\" && length(seq) >= min_frag_len) print header \"\\n\" seq}} "
        "/^>/ {{emit(); header=$0; seq=\"\"; next}} {{seq=seq $0}} END {{emit()}}' "
        "{output}.tmp > {output} && rm -f {output}.tmp"

rule merge_consensus_fasta:
    input:
        fa = expand("consensus/tmp/{chr}/{id}_{chr}_{i}.noN.fix.fasta", id=config['samples']['id'],chr=CHROMS, i=split_cnt),
    output:
        "consensus/consensus.fasta"
    benchmark:
        "Benchmarks/consensus_fasta.merge_consensus_fasta.txt"
    shell:
        "cat {input.fa} > {output} "

rule fasta_frag_len_distribution_consensus:
    input:
        fasta = "consensus/consensus.fasta",
        stats = expand(
            "consensus/tmp/{chr}/{id}_{chr}_{i}.assembly_counts.tsv",
            id=config['samples']['id'], chr=CHROMS, i=split_cnt
        )
    output:
        png = "consensus/consensus_frag_length_distribution.png",
        stats = "consensus/consensus_assembly_counts.tsv"
        # "Align/frag_length_distribution.txt"
    params:
        python = config['params']['general_python'],
        src_dir = config['params']['src_dir'], 
        minreads_fasta = config['calc_frag']['minreads_fasta'], 
    benchmark:
        "Benchmarks/consensus_fasta.fasta_frag_len_distribution.txt"
    run:
        assembled_reads = 0
        assembled_umi_groups = 0
        with open(output.stats, "w") as merged_stats:
            merged_stats.write(
                "chrom\tsplit_index\tconsensus_assembled_reads\t"
                "consensus_assembled_umi_groups\n"
            )
            for stats_file in input.stats:
                with open(stats_file) as stats_handle:
                    header = next(stats_handle).rstrip("\n")
                    expected_header = (
                        "chrom\tsplit_index\tconsensus_assembled_reads\t"
                        "consensus_assembled_umi_groups"
                    )
                    if header != expected_header:
                        raise ValueError(f"Unexpected consensus stats header in {stats_file}: {header}")
                    line = next(stats_handle).rstrip("\n")
                    chrom, split_idx, reads, umi_groups = line.split("\t")
                    assembled_reads += int(reads)
                    assembled_umi_groups += int(umi_groups)
            merged_stats.write(
                f"all\tall\t{assembled_reads}\t{assembled_umi_groups}\n"
            )
        print(f"Consensus assembled reads: {assembled_reads}")
        print(f"Consensus assembled UMI groups: {assembled_umi_groups}")
        command = ["{params.python}",
                    "{params.src_dir}/modules/clfr/consensus_fasta/exon2fasta.py",
                    "--fasta {input.fasta}",
                    "--outdir consensus/",
                    "--name consensus",
                    "--minreads_fasta {params.minreads_fasta}",
                    "--counts_file {output.stats}",
                    "--module metrics_basic"] 
        shell(" ".join(command))

rule consensus_summary_report:
    input:
        counts = "consensus/consensus_assembly_counts.tsv",
        idx = "Align/samtools_idx.txt",
        split_log = "split_stat_read1.log"
    output:
        "consensus/summary_report.txt"
    params:
        python = config['params']['general_python'],
        src_dir = config['params']['src_dir'],
        cutoff = config['params']['reads_per_BC']
    benchmark:
        "Benchmarks/consensus_fasta.consensus_summary_report.txt"
    shell:
        """
        set -euo pipefail
        tmpdir=$(mktemp -d)
        {params.python} {params.src_dir}/modules/shared/metrics/count_umi.py {input.split_log} \\
            --cutoff {params.cutoff} --outdir $tmpdir
        {{
            echo "## consensus_assembly_counts"
            cat {input.counts}
            echo
            echo "## total mapped reads"
            awk '{{n+=$3}} END{{print n}}' {input.idx}
            echo
            echo "## reads/bc bin distribution (cutoff N>={params.cutoff})"
            cat $tmpdir/N{params.cutoff}.per_umi_read_count.bins.tsv
            echo
            echo "## per-bc read count stats (N>={params.cutoff})"
            cat $tmpdir/N{params.cutoff}.per_umi_read_count.stats.tsv
        }} > {output}
        rm -rf $tmpdir
        """

rule map_fasta_consensus:
    input:
        "consensus/consensus.fasta"
    output:
        "consensus/consensus.paf"
    threads:
        config['threads'].get('consensus_minimap', 20)
    params:
        minimap = config['params']['minimap2'],
        refgenome = config['params'].get('minimap2_mrna_index') or config['params']['ref_fa_mrna']
    benchmark:
        "Benchmarks/consensus_fasta.map_fasta_consensus.txt"
    run:
        command = ["{params.minimap} -x asm20 -t {threads} ",
                   "{params.refgenome} ",
                   "{input} > {output} "
                   ] 
        shell(" ".join(command))

rule correc_direction_consensus:
    input:
        consensus_fa = "consensus/consensus.fasta",
        consensus_paf = "consensus/consensus.paf"
    output:
        "consensus/consensus.fixRC.fasta"
    params:
        python = config['params']['general_python'],
        flank_end = config['frag_de_novo']['flank_end'],
        adapter_seq = config['frag_de_novo']['adapter_seq'],
        src_dir = config['params']['src_dir']
    benchmark:
        "Benchmarks/consensus_fasta.correc_direction_consensus.txt"
    run:
        if config['params']['correct_rc']==True:
            command = ["{params.python} ",
                    "{params.src_dir}/modules/clfr/denovo/denovo_supp.py ",
                    "--adapter_seq {params.adapter_seq} ",
                    "--flank_end {params.flank_end} ",
                    "--fasta {input.consensus_fa}",
                    "--module fix_fa_rc ",
                    "--outdir ./ "
                    ] 
            shell(" ".join(command))
        else:
            shell("cd consensus && ln -s consensus.fasta consensus.fixRC.fasta ")

# rule eval_Mandalorion:
#     input:
#         "consensus/consensus.fixRC.fasta"
#     output:
#         "Align/tmp/mando/Isoforms.filtered.clean.gtf"
#     params:
#         gtf = config['params']['mandalorion_gtf'],
#         ref = config['params']['mandalorion_ref'],
#         src = config['params']['mandalorion'],
#         python_exec = config['params']['mandalorion_python'],
#         cpus = config['threads']['metrics']
#     shell:
#         """
#         {params.python_exec} {params.src} -g {params.gtf} -G {params.ref} -f {input} -R 1 -t {params.cpus} -p Align/tmp/mando
#         """




# rule eval_SQANTI3:
#     input:
#         "Align/tmp/mando/Isoforms.filtered.clean.gtf"
#     output:
#         "Align/tmp/SQANTI3_QC_output/data_SQANTI3_report.pdf"
#     params:
#         python_exec = config['params']['sqanti3_python'],
#         sqanti3_script = config['params']['sqanti3_script'],
#         ref_gtf = config['params']['mandalorion_gtf'],
#         ref_fasta = config['params']['mandalorion_ref'],
#         cage_peak_bed = config['params']['sqanti3_cage_peak'],
#         polyA_motif_list = config['params']['sqanti3_polyA_motif'],
#         cpus = config['threads']['metrics']
#     shell:
#         """
#         {params.python_exec} {params.sqanti3_script} \
#             --isoforms {input} \
#             --refGTF {params.ref_gtf} \
#             --refFasta {params.ref_fasta} \
#             --polyA_motif_list {params.polyA_motif_list} \
#             --CAGE_peak {params.cage_peak_bed} \
#             --force_id_ignore \
#             --skipORF \
#             -o data \
#             --cpus {params.cpus} \
#             -d Align/tmp/SQANTI3_QC_output \
#             --report both
#         """

rule report_consensus:
    input:
        "Align/tmp/SQANTI3_QC_output/data_SQANTI3_report.pdf"
    output:
        "consensus/consensus.fixRC_SQANTI3_report.pdf"
    run:
        shell("cp {input} {output} ")
