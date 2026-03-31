rule reformat_tag:
    input:
        "Align/data.sort.removedup_rm000.bam"
    output:
        "Align/data.sort.removedup_rm000.tag.bam"
    params:
        jvarkit_samjdk = config['params']['jvarkit_samjdk']
    run:
        if config.get('params', {}).get('bc_condition') == 'combined':
            shell("java -jar {params.jvarkit_samjdk} "
                  "--samoutputformat BAM -e "
                  "'final String readId = record.getReadName(); "
                  "final Object mi = record.getAttribute(\"MI\"); "
                  "if (mi != null) {{ String newReadId = readId + \"_\" + mi.toString(); record.setReadName(newReadId); }} "
                  "return record;' "
                  "{input} > {output} 2> reformat_tag.log")
        else:
            shell("cd Align && ln -s data.sort.removedup_rm000.bam data.sort.removedup_rm000.tag.bam")

rule index_tag_bam:
    input:
        "Align/data.sort.removedup_rm000.tag.bam"
    output:
        "Align/data.sort.removedup_rm000.tag.bam.bai"
    shell:
        "samtools index -@ 10 {input}"

rule get_bed:
    input:
        bam= "Align/data.sort.removedup_rm000.tag.bam",
        bai = "Align/data.sort.removedup_rm000.tag.bam.bai",
    output:
        "Align/{chr}_bed"
    group: 
        "group_10"
    threads: 
        1
    shell:
        "samtools view -b -F 0x900 {input.bam} {wildcards.chr} | bedtools bamtobed -split -i - > {output} "
        

rule calc_frag_bed:
    input:
        "Align/{chr}_bed"
    output:
        "Align/{chr}.chr2bx"
    threads:
        config['threads']['calc_frag']
    params:
        min_frag = MIN_FRAG,
        read_len = config['params']['read_len'],
        n_tolerance = config['calc_frag']['n_tolerance'],
        toolsdir = config['params']['toolsdir'],
        include_dups = config['calc_frag']['include_dups'],
        minreads = config['calc_frag']['minreads_fasta'],
        src_dir = config['params']['src_dir'],
        python = config['params']['general_python'],
        # additional_bc_len = config['params']['additional_bc_len'],
        # bc_condition = config['params']['bc_condition'],
    run:
        command = ["{params.python}",
                   "{params.src_dir}/modules/clfr/filterN100/calc_frag_len.bed.py",
                   "--minfrag {params.min_frag}",
                   "--splitdist 50000",
                   "--readlen {params.read_len}",
                   "--chroms {wildcards.chr} ",
                   "--threads {threads}",
                   "--minreads {params.minreads}",
                   "--n_tolerance {params.n_tolerance}", 
                #    "--additional_bc_len {params.additional_bc_len}",
                #    "--bc_condition {params.bc_condition}",
                   "--outdir Calc_Frag_Length_50000"] 
        shell(" ".join(command))

rule merge_by_bx:
    input:
        "Align/{chr}.chr2bx"
    output:
        "Align/{chr}.N100_merged"
    shell:
        "bedtools sort -i {input} | bedtools merge -s -c 6 -o distinct -i - > {output} "

rule bc_bed:
    input:
        "Align/{chr}.N100_merged"
    output:
        "Align/{chr}_bc.bed"
    shell:
        "bedtools sort -i {input} | bedtools merge -s -c 6 -o distinct -i - > {output} "
  
       
rule replace_BX_with_chr:
    input:
        inputfile="Align/{chr}.N100_merged"
    output:
        outputfile="Align/{chr}.N100_intChr"
    params:
        chr="{chr}"
    run:
        outfile = open(output.outputfile, 'w')
        with open(input.inputfile, 'r') as f:
            for line in f:
                info = line.strip().split('\t')
                bx = info[0]
                ## need to be BED6 format
                new_line = '\t'.join([params.chr]+ info[1:-1]+[bx, info[-1]]+[info[-1]])
                outfile.write(new_line+'\n')
        outfile.close()

rule bed_fasta:
    input:
        "Align/{chr}.N100_intChr"
    output:
        "Align/{chr}.fa"
    params:
        ref=REF
    shell:
        "bedtools getfasta -fi {params.ref} -bed {input} -fo {output} -name -s "


rule merge_fasta:
    input:
        fa = expand("Align/{chr}.fa", chr=CHROMS),
    output:
        "Align/all_tmp.fasta"
    shell:
        "cat {input.fa} > {output} "
        
rule bc_reads_dict:
    input:
        "split_stat_read1.log"
    output:
        "Align/bc_reads_dict"
    params:
        python = config['params']['general_python'],
        src_dir = config['params']['src_dir'],
        minreads_fasta = config['calc_frag']['minreads_fasta']
    run:
        command = ["{params.python} ",
                    "{params.src_dir}/modules/clfr/filterN100/filterN100.latest.py ",
                    "--minreads_fasta {params.minreads_fasta} "
                    "--module filter_2bc_count"] 
        shell(" ".join(command))   

rule merge_exon:
    input:
        "Align/all_tmp.fasta",
        "Align/bc_reads_dict"
    output:
        "Align/all_N100.fasta"
    params:
        python = config['params']['general_python'],
        src_dir = config['params']['src_dir'],
    run:
        command = ["{params.python}",
                    "{params.src_dir}/modules/clfr/filterN100/filterN100.latest.py",
                    "--module merge_exon"] 
        shell(" ".join(command))

rule plot_frag_len_distribution:
    input:
        "Align/all_N100.fasta"
    output:
        "Align/frag_length_distribution_N100.pdf"
    params:
        python = config['params']['general_python'],
        src_dir = config['params']['src_dir'], 
        minreads_fasta = config['calc_frag']['minreads_fasta'], 
    run:
        command = ["{params.python}",
                    "{params.src_dir}/modules/clfr/rna_16s/rna_16s.py",
                    "--fasta {input}",
                    "--outdir Align/",
                    "--minreads_fasta {params.minreads_fasta}",
                    "--module metrics_basic"] 
        shell(" ".join(command))

# plot per_base_density, +'ERCC-00130', 'ERCC-00096', 'ERCC-00116', 'ERCC-00025'
rule plot_frag_base_coverage:
    input:
        "Align/all_N100.fasta"
    output:
        "Align/frag_coverage_done"
    params:
        python = config['params']['general_python'],
        src_dir = config['params']['src_dir'],
    run:
        command = ["{params.python}",
                    "{params.src_dir}/modules/clfr/rna_16s/rna_16s.py",
                    "--chrom chr21",
                    "--module frag_coverage"] 
        shell(" ".join(command))