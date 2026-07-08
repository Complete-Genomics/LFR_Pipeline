from pathlib import Path
import re
import glob
complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}

SEQUENCE_TYPE = config['params']['sequence_type'].lower()
MRNA_MAPPER = config['params']['mrna_mapper'].lower()


# We've gotten some inconsistent fastq names in the past
# This code scrapes the directories for a file matching "_1.fq.gz"
# This becomes the input for aggregating the reads
def get_fastqs_one(wildcards):
    # get fq base path from config
    fq_path = Path(config['samples']['fq_path'])
    # get lanes to add to the path
    lanes = config['samples']['lanes']
    fq_files = []
    for lane in lanes:
        # set aggregate path for each lane
        loc_path = fq_path / lane
        # iterate through the files in the filepath
        for fp in glob.glob(f"{loc_path}/*_read_1.fq.gz"):
            if "discard" not in str(fp):
                # print(fp)
                fq_files.append(str(fp))
    return fq_files            


def get_fastqs_two(wildcards):
    # get fq base path from config
    fq_path = Path(config['samples']['fq_path'])
    # get lanes to add to the path
    lanes = config['samples']['lanes']
    fq_files= []
    for lane in lanes:
        # set aggregate path for each lane
        loc_path = fq_path / lane
        # iterate through the files in the filepath
        
        if SEQUENCE_TYPE=='pe':
            for fp in glob.glob(f"{loc_path}/*_read_2.fq.gz"):
                if "discard" not in str(fp):
                    # print(fp)
                    fq_files.append(str(fp))
        elif SEQUENCE_TYPE=='se':
            for fp in glob.glob(f"{loc_path}/*_read*.fq.gz"):
                if "discard" not in str(fp):
                    # print(fp)
                    fq_files.append(str(fp))
    return fq_files            


# aggregate all fq1 files in data/
rule cat_read_one:
    input:
        get_fastqs_one
    output:
        "data/read_1.fq.gz"
    run:
        fq_path = Path(config['samples']['fq_path'])
        lane = config['samples']['lanes']
        if SEQUENCE_TYPE=='pe':
            # if their are multiple files concatenate them
            if len(lane) > 1:
                shell("cat {input} > {output}")
            # otherwise just link the one path
            else:
                shell("ln -s ../{input} {output}")
        elif SEQUENCE_TYPE=='se':
            shell("touch data/read_1.fq.gz")


# aggregate all fq2 files in data/
rule cat_read_two:
    input:
        get_fastqs_two
    output:
        "data/read_2.fq.gz"
    run:
        fq_path = Path(config['samples']['fq_path'])
        lane = config['samples']['lanes']
        if SEQUENCE_TYPE=='se':
            shell("ln -s ../{input} {output}")
        elif SEQUENCE_TYPE=='pe':
            # if their are multiple files concatenate them
            if len(lane) > 1:
                shell("cat {input} > {output}")
            # otherwise just link the one path
            else:
                shell("ln -s ../{input} {output}")
        

# barcodes can also get messy
# sometimes they're from the barcode list and sometimes they're from the RC list
# This function samples the barcodes and attempts to determine the appropriate list
def determine_barcode_list(n_samples):
    import gzip
    import sys
    base_swap = config['params']['bc_condition'].lower()
    _dir= config['params']['toolsdir'] 
    if base_swap=='standard_gc_swap':
        # GC dye swapped
        barcode_file = f"{_dir}/barcode.GCswap.list"
        barcode_rc_file = f"{_dir}/barcode_RC.GCswap.list"
    elif base_swap=='standard_ac_swap':
        barcode_file = f"{_dir}/barcode.ACswap.list"
        barcode_rc_file = f"{_dir}/barcode_RC.ACswap.list"
    else:
        barcode_file = f"{_dir}/barcode.list"
        barcode_rc_file = f"{_dir}/barcode_RC.list"
        # set potential barcode files based on the config file, also set the fastq path
        # barcode_file = config['params']['toolsdir'] + config['params']['barcode']
        # barcode_rc_file = config['params']['toolsdir'] + config['params']['barcode_RC']

    fastq_path = "data/read_2.fq.gz"

    # Read in barcodes fro, the barcodes file
    # We allow for one mismatch for every barcode
    # so we add those as well
    def get_barcodes(barcodes_file):
        print(f"reading in {barcodes_file}", file=sys.stderr)
        nucs = ['A', 'C', 'G', 'T']
        barcodes = []
        with open(barcodes_file, "r") as bcs_file:
            for line in bcs_file:
                bc = line.strip().split()[0]
                # iterate through each nucleotide in the barcode
                for i in range(0, len(bc)):
                    # iterate through potential mismatch nucleotides
                    for nuc in nucs:
                        bc_alt = list(bc)
                        bc_alt[i] = nuc
                        # add mismatched barcode to the barcodes list
                        barcodes.append("".join(bc_alt))
            
            # return the list as a set, their shouldn't be any duplicates
            # this is mostly to make the search faster
            barcodes_set = set(barcodes)
            return barcodes_set

    
    # scrape the fastq file for barcodes
    def get_fq_barcodes(fastq_path, n_samples):
        
        with gzip.open(fastq_path, "rb") as fastq:
            fq_bcs = []
            counter = 1
            print("Getting fastq barcodes", file=sys.stderr)
            for line in fastq:
                # keep track of number of barcodes we've looked at
                counter += 1
                # break after we've sampled enough barcodes
                if counter == n_samples + 1:
                    break
                if counter % 400000 == 0:
                    print(f"Read in {counter/4} lines", file=sys.stderr)
                # append barcodes for each read
                if counter % 4 == 3:
                    seq = line.decode('utf-8').strip()
                    bc_start_idx = config['params']['bc_start']-1
                    _bc_len, _bc_gap = 10, 6
                    bc1=seq[bc_start_idx: bc_start_idx+_bc_len]
                    bc2=seq[bc_start_idx+_bc_len+_bc_gap: bc_start_idx+_bc_gap+_bc_len*2]
                    bc3=seq[bc_start_idx+_bc_gap*2+_bc_len*2: bc_start_idx+_bc_gap*2+_bc_len*3]
                    fq_bcs.append(bc1)
                    fq_bcs.append(bc2)
                    fq_bcs.append(bc3)
            return fq_bcs

    # get barcodes, rc_barcodes, and fq_barcodes
    barcodes = get_barcodes(barcode_file)
    barcodes_rc = get_barcodes(barcode_rc_file)
    fq_bcs = get_fq_barcodes(fastq_path, n_samples)

    bcs_found_lst = []
    rc_bcs_found_lst = []
    # check fq barcodes against the barcodes and rc_barcodes
    
    for bc in fq_bcs:
        if bc in barcodes:
            bcs_found_lst.append(1)
        else:
            bcs_found_lst.append(0)
        if bc in barcodes_rc:
            rc_bcs_found_lst.append(1)
        else:
            rc_bcs_found_lst.append(0)
    bcs_found = sum(bcs_found_lst)
    rc_bcs_found = sum(rc_bcs_found_lst)
    print(f"Barcodes: {bcs_found}\nRC_Barcodes: {rc_bcs_found}\n",file=sys.stderr)
    # if more barcodes were found than RC barcodes return the barcodes path
    if bcs_found>rc_bcs_found:
        print("Barcode_list: {}".format(config['params']['barcode']), file=sys.stderr)
        return f"{_dir}/barcode.list"
    # if more RC barcodes found retunr RC barcodes path
    elif rc_bcs_found>bcs_found:
        print("Barcode_list: {}".format(config['params']['barcode_RC']), file=sys.stderr)
        return f"{_dir}/barcode_RC.list"
    # if the same number of barcodes were found try again with more samples (unlikely)
    # if we sample more than 10 mil reads, exit - something is wrong
    else:
        if n_samples > 4000000:
            print("Can't determine correct barcode file for sample", file=sys.stderr)
            sys.exit(1)
        else:
            return determine_barcode_list(n_samples * 2)

        

# rule for splitting reads
rule split_reads:
    input:
        expand("data/read_{i}.fq.gz", i=range(1,3))
    output:
        expand("data/split_read.{i}.fq.gz", i=range(1,3)),
        "split_stat_read1.log"
    params:
        r2_len = config['params']['read_len'],
        r1_len = config['params']['read_len_r1'], 
        toolsdir = config['params']['toolsdir'],
        src_dir = config['params']['src_dir'],
        # randomBC_dir = config['params']['randomBC_dir'],
        general_python = config['params']['general_python'],
        bc_len = config['params']['sbc_len'],
        bc_len_redundant = config['params']['cbc_len'],
        cbc_len = config['params']['cbc_len'],
        bc_start = config['params']['bc_start'],
        gdna_start = config['params']['gdna_start'],
        additional_bc_start = config['params']['additional_bc_start'],
        additional_bc_len = config['params']['additional_bc_len'],
        gdna_start_r1 = config['params']['gdna_start_r1'], 
        additional_bc_len_r1 = config['params']['additional_bc_len_r1'],
        adapter_len = config['params']['adapter_len'], 
    threads:
        config['threads']['split_group']
    run:
        if config['modules']['BCsplit']==True:
            if config['params']['bc_condition'] == 'random_bc_umi_rc':
                shell("{params.general_python} {params.src_dir}/modules/shared/splitreads/split_barcode_LFR.py "
                "--r1 {input[0]} --r2 {input[1]} --read_len {params.r2_len} --output data/split_read "
                "--cbc_len {params.cbc_len} --bc_len_redundant {params.bc_len_redundant} "
                "--bc_start {params.bc_start} --gdna_start {params.gdna_start} "
                "--additional_bc_start {params.additional_bc_start} --additional_bc_len {params.additional_bc_len} "
                "--gdna_start_r1 {params.gdna_start_r1} --read_len_r1 {params.r1_len} "
                "--additional_bc_len_r1 {params.additional_bc_len_r1} --reverse_complement --threads {threads} "
                "2> data/split_stat_read.err")
            elif config['params']['bc_condition'] == 'random_bc':
                shell("{params.general_python} {params.src_dir}/modules/shared/splitreads/split_barcode_LFR.py "
                "--r1 {input[0]} --r2 {input[1]} --read_len {params.r2_len} --output data/split_read "
                "--cbc_len {params.cbc_len} --bc_len_redundant {params.bc_len_redundant} "
                "--bc_start {params.bc_start} --gdna_start {params.gdna_start} "
                "--additional_bc_start {params.additional_bc_start} --additional_bc_len {params.additional_bc_len} "
                "--gdna_start_r1 {params.gdna_start_r1} --read_len_r1 {params.r1_len} "
                "--additional_bc_len_r1 {params.additional_bc_len_r1} --threads {threads} "
                "2> data/split_stat_read.err")
            elif config['params']['bc_condition'] == 'standard':
                # determine which barcode list to use
                if config['params']['enforced_bc_list']=='barcode':
                    params.barcode = '{params.toolsdir}/barcode.list'
                elif config['params']['enforced_bc_list']=='barcode_RC':
                    params.barcode = '{params.toolsdir}/barcode_RC.list'
                else:
                    params.barcode = determine_barcode_list(400000)

                if params.bc_len != 42:
                    raise ValueError("Only 42 bp stLFR barcodes are supported; bc_len == 30 is retired.")
                params.barcode = determine_barcode_list(400000)
                shell("{params.general_python} {params.src_dir}/modules/shared/splitreads/split_barcode_LFR.py "
                    "--barcode {params.barcode} "
                    "--r1 {input[0]} --r2 {input[1]} --read_len {params.r2_len} --output data/split_read "
                    "--bc_start {params.bc_start} --gdna_start {params.gdna_start} "
                    "--additional_bc_start {params.additional_bc_start} --additional_bc_len {params.additional_bc_len} "
                    "--gdna_start_r1 {params.gdna_start_r1} --read_len_r1 {params.r1_len} "
                    "--swap none --output_mode separate --threads {threads} "
                    "2> data/split_stat_read.err")
            elif config['params']['bc_condition'] == 'BCgDNA':
                params.barcode = determine_barcode_list(400000)
                shell("{params.general_python} {params.src_dir}/modules/shared/splitreads/split_barcode_LFR.py "
                        "--barcode {params.barcode} "
                        "--r1 {input[0]} --r2 {input[1]} --read_len {params.r2_len} --output data/split_read "
                        "--bc_start {params.bc_start} --gdna_start {params.gdna_start} "
                        "--additional_bc_start {params.additional_bc_start} --additional_bc_len {params.additional_bc_len} "
                        "--gdna_start_r1 {params.gdna_start_r1} --read_len_r1 {params.r1_len} "
                        "--adapter_len {params.adapter_len} --swap none --output_mode stratified --threads {threads} "
                        "2> data/split_stat_read.err")
            elif config['params']['bc_condition'].lower() == 'standard_ac_swap':
                params.barcode = determine_barcode_list(400000)
                shell("{params.general_python} {params.src_dir}/modules/shared/splitreads/split_barcode_LFR.py "
                        "--barcode {params.barcode} "
                        "--r1 {input[0]} --r2 {input[1]} --read_len {params.r2_len} --output data/split_read "
                        "--bc_start {params.bc_start} --gdna_start {params.gdna_start} "
                        "--additional_bc_start {params.additional_bc_start} --additional_bc_len {params.additional_bc_len} "
                        "--gdna_start_r1 {params.gdna_start_r1} --read_len_r1 {params.r1_len} "
                        "--swap ac --output_mode single --threads {threads} "
                        "2> data/split_stat_read.err")
            elif config['params']['bc_condition'].lower() == 'standard_gc_swap':
                params.barcode = determine_barcode_list(400000)
                shell("{params.general_python} {params.src_dir}/modules/shared/splitreads/split_barcode_LFR.py "
                        "--barcode {params.barcode} "
                        "--r1 {input[0]} --r2 {input[1]} --read_len {params.r2_len} --output data/split_read "
                        "--bc_start {params.bc_start} --gdna_start {params.gdna_start} "
                        "--additional_bc_start {params.additional_bc_start} --additional_bc_len {params.additional_bc_len} "
                        "--gdna_start_r1 {params.gdna_start_r1} --read_len_r1 {params.r1_len} "
                        "--swap gc --output_mode single --threads {threads} "
                        "2> data/split_stat_read.err")
            else:
                print("unknown type")
        else:
            shell("touch split_stat_read1.log && cd data && ln -s read_1.fq.gz split_read.1.fq.gz && ln -s read_2.fq.gz split_read.2.fq.gz && cd ..")

rule trim_reads:
    input:
        f1="data/split_read.1.fq.gz",
        f2="data/split_read.2.fq.gz"
    output:
        t1= "data/split_read_1_trimmed.fastq.gz",
        t2= "data/split_read_2_trimmed.fastq.gz"
    params:
        bbduk = config['params']['bbduk'],
        sequence_type= config['params']['sequence_type'].lower()
    shell:
        """
        if [[ "{params.sequence_type}" == "pe" ]]; then
            {params.bbduk} in1={input.f1} in2={input.f2} out1={output.t1} out2={output.t2} qtrim=rl 
        elif [[ "{params.sequence_type}" == "se" ]]; then
            {params.bbduk} in={input.f2} out={output.t2} qtrim=rl && touch data/split_read_1_trimmed.fastq.gz
        else
            echo "Unknown type {params.sequence_type}" >&2;
            exit 1;
        fi
        """
