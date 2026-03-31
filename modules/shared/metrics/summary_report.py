#!/usr/bin/env python3
import subprocess
import sys
import yaml
import re
import os

def get_split_read_log(log_path):
    '''
    parse split log for pertinent output
    '''
    with open(log_path, "r") as log:
        for i in range(0,4):
            line = log.readline()
            if i == 1:
                real_bc = line.strip().split()[2]
            if i == 2:
                read_pair = line.strip().split()[2]
            if i == 3:
                per_split = line.strip().split()[4][1:] + "%"
                
    return real_bc, read_pair, per_split


def get_picard_metrics(log_path):
    '''
    parse picard metrics for pertinent output
    '''
    with open(log_path, "r") as metrics:
        for line in metrics:
            if line.startswith("Mapping rate"):
                map_rate = line.strip().split()[2]
            if line.startswith("Duplicate rate"):
                dup_rate = line.strip().split()[2]
                
    return map_rate, dup_rate
                

def get_coverage_depth(log_path):
    '''
    parse coverage depth for pertinent output
    '''
    with open(log_path, "r") as coverage:
        for i in range(0,2):
            line = coverage.readline()
            if i == 0:
                dep = line.strip().split()[3]
            if i == 1:
                cov = line.strip().split()[1]
                
    return dep, cov


def get_insert_size(log_path):
    '''
    parse insert size metrics for pertinent output
    '''
    with open(log_path, "r") as insert_size:
        for i in range(0,3):
            line = insert_size.readline()
            if i == 2:
                median = line.strip().split()[0]
                mean = line.strip().split()[4]
    
    return median, mean


def get_insert_size_gatk(log_path):
    '''
    parse a gatk insert size metrics file
    '''
    with open(log_path, "r") as insert_size:
        for i in range(0,8):
            line = insert_size.readline()
            if i == 7:
                median = line.strip().split()[0]
                mean = line.strip().split()[5]
    
    return median, mean


def parse_report(path):
    '''
    parse fragment calculations output for pertinent output
    '''
    with open(path, "r") as barcodes:
        frag_stats = []
        for line in barcodes:
            frag_stats.append(line.strip())
    return frag_stats


def get_hapcut_results(path):
    '''
    parse hapcut output for results
    '''
    with open(path, "r") as hapcut:
        hapcut_results = {}
        for line in hapcut:
            if not line.strip():
                continue
            parts = re.split(r'\s{2,}', line.strip())
            hapcut_results[parts[0]] = parts[1]

    return hapcut_results['phased count:'], hapcut_results['AN50:'], hapcut_results['N50:']
            
                
def get_longhap_results(path):
    '''
    parse longhap output for results
    '''
    with open(path, "r") as longhap:
        for line in longhap:
            if line.startswith("N50:"):
                n_fifty = line.strip().split(":")[1]
                n_fifty = n_fifty.strip()
            if line.startswith("AN50:"):
                an_fifty = line.strip().split(":")[1]
                an_fifty = an_fifty.strip()
            if line.startswith("phased_snp:"):
                phased = line.strip().split()[0].split(":")[1]
    
    return phased, an_fifty, n_fifty


if __name__ == "__main__":
    import sys
    import yaml
    
    try:
        # load config file to get sample id, read length, and split distances
        with open("config.yaml") as cnf_path:
            config = yaml.load(cnf_path, yaml.SafeLoader)
            sample_id = config['samples']['id']
            read_len = config['params']['read_len']
            try:
                # if r1,r2 with different length
                read1_len = config['params']['read1_len']
            except:
                read1_len = read_len
            split_dist = config['calc_frag']['split_dist']
    except:
        print("config.yaml not found", file=sys.stderr)

    
    # define expected paths
    split_read_log = "split_stat_read1.log"
    picard_metrics = "Align/picard_align_metrics.txt"
    # picard_metrics2 = "Align/picard_align_metrics_BCreads_only.txt"
    coverage_depth = "Align/coverage_depth.txt"
    longhap = "Make_Vcf/step4_longhap/longhap_results.txt"
    hapcut = "Make_Vcf/step3_hapcut/step4_compare_with_refphasing/hapcut_eval.txt"
    
    
    try:
        # attempt to parse slit read log and output results
        real_bc, read_pair, per_split = get_split_read_log(split_read_log)
        print("{}\t{}".format("Sample ID:", sample_id))
        print("{}\t{}".format("Barcode Count:", real_bc))
        print("{}\t{}".format("Read Pairs:", read_pair))
        print("{}\t{}".format("Total Bases:", int(read_pair)*(read_len+read1_len)))
        print("{}\t{}".format("Barcode Split Rate:", per_split))
    except:
        print("Couldn't get split log stats", file=sys.stderr)

    try:
        # attempt to parse picard metrics and coverage depth
        # output results
        map_rate, dup_rate = get_picard_metrics(picard_metrics)
        dep, cov = get_coverage_depth(coverage_depth)
        print("{}\t{}".format("Mapping Rate:", map_rate))
        print("{}\t{}".format("Coverage:", cov))
        print("{}\t{}".format("Depth:", dep))
        print("{}\t{}".format("Duplicate Rate(BC-independent):", ""))
        print("{}\t{}".format("Duplicate Rate(BC-dependent):", dup_rate))
    except:
        print("Couldn't get mapping and depth stats", file=sys.stderr)
    # try:
    #     # attempt to parse picard metrics and coverage depth
    #     # output results
    #     map_rate, dup_rate = get_picard_metrics(picard_metrics2)
    #     print("{}\t{}".format("Mapping Rate (BCreads):", map_rate))
    # except:
    #     print("Couldn't get mapping and depth stats BConly", file=sys.stderr)

    try:
        # attempt to parse insert size metrics
        insert_size_metrics = "Align/sentieon_is_{}_metric.txt".format(sample_id)
        median, mean = get_insert_size(insert_size_metrics)
        print("{}\t{}".format("Median Insert Size:", median))
        print("{}\t{}".format("Mean Insert Size:", mean))
    except:
        print("Couldn't get insert size stats, trying gatk file path", file=sys.stderr)
        try:
            # if the sentieon file doesn't work try gatk
            insert_size_metrics = "Align/gatk_metrics_data.insert_size_metrics"
            median, mean = get_insert_size_gatk(insert_size_metrics)
            print("{}\t{}".format("Median Insert Size:", median))
            print("{}\t{}".format("Mean Insert Size:", mean))
        except:
            print("no GATK median insert size stats")
            print("no GATK mean insert size stats")
            

    for dist in split_dist:
        # for all split distances
        # attempt to parse the appropriate frag_and_bc_summary file
        fragment_count = "Calc_Frag_Length_" + str(dist) + "/frag_and_bc_summary.txt"
        try:
            fragment_stats = parse_report(fragment_count)
            for i in range(13):
                print(fragment_stats[i])
        except:
            print("Couldn't get fragment stats", file=sys.stderr)

    try: 
        path = "Calc_Frag_Length_50000/reads_per_bc_bins.txt"
        reads_per_bins = parse_report(path)
        for i in range(len(reads_per_bins)):
            print(reads_per_bins[i])
    except:
        print("Couldn't get fragment stats", file=sys.stderr)

    try:
        # attempt to parse hapcut results
        phased_hap, an_fifty_hap, n_fifty_hap = get_hapcut_results(hapcut)
        print("{}\t{}".format("Hapcut AN50:", an_fifty_hap))
        print("{}\t{}".format("Hapcut N50:", n_fifty_hap))
        print("{}\t{}".format("Hapcut Phased Count:", phased_hap))
    except:
        print("Couldn't get hapcut results", file=sys.stderr)
    
    try:
        # attempt to parse longhap results
        phased_lh, an_fifty_lh, n_fifty_lh = get_longhap_results(longhap)
        print("{}\t{}".format("LongHap AN50:", an_fifty_lh))
        print("{}\t{}".format("LongHap N50:", n_fifty_lh))
        print("{}\t{}".format("LongHap Phased SNPs:", phased_lh))
    except:
        print("Couldn't get longhap results", file=sys.stderr)

    file = 'Align/mapped_uniq_bc_bases_count.txt'
    if os.path.exists(file):
        try:
            cmd= f"cat {file}"
            subprocess.call(cmd, shell=True)
        except:
            print("Couldn't get mapped_uniq_bc_bases_count", file=sys.stderr)
