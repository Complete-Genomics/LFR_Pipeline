from rna_16s import run_write2fq_index
import subprocess
import argparse
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--dirname", type=str)
    parser.add_argument("--r1r2", type=int)
    parser.add_argument("--bc_list", type=str)
    args = parser.parse_args() 
    dirname = args.dirname
    r1r2 = args.r1r2
    bc_list = args.bc_list

    n=3
    ### get read/bc distribution
    # reads_per_bc_distribution_regardless_mapping(dirname)

    ### get bc.2.fq.gz
    # bc_list ='bc4denovo_n.txt'
    run_write2fq_index(r1r2, bc_list)
    subprocess.call(f'touch fq{r1r2}_done', shell=True)
