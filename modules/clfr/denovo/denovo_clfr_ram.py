'''
intermediate file to RAM, reduce I/O

usage:
# as megahit needs >=2 cpu, 2*30=60 cpu
num_processes=30
python $src --module create_cmd --num_node 5
python $src --module denovo_parallel --num_processes ${num_processes} --n_line_chunk 2000000 --start_idx 5 --end_idx 200 --sequence_type pe --nth_of_nodes 0

qsub -cwd -l vf=200G,num_proc=${num_processes} -P P21Z18000N0016 -q mgi_supermem.q denovo.sh

TODO:
1/ bc list, 
2/ meta_data list filter 100
3/ multiple assembly

'''
import os, io
import subprocess
import multiprocessing as mp
from multiprocessing import Pool
import argparse
from collections import defaultdict
import itertools
import fcntl
import datetime
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("--num_processes", type=int, required=False)
parser.add_argument("--bc_list", type=str, required=False)
parser.add_argument("--sequence_type", type=str, required=False)
parser.add_argument("--k_min", type=int, required=False)
parser.add_argument("--k_max", type=int, required=False)
parser.add_argument("--min_ctg_len", type=int, required=False)
parser.add_argument("--n_line_chunk", type=int, required=False)
parser.add_argument("--start_idx", type=int, required=False)
parser.add_argument("--end_idx", type=int, required=False)
parser.add_argument("--module", type=str)
parser.add_argument("--num_node", type=int, required=False)
parser.add_argument("--nth_of_nodes", type=int, required=False)
parser.add_argument('--debug', action='store_true', default=False, help='Enable debug mode')
parser.add_argument("--megahit", type=str, default='megahit', help='Path to megahit binary')
parser.add_argument("--rg", type=str, default='rg', help='Path to rg (ripgrep) binary')


args = parser.parse_args()
num_processes = args.num_processes
# bc with >50 reads
bc_list = args.bc_list
sequence_type = args.sequence_type
K_MIN = args.k_min
K_MAX = args.k_max
MIN_CTG_LEN = args.min_ctg_len
n_line_chunk = args.n_line_chunk
num_node = args.num_node
module = args.module
start_idx = args.start_idx
end_idx = args.end_idx
ID = args.nth_of_nodes
MEGAHIT = args.megahit
RG = args.rg

def process_barcode_se(barcode, shared_meta_data2, lock):
    K_MIN = 41
    K_MAX = 41
    MIN_CTG_LEN = 400
    megahit = MEGAHIT
    rg = RG

    try:
        # Create in-memory files for R1 and R2 using io.BytesIO
        # r1_fasta = io.BytesIO()
        r2_fasta = io.BytesIO()
        # r1_fasta.writelines(f"{line}\n".encode() for line in shared_meta_data1.get(barcode))
        r2_fasta.writelines(f"{line}\n".encode() for line in shared_meta_data2.get(barcode))
        # r1_fasta.seek(0)
        r2_fasta.seek(0)
        
        num_cpu =2 # at least 2
        # Command for megahit
        megahit_command = (
            f"{megahit} -r /dev/stdin -t {num_cpu} "
            f"-o /dev/shm/{BATCH_LANE}_{ID}/{barcode} --out-prefix {barcode} --k-min {K_MIN} --k-max {K_MAX} --force "
            f"--min-contig-len={MIN_CTG_LEN}"
        )

        # Run megahit using pipes
        process = subprocess.Popen(
            megahit_command,
            shell=True,
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
        )
        stdout, stderr = process.communicate(input=r2_fasta.read())

        bc = f'>{barcode}'   
        # if DEBUG==True:
        #     command = f"sed '1s/^/{bc}/; 3s/^/{bc}/; 5s/^/{bc}/; 7s/^/{bc}/' /dev/shm/{BATCH_LANE}_{ID}/{barcode}/{barcode}.contigs.fa >> denovo/final_contigs_{ID}.fa "
        # else:
        command = f"sed '1s/^/{bc}/; 3s/^/{bc}/; 5s/^/{bc}/; 7s/^/{bc}/' /dev/shm/{BATCH_LANE}_{ID}/{barcode}/{barcode}.contigs.fa >> denovo/final_contigs_{ID}.fa && rm -r /dev/shm/{BATCH_LANE}_{ID}/{barcode} "
        with lock:
            res = subprocess.run(command, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    except Exception as e: print(e)
    return


def process_barcode_pe(barcode, shared_meta_data1, shared_meta_data2, lock):
    K_MIN = 41
    K_MAX = 41
    MIN_CTG_LEN = 400
    megahit = MEGAHIT
    rg = RG

    try:
        # Create in-memory files for R1 and R2 using io.BytesIO
        r1_fasta = io.BytesIO()
        r2_fasta = io.BytesIO()
        r1_fasta.writelines(f"{line}\n".encode() for line in shared_meta_data1.get(barcode))
        r2_fasta.writelines(f"{line}\n".encode() for line in shared_meta_data2.get(barcode))
        r1_fasta.seek(0)
        r2_fasta.seek(0)
        
        num_cpu =2 # at least 2
        # Command for megahit
        megahit_command = (
            f"{megahit} -1 /dev/stdin -2 /dev/stdin -t {num_cpu} "
            f"-o /dev/shm/{BATCH_LANE}_{ID}/{barcode} --out-prefix {barcode} --k-min {K_MIN} --k-max {K_MAX} --force "
            f"--min-contig-len={MIN_CTG_LEN}"
        )

        # Run megahit using pipes
        process = subprocess.Popen(
            megahit_command,
            shell=True,
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
        )
        stdout, stderr = process.communicate(input=r1_fasta.read() + r2_fasta.read())
        # Check if the megahit command executed successfully
        # if process.returncode != 0:
        #     print(f"Megahit failed: {stderr.decode()}")
        #     return
        # else:
        #     print(f"Megahit success: {stdout.decode()}")

        ## add bc to fasta header line
        bc = f'>{barcode}'
        
        command = f"sed '1s/^/{bc}/; 3s/^/{bc}/; 5s/^/{bc}/; 7s/^/{bc}/' /dev/shm/{BATCH_LANE}_{ID}/{barcode}/{barcode}.contigs.fa >> denovo/final_contigs_{ID}.fa && rm -r /dev/shm/{BATCH_LANE}_{ID}/{barcode} "
        with lock:
            res = subprocess.run(command, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            # Check if the sed command executed successfully
            # if res.returncode != 0:
            #     print(f"Sed command failed: {res.stderr.decode()}")
            # else:
            #     print(f"Sed command success: {res.stdout.decode()}")

    except Exception as e: print(e)
    return

def denovo_pe(n_line_chunk, start_idx):
    K_MIN = 41
    K_MAX = 41
    MIN_CTG_LEN = 400
    LOC = 'sj'
    bc_len = 15

    meta_data1 = defaultdict(list)
    meta_data2 = defaultdict(list)

    ## load reads to mem to speedup, *sgrep.tsv with diff len after trim, unable to use seek to get idx
    with open(f'denovo/{NAME}1_sgrep.tsv', 'r') as f:
        for line in itertools.islice(f, start_idx, start_idx+n_line_chunk):
            info = line.strip().split('\t')
            bc = info[0][5:5+bc_len]
            id = '>'+info[0][22:]
            seq = info[1]
            meta_data1[bc].append(id)
            meta_data1[bc].append(seq)

    with open(f'denovo/{NAME}2_sgrep.tsv', 'r') as f:
        for line in itertools.islice(f, start_idx, start_idx+n_line_chunk):
            info = line.strip().split('\t')
            bc = info[0][5:5+bc_len]
            id = '>'+info[0][22:]
            seq = info[1]
            meta_data2[bc].append(id)
            meta_data2[bc].append(seq)

    # subprocess.call(f'mkdir -p /dev/shm/{BATCH_LANE}_{ID}', shell=True)

    with mp.Manager() as manager:
        shared_meta_data1 = manager.dict(meta_data1)
        shared_meta_data2 = manager.dict(meta_data2)
        lock = manager.Lock()

        with mp.Pool(num_processes) as pool:
            pool.starmap(process_barcode_pe, [(barcode, shared_meta_data1, shared_meta_data2, lock) for barcode in meta_data2.keys()])

    print(f'start_idx={start_idx}')
    print(f'denovo_BC_counts={len(meta_data2)}')


    ## if no results/bc/bc.fa, bc folder will not be removed
    # try:
    #     cmd=f'rm -r /dev/shm/{BATCH_LANE}_{ID}/*'
    #     subprocess.call(cmd, shell=True)
    #     # cmd=f'rm -r denovo/splits_{ID}/*'
    #     # subprocess.call(cmd, shell=True)
    # except Exception as e: print(e)

    return

def denovo_se(n_line_chunk, start_idx):
    K_MIN = 41
    K_MAX = 41
    MIN_CTG_LEN = 400
    LOC = 'sj'
    bc_len = 15
    meta_data2 = defaultdict(list)

    ## load reads to mem to speedup, *sgrep.tsv with diff len after trim, unable to use seek to get idx
    with open(f'denovo/{NAME}2_sgrep.tsv', 'r') as f:
        for line in itertools.islice(f, start_idx, start_idx+n_line_chunk):
            info = line.strip().split('\t')
            bc = info[0][5:5+bc_len]
            id = '>'+info[0][22:]
            seq = info[1]
            meta_data2[bc].append(id)
            meta_data2[bc].append(seq)

    with mp.Manager() as manager:
        # shared_meta_data1 = manager.dict(meta_data1)
        shared_meta_data2 = manager.dict(meta_data2)
        lock = manager.Lock()

        with mp.Pool(num_processes) as pool:
            pool.starmap(process_barcode_se, [(barcode, shared_meta_data2, lock) for barcode in meta_data2.keys()])

    print(f'start_idx={start_idx}')
    print(f'denovo_BC_counts={len(meta_data2)}')

    return

def read_fraction_of_file(file_path, num_node):
    '''
    output cmd to run on different nodes
    --start_idx 0 --end_idx 142 --nth_of_nodes 0, next run on basm02
    --start_idx 142 --end_idx 284 --nth_of_nodes 1, next run on basm08
    '''
    with open(file_path, 'r') as f:
        total_lines = sum(1 for line in f)

    lines_each_chunk = total_lines // num_node

    for i in range(0, total_lines, lines_each_chunk):
        start_idx = i
        end_idx = min(i + lines_each_chunk, total_lines)
        nth_of_nodes = i//lines_each_chunk
        print(f'--start_idx {start_idx} --end_idx {end_idx} --nth_of_nodes {nth_of_nodes}')

    return 

def create_bins(start_idx, end_idx, bin_size):
    bins = []
    for i in range(start_idx, end_idx, bin_size):
        bins.append([i, min(i + bin_size, end_idx)])
    return bins

def parse_bc_fasta(in_fasta):
    bc_len = 15

    meta_data2 = defaultdict(list)
    with open(in_fasta, 'r') as f:
        for line in f:
            info = line.strip().split('\t')
            bc = info[0][5:5+bc_len]
            id = '>'+info[0][22:]
            seq = info[1]
            meta_data2[bc].append(id)
            meta_data2[bc].append(seq)

    for barcode in meta_data2.keys():
        rm = 'AAAAAAAAAAAAAAA'
        if barcode !=rm:
            with open(f'splits/{barcode}_R2.fasta', 'w') as f:
                for line in meta_data2.get(barcode):
                    f.write(f"{line}\n")

if __name__ == "__main__":
    NAME = 'data_R'
    LOC = 'sj'
    current_path = Path.cwd()
    info =str(current_path).split('/')
    BATCH_LANE= info[-3]+'_'+info[-2]
    if args.debug:
        DEBUG=True
    else:
        DEBUG=False

    if module =='denovo_parallel':
        print(f'start={datetime.datetime.now()}')
        # if not os.path.exists(f'denovo/splits_{ID}'):
        #     os.mkdir(f'denovo/splits_{ID}')
        # if not os.path.exists(f'denovo/{BATCH_LANE}_{ID}'):
        #     os.mkdir(f'denovo/{BATCH_LANE}_{ID}')

        subprocess.call(f'mkdir -p /dev/shm/{BATCH_LANE}_{ID}', shell=True)
        bins = create_bins(start_idx, end_idx, n_line_chunk)

        if sequence_type == 'pe':
            for start_idx_each_chunk,end_idx_each_chunk in bins:
                n_line_each_chunk = end_idx_each_chunk-start_idx_each_chunk
                # print(f'{start_idx_each_chunk} , {end_idx_each_chunk}')
                denovo_pe(n_line_each_chunk, start_idx_each_chunk)
        elif sequence_type == 'se':
            for start_idx_each_chunk,end_idx_each_chunk in bins:
                n_line_each_chunk = end_idx_each_chunk-start_idx_each_chunk
                denovo_se(n_line_each_chunk, start_idx_each_chunk)
                # print(f'{start_idx_each_chunk} , {n_line_each_chunk}, {end_idx_each_chunk}')
        else:
            print('sequence_type error')
        print(f'end={datetime.datetime.now()}')

        cmd = f'touch denovo/frag_denovo_done'
        subprocess.call(cmd, shell=True)

    elif module == 'create_cmd':
        # create command to split across nodes
        file_path = f'denovo/{NAME}2_sgrep.tsv'
        read_fraction_of_file(file_path, num_node)

    elif module =='parse_bc_fasta':
        n_bc = 1000
        in_fasta = 'test_1k_sgrep.tsv'
        parse_bc_fasta(in_fasta)
    else:
        print('module not found')


