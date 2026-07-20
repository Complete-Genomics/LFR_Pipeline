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
import shutil
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
parser.add_argument("--start_idx", type=int, default=0, required=False)
parser.add_argument("--end_idx", type=int, required=False)
parser.add_argument("--module", type=str)
parser.add_argument("--num_node", type=int, required=False)
parser.add_argument("--nth_of_nodes", type=int, required=False)
parser.add_argument('--debug', action='store_true', default=False, help='Enable debug mode')
parser.add_argument("--megahit", type=str, default='megahit', help='Path to megahit binary')
parser.add_argument("--tmp_dir", type=str, default="/dev/shm", help="Root directory for megahit temporary outputs")
# parser.add_argument("--rg", type=str, default='rg', help='Path to rg (ripgrep) binary')


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
TMP_DIR = Path(args.tmp_dir)
# RG = args.rg

class NullLock:
    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        return False


def tmp_root():
    return TMP_DIR / f"{BATCH_LANE}_{ID}"


def barcode_tmp_dir(barcode):
    return tmp_root() / barcode


def cleanup_path(path):
    if DEBUG:
        return
    shutil.rmtree(path, ignore_errors=True)


def append_contigs_and_cleanup(barcode, returncode, stderr, lock):
    barcode_dir = barcode_tmp_dir(barcode)
    contigs = barcode_dir / f"{barcode}.contigs.fa"
    try:
        if returncode != 0:
            print(f"WARNING: megahit failed for barcode {barcode}: {stderr.decode(errors='replace').strip()}", flush=True)
            return
        if not contigs.exists():
            print(f"WARNING: megahit output missing for barcode {barcode}: {contigs}", flush=True)
            return

        bc = f'>{barcode}'
        with lock:
            with open(contigs) as source, open(f"denovo/final_contigs_{ID}.fa", "a") as out:
                for line_no, line in enumerate(source, 1):
                    if line_no in (1, 3, 5, 7):
                        out.write(bc + line)
                    else:
                        out.write(line)
    finally:
        cleanup_path(barcode_dir)


def process_barcode_se(barcode, shared_meta_data2, lock):
    K_MIN = 41
    K_MAX = 41
    MIN_CTG_LEN = 400
    megahit = MEGAHIT
    # rg = RG

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
            f"-o {barcode_tmp_dir(barcode)} --out-prefix {barcode} --k-min {K_MIN} --k-max {K_MAX} --force "
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

        append_contigs_and_cleanup(barcode, process.returncode, stderr, lock)

    except Exception as e: print(e)
    return


def process_barcode_pe(barcode, shared_meta_data1, shared_meta_data2, lock):
    K_MIN = 41
    K_MAX = 41
    MIN_CTG_LEN = 400
    megahit = MEGAHIT

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
            f"-o {barcode_tmp_dir(barcode)} --out-prefix {barcode} --k-min {K_MIN} --k-max {K_MAX} --force "
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

        append_contigs_and_cleanup(barcode, process.returncode, stderr, lock)

    except Exception as e: print(e)
    return

def add_sgrep_line(meta_data, line):
    bc_len = 15
    info = line.strip().split('\t')
    if len(info) < 2:
        return False
    bc = info[0][5:5+bc_len]
    id = '>'+info[0][22:]
    seq = info[1]
    meta_data[bc].append(id)
    meta_data[bc].append(seq)
    return True


def process_pe_metadata(meta_data1, meta_data2, start_idx):

    # subprocess.call(f'mkdir -p /dev/shm/{BATCH_LANE}_{ID}', shell=True)

    if num_processes == 1:
        lock = NullLock()
        for barcode in meta_data2.keys():
            process_barcode_pe(barcode, meta_data1, meta_data2, lock)
    else:
        with mp.Manager() as manager:
            shared_meta_data1 = manager.dict(meta_data1)
            shared_meta_data2 = manager.dict(meta_data2)
            lock = manager.Lock()

            with mp.Pool(num_processes) as pool:
                pool.starmap(process_barcode_pe, [(barcode, shared_meta_data1, shared_meta_data2, lock) for barcode in meta_data2.keys()])

    print(f'start_idx={start_idx}')
    print(f'denovo_BC_counts={len(meta_data2)}')
    return sum(len(v) // 2 for v in meta_data2.values())


def process_se_metadata(meta_data2, start_idx):
    if num_processes == 1:
        lock = NullLock()
        for barcode in meta_data2.keys():
            process_barcode_se(barcode, meta_data2, lock)
    else:
        with mp.Manager() as manager:
            # shared_meta_data1 = manager.dict(meta_data1)
            shared_meta_data2 = manager.dict(meta_data2)
            lock = manager.Lock()

            with mp.Pool(num_processes) as pool:
                pool.starmap(process_barcode_se, [(barcode, shared_meta_data2, lock) for barcode in meta_data2.keys()])

    print(f'start_idx={start_idx}')
    print(f'denovo_BC_counts={len(meta_data2)}')
    return sum(len(v) // 2 for v in meta_data2.values())


def denovo_pe(n_line_chunk, start_idx):
    meta_data1 = defaultdict(list)
    meta_data2 = defaultdict(list)

    ## load reads to mem to speedup, *sgrep.tsv with diff len after trim, unable to use seek to get idx
    with open(f'denovo/{NAME}1_sgrep.tsv', 'r') as f:
        for line in itertools.islice(f, start_idx, start_idx+n_line_chunk):
            add_sgrep_line(meta_data1, line)

    with open(f'denovo/{NAME}2_sgrep.tsv', 'r') as f:
        for line in itertools.islice(f, start_idx, start_idx+n_line_chunk):
            add_sgrep_line(meta_data2, line)

    return process_pe_metadata(meta_data1, meta_data2, start_idx)


def denovo_se(n_line_chunk, start_idx):
    meta_data2 = defaultdict(list)

    ## load reads to mem to speedup, *sgrep.tsv with diff len after trim, unable to use seek to get idx
    with open(f'denovo/{NAME}2_sgrep.tsv', 'r') as f:
        for line in itertools.islice(f, start_idx, start_idx+n_line_chunk):
            add_sgrep_line(meta_data2, line)

    return process_se_metadata(meta_data2, start_idx)


def iter_denovo_se_chunks(start_idx, n_line_chunk):
    with open(f'denovo/{NAME}2_sgrep.tsv', 'r') as f:
        for _ in itertools.islice(f, start_idx):
            pass
        chunk_start = start_idx
        while True:
            meta_data2 = defaultdict(list)
            n_lines = 0
            for line in itertools.islice(f, n_line_chunk):
                if add_sgrep_line(meta_data2, line):
                    n_lines += 1
            if n_lines == 0:
                break
            yield chunk_start, meta_data2
            chunk_start += n_lines


def iter_denovo_pe_chunks(start_idx, n_line_chunk):
    with open(f'denovo/{NAME}1_sgrep.tsv', 'r') as f1, open(f'denovo/{NAME}2_sgrep.tsv', 'r') as f2:
        for _ in itertools.islice(f1, start_idx):
            pass
        for _ in itertools.islice(f2, start_idx):
            pass
        chunk_start = start_idx
        while True:
            meta_data1 = defaultdict(list)
            meta_data2 = defaultdict(list)
            n_lines = 0
            for line1, line2 in itertools.islice(zip(f1, f2), n_line_chunk):
                ok1 = add_sgrep_line(meta_data1, line1)
                ok2 = add_sgrep_line(meta_data2, line2)
                if ok1 and ok2:
                    n_lines += 1
            if n_lines == 0:
                break
            yield chunk_start, meta_data1, meta_data2
            chunk_start += n_lines

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

def count_sgrep_lines():
    file_path = f'denovo/{NAME}2_sgrep.tsv'
    with open(file_path, 'r') as f:
        return sum(1 for _ in f)

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

        tmp_root().mkdir(parents=True, exist_ok=True)
        try:
            if end_idx is None:
                print(
                    f'end_idx not specified; streaming all reads from start_idx={start_idx}',
                    flush=True,
                )
                if sequence_type == 'pe':
                    for chunk_start, meta_data1, meta_data2 in iter_denovo_pe_chunks(start_idx, n_line_chunk):
                        print(f'processing chunk start_idx={chunk_start} reads={sum(len(v) // 2 for v in meta_data2.values())}', flush=True)
                        process_pe_metadata(meta_data1, meta_data2, chunk_start)
                elif sequence_type == 'se':
                    for chunk_start, meta_data2 in iter_denovo_se_chunks(start_idx, n_line_chunk):
                        print(f'processing chunk start_idx={chunk_start} reads={sum(len(v) // 2 for v in meta_data2.values())}', flush=True)
                        process_se_metadata(meta_data2, chunk_start)
                else:
                    raise SystemExit(f'sequence_type error: {sequence_type}')
            else:
                if end_idx <= start_idx:
                    raise SystemExit(f"Invalid denovo range: start_idx={start_idx}, end_idx={end_idx}")
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
                    raise SystemExit(f'sequence_type error: {sequence_type}')
        finally:
            cleanup_path(tmp_root())
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
