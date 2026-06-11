"""
cwgs_preprocess.py — barcode split for cWGS data

Usage:
    python cwgs_preprocess.py config.yaml

Output:
    split_read.1.fq.gz
    split_read.2.fq.gz
    split_stat_read1.log
    split_stat_read.err
input:
    read_1.fq.gz
    read_2.fq.gz
"""

import sys
import gzip
import subprocess
import yaml


# ── config loading ────────────────────────────────────────────────────────────

def load_config(path):
    with open(path) as f:
        return yaml.safe_load(f)


# ── barcode list detection (ported from splitreads.smk) ───────────────────────

def determine_barcode_list(config, n_samples=400000):
    """Auto-detect whether barcode.list or barcode_RC.list matches R2 reads.
    Returns the absolute path to the correct barcode list file.
    """
    toolsdir = config['params']['toolsdir']
    # bc_condition = config['params']['bc_condition'].lower()

    barcode_file    = f"{toolsdir}/barcode.list"
    barcode_rc_file = f"{toolsdir}/barcode_RC.list"

    fastq_path = "read_2.fq.gz"

    def get_barcodes(barcodes_file):
        nucs = ['A', 'C', 'G', 'T']
        barcodes = []
        with open(barcodes_file) as f:
            for line in f:
                bc = line.strip().split()[0]
                for i in range(len(bc)):
                    for nuc in nucs:
                        bc_alt = list(bc)
                        bc_alt[i] = nuc
                        barcodes.append("".join(bc_alt))
        return set(barcodes)

    def get_fq_barcodes(fastq_path, n_samples):
        bc_start_idx = config['params']['bc_start'] - 1
        _bc_len, _bc_gap = 10, 6
        fq_bcs = []
        counter = 1
        print("Sampling R2 barcodes...", file=sys.stderr)
        with gzip.open(fastq_path, "rb") as fq:
            for line in fq:
                counter += 1
                if counter == n_samples + 1:
                    break
                if counter % 4 == 3:
                    seq = line.decode('utf-8').strip()
                    b1 = seq[bc_start_idx                              : bc_start_idx + _bc_len]
                    b2 = seq[bc_start_idx + _bc_len + _bc_gap          : bc_start_idx + _bc_gap + _bc_len * 2]
                    b3 = seq[bc_start_idx + _bc_gap * 2 + _bc_len * 2  : bc_start_idx + _bc_gap * 2 + _bc_len * 3]
                    fq_bcs += [b1, b2, b3]
        return fq_bcs

    barcodes    = get_barcodes(barcode_file)
    barcodes_rc = get_barcodes(barcode_rc_file)
    fq_bcs      = get_fq_barcodes(fastq_path, n_samples)

    bcs_found    = sum(1 for bc in fq_bcs if bc in barcodes)
    rc_bcs_found = sum(1 for bc in fq_bcs if bc in barcodes_rc)
    print(f"Barcodes matched: {bcs_found}  RC matched: {rc_bcs_found}", file=sys.stderr)

    if bcs_found > rc_bcs_found:
        print(f"Using {barcode_file}", file=sys.stderr)
        return barcode_file
    elif rc_bcs_found > bcs_found:
        print(f"Using {barcode_rc_file}", file=sys.stderr)
        return barcode_rc_file
    else:
        if n_samples > 4_000_000:
            print("ERROR: cannot determine barcode list", file=sys.stderr)
            sys.exit(1)
        return determine_barcode_list(config, n_samples * 2)


# ── main ──────────────────────────────────────────────────────────────────────

def main():
    if len(sys.argv) != 2:
        print("Usage: python cwgs_preprocess.py config.yaml")
        sys.exit(1)

    config = load_config(sys.argv[1])
    p = config['params']

    src_dir   = p['src_dir'].rstrip("/")
    python    = p.get('general_python') or sys.executable
    r2_len    = p['read_len']
    r1_len    = p['read_len_r1']
    bc_start  = p['bc_start']
    gdna_start          = p['gdna_start']
    additional_bc_start = p['additional_bc_start']
    additional_bc_len   = p['additional_bc_len']
    gdna_start_r1       = p['gdna_start_r1']

    print("Determining barcode list...", file=sys.stderr)
    barcode = determine_barcode_list(config)

    threads = config.get('threads', {}).get('split_group', 1)
    cmd = [
        python,
        f"{src_dir}/modules/shared/splitreads/split_barcode_stLFR.py",
        "--barcode", barcode,
        "--r1", "read_1.fq.gz",
        "--r2", "read_2.fq.gz",
        "--read_len", str(r2_len),
        "--output", "split_read",
        "--bc_start", str(bc_start),
        "--gdna_start", str(gdna_start),
        "--additional_bc_start", str(additional_bc_start),
        "--additional_bc_len", str(additional_bc_len),
        "--gdna_start_r1", str(gdna_start_r1),
        "--read_len_r1", str(r1_len),
        "--swap", "none",
        "--output_mode", "cwgs",
        "--threads", str(threads),
    ]
    print("Running:\n  %s" % " ".join(cmd), file=sys.stderr)
    with open("split_stat_read.err", "w") as err:
        ret = subprocess.call(cmd, stderr=err)
    if ret != 0:
        print(f"ERROR: split_barcode_stLFR.py exited with code {ret}", file=sys.stderr)
        sys.exit(ret)
    print("Done.", file=sys.stderr)


if __name__ == "__main__":
    main()
