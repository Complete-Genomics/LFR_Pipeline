#!/usr/bin/env python3
import argparse
import gzip
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--dirname", type=str, default=".")
    parser.add_argument("--r1r2", type=int, choices=(1, 2), required=True)
    parser.add_argument("--bc_list", type=str, required=True)
    parser.add_argument("--outdir", type=str, default="fq")
    return parser.parse_args()


def load_barcodes(path):
    with open(path) as handle:
        return {line.strip() for line in handle if line.strip()}


def barcode_from_header(header):
    fields = header.strip().split()
    for field in fields[1:]:
        if field.startswith("BX:Z:"):
            return field[5:]

    read_id = fields[0]
    if "#" not in read_id:
        return None
    barcode = read_id.split("#", 1)[1]
    barcode = barcode.split("#", 1)[0]
    barcode = barcode.split("/", 1)[0]
    return barcode or None


def read_fastq_record(handle):
    header = handle.readline()
    if not header:
        return None
    seq = handle.readline()
    plus = handle.readline()
    qual = handle.readline()
    if not qual:
        raise ValueError("Truncated FASTQ record")
    return header, seq, plus, qual


def split_fastq_by_barcode(dirname, r1r2, bc_list, outdir):
    fastq = Path(dirname) / "data" / f"split_read.{r1r2}.fq.gz"
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    handles = {}
    counts = {bc: 0 for bc in bc_list}
    try:
        with gzip.open(fastq, "rt") as handle:
            while True:
                record = read_fastq_record(handle)
                if record is None:
                    break
                barcode = barcode_from_header(record[0])
                if barcode not in bc_list:
                    continue
                if barcode not in handles:
                    handles[barcode] = gzip.open(outdir / f"{barcode}.{r1r2}.fq.gz", "wt")
                handles[barcode].writelines(record)
                counts[barcode] += 1
    finally:
        for handle in handles.values():
            handle.close()

    for barcode in bc_list:
        path = outdir / f"{barcode}.{r1r2}.fq.gz"
        if not path.exists():
            with gzip.open(path, "wt"):
                pass

    with open(outdir / f"fq{r1r2}_done", "w") as out:
        for barcode in sorted(counts):
            out.write(f"{barcode}\t{counts[barcode]}\n")


def main():
    args = parse_args()
    split_fastq_by_barcode(
        dirname=args.dirname,
        r1r2=args.r1r2,
        bc_list=load_barcodes(args.bc_list),
        outdir=args.outdir,
    )


if __name__ == "__main__":
    main()
