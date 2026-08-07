#!/usr/bin/env python3
"""Combine read-back self-consistency + reference-free chimera detection into
one per-contig QC flag, and split denovo.longest.fasta accordingly.

No reference database required (production soil/env samples have none) --
both signals are computed from the assembly's own contigs/reads:
  - read-back: does the contig's own UMI reads cover and bracket it?
    (readback_qc_single.py output)
  - chimera: vsearch --uchime_denovo on the contig set itself (reference-free;
    see denovo.md sec 26/29 for why ref-based detection isn't available here
    and why reference-free alone is weaker but still a real signal)

Flag rule (denovo.md sec 26 conclusion 3): low_confidence if
  breadth_1x < 0.95  OR  two_end_supported == 0  OR  chimera == 'Y'
"""
import argparse
import csv


def load_readback(path):
    rows = {}
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            rows[row["barcode"]] = row
    return rows


def load_junction_flags(path):
    """barcode -> (junction_suspect, span_cov_ratio) from denovo_junction_qc.py."""
    flags = {}
    if not path:
        return flags
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            flags[row["barcode"]] = (int(row["junction_suspect"]),
                                      float(row["span_cov_ratio"]))
    return flags


def load_chimera_flags(uchimeout_path):
    flags = {}
    if not uchimeout_path:
        return flags
    with open(uchimeout_path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 18:
                continue
            barcode = f[1].split(">")[0]
            flags[barcode] = f[17]
    return flags


def load_fasta(path):
    seqs = {}
    name = None
    chunks = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(chunks)
                name = line[1:].split(">", 1)[0]
                chunks = []
            else:
                chunks.append(line)
        if name is not None:
            seqs[name] = "".join(chunks)
    return seqs


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--contigs", required=True, help="denovo.longest.fasta")
    ap.add_argument("--readback", required=True, help="readback_qc_single.py output tsv")
    ap.add_argument("--uchimeout", required=True, help="vsearch --uchime_denovo --uchimeout tsv")
    ap.add_argument("--junction", help="denovo_junction_qc.py output tsv")
    ap.add_argument("--out-report", required=True)
    ap.add_argument("--out-highconf-fasta", required=True)
    ap.add_argument("--out-flagged-fasta", required=True)
    args = ap.parse_args()

    contigs = load_fasta(args.contigs)
    readback = load_readback(args.readback)
    chimera = load_chimera_flags(args.uchimeout)
    junction = load_junction_flags(args.junction)

    fields = ["barcode", "contig_len", "breadth_1x", "two_end_supported",
              "chimera", "span_cov_ratio", "junction_suspect", "low_confidence"]
    n_low = 0
    with open(args.out_report, "w", newline="") as out_report, \
         open(args.out_highconf_fasta, "w") as out_hi, \
         open(args.out_flagged_fasta, "w") as out_lo:
        writer = csv.DictWriter(out_report, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for barcode, seq in contigs.items():
            rb = readback.get(barcode)
            breadth_1x = float(rb["breadth_1x"]) if rb else None
            two_end = int(rb["two_end_supported"]) if rb else None
            chim = chimera.get(barcode, "N")
            jflag, jratio = junction.get(barcode, (0, None))
            # junction_suspect is the only one of these with validated chimera
            # discrimination (denovo.md sec 31); the others cover read-support
            # failures, which is a different defect.
            # two_end_supported is reported but deliberately NOT part of the
            # gate: on the Zymo control it removed 34% of contigs while moving
            # mean identity only 94.31 -> 94.15 (denovo.md sec 30), i.e. it
            # costs a third of the yield for no measurable chimera benefit.
            low_confidence = (
                rb is None
                or breadth_1x < 0.95
                or chim == "Y"
                or jflag == 1
            )
            writer.writerow({
                "barcode": barcode,
                "contig_len": len(seq),
                "breadth_1x": breadth_1x,
                "two_end_supported": two_end,
                "chimera": chim,
                "span_cov_ratio": jratio,
                "junction_suspect": jflag,
                "low_confidence": int(low_confidence),
            })
            target = out_lo if low_confidence else out_hi
            target.write(f">{barcode}\n{seq}\n")
            n_low += int(low_confidence)

    n = len(contigs)
    print(f"contigs_total={n}")
    print(f"low_confidence={n_low} ({100*n_low/n:.2f}%)" if n else "low_confidence=0")
    print(f"high_confidence={n - n_low} ({100*(n - n_low)/n:.2f}%)" if n else "high_confidence=0")


if __name__ == "__main__":
    main()