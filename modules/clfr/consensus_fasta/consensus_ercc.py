#!/usr/bin/env python3
"""Summarize ERCC consensus FASTA length recovery and variants from PAF.

Use an existing minimap2 PAF when available. For exact SNP/indel positions,
generate PAF with minimap2 --cs=long. Without cs tags, the script still reports
alignment coverage and NM mismatch/edit distance when present, but cannot list
per-base variants.
"""

from __future__ import print_function

import argparse
import csv
import os
import re
import sys


ERCC_RE = re.compile(r"(ERCC-\d+)")
CS_RE = re.compile(r"(:[0-9]+|\*[a-z][a-z]|\+[a-z]+|-[a-z]+|=[a-z]+|~[a-z]{2}[0-9]+[a-z]{2})")


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--fasta", required=True, help="Consensus FASTA file.")
    parser.add_argument("--paf", default="consensus/consensus.paf", help="minimap2 PAF for the consensus FASTA.")
    parser.add_argument("--ercc_ref", default=None, help="ERCC truth table with 'ERCC ID' and 'Sequence'.")
    parser.add_argument("--summary", default="consensus/consensus_ercc_summary.tsv")
    parser.add_argument("--variants", default="consensus/consensus_ercc_variants.tsv")
    return parser.parse_args()


def default_ercc_ref():
    here = os.path.dirname(os.path.abspath(__file__))
    return os.path.abspath(os.path.join(here, "..", "calc_frag_len", "ercc_truth.txt"))


def read_fasta_lengths(path):
    lengths = {}
    name = None
    length = 0
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    lengths[name] = length
                name = line[1:].split()[0]
                length = 0
            else:
                length += len(line)
    if name is not None:
        lengths[name] = length
    return lengths


def read_ercc_truth(path):
    refs = {}
    with open(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if "ERCC ID" not in reader.fieldnames or "Sequence" not in reader.fieldnames:
            raise ValueError("ERCC reference must contain 'ERCC ID' and 'Sequence'")
        for row in reader:
            ercc_id = row["ERCC ID"].strip()
            seq = row["Sequence"].strip().upper()
            if ercc_id and seq:
                refs[ercc_id] = seq
    return refs


def ercc_from_text(text):
    match = ERCC_RE.search(text)
    return match.group(1) if match else None


def parse_tags(fields):
    tags = {}
    for field in fields[12:]:
        parts = field.split(":", 2)
        if len(parts) == 3:
            tags[parts[0]] = parts[2]
    return tags


def read_best_ercc_paf(path, ercc_refs):
    best = {}
    with open(path) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                continue
            query = fields[0]
            target = fields[5]
            ercc_id = ercc_from_text(target) or ercc_from_text(query)
            if ercc_id not in ercc_refs:
                continue
            record = {
                "query": query,
                "query_len": int(fields[1]),
                "query_start": int(fields[2]),
                "query_end": int(fields[3]),
                "strand": fields[4],
                "target": target,
                "target_len": int(fields[6]),
                "target_start": int(fields[7]),
                "target_end": int(fields[8]),
                "matches": int(fields[9]),
                "block_len": int(fields[10]),
                "mapq": int(fields[11]),
                "tags": parse_tags(fields),
                "ercc_id": ercc_id,
            }
            old = best.get(query)
            if old is None or (record["matches"], record["block_len"], record["mapq"]) > (
                old["matches"], old["block_len"], old["mapq"]
            ):
                best[query] = record
    return best


def count_variants_from_cs(query, ercc_id, cs, target_start, strand):
    variants = []
    snp_count = 0
    ins_count = 0
    del_count = 0
    ref_pos = target_start + 1

    for token in CS_RE.findall(cs):
        op = token[0]
        payload = token[1:]
        if op == ":":
            ref_pos += int(payload)
        elif op == "=":
            ref_pos += len(payload)
        elif op == "*":
            ref_base = payload[0].upper()
            query_base = payload[1].upper()
            snp_count += 1
            variants.append({
                "consensus_id": query,
                "ercc_id": ercc_id,
                "variant_type": "snp",
                "ref_pos_1based": ref_pos,
                "ref_base": ref_base,
                "consensus_base": query_base,
                "strand": strand,
            })
            ref_pos += 1
        elif op == "+":
            ins_count += len(payload)
            variants.append({
                "consensus_id": query,
                "ercc_id": ercc_id,
                "variant_type": "insertion",
                "ref_pos_1based": ref_pos,
                "ref_base": "-",
                "consensus_base": payload.upper(),
                "strand": strand,
            })
        elif op == "-":
            del_count += len(payload)
            variants.append({
                "consensus_id": query,
                "ercc_id": ercc_id,
                "variant_type": "deletion",
                "ref_pos_1based": ref_pos,
                "ref_base": payload.upper(),
                "consensus_base": "-",
                "strand": strand,
            })
            ref_pos += len(payload)
        elif op == "~":
            match = re.match(r"[a-z]{2}([0-9]+)[a-z]{2}", payload)
            if match:
                ref_pos += int(match.group(1))

    return snp_count, ins_count, del_count, variants


def write_table(path, rows, fieldnames):
    directory = os.path.dirname(path)
    if directory:
        os.makedirs(directory, exist_ok=True)
    with open(path, "w") as out:
        writer = csv.DictWriter(out, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def main():
    args = parse_args()
    ercc_ref = args.ercc_ref or default_ercc_ref()
    ercc_refs = read_ercc_truth(ercc_ref)
    fasta_lengths = read_fasta_lengths(args.fasta)
    best_paf = read_best_ercc_paf(args.paf, ercc_refs)

    summaries = []
    variants = []
    missing_cs = 0
    for query, record in sorted(best_paf.items()):
        ercc_id = record["ercc_id"]
        ref_len = len(ercc_refs[ercc_id])
        consensus_len = fasta_lengths.get(query, record["query_len"])
        aligned_ref_bases = record["target_end"] - record["target_start"]
        aligned_query_bases = record["query_end"] - record["query_start"]
        nm = record["tags"].get("NM", "")
        cs = record["tags"].get("cs")

        snp_count = ""
        ins_count = ""
        del_count = ""
        if cs:
            snp_count, ins_count, del_count, record_variants = count_variants_from_cs(
                query, ercc_id, cs, record["target_start"], record["strand"]
            )
            variants.extend(record_variants)
        else:
            missing_cs += 1

        summaries.append({
            "consensus_id": query,
            "ercc_id": ercc_id,
            "ref_len": ref_len,
            "consensus_len": consensus_len,
            "consensus_len_pct_ref": 100.0 * consensus_len / ref_len if ref_len else 0,
            "aligned_ref_bases": aligned_ref_bases,
            "aligned_ref_pct": 100.0 * aligned_ref_bases / ref_len if ref_len else 0,
            "aligned_query_bases": aligned_query_bases,
            "query_aligned_pct": 100.0 * aligned_query_bases / consensus_len if consensus_len else 0,
            "matches": record["matches"],
            "block_len": record["block_len"],
            "mapq": record["mapq"],
            "strand": record["strand"],
            "target_start_1based": record["target_start"] + 1,
            "target_end_1based": record["target_end"],
            "nm": nm,
            "snp_count": snp_count,
            "insertion_count": ins_count,
            "deletion_count": del_count,
            "has_snp": "yes" if snp_count != "" and snp_count > 0 else ("unknown" if snp_count == "" else "no"),
            "has_indel": "yes" if ins_count != "" and (ins_count + del_count) > 0 else ("unknown" if ins_count == "" else "no"),
        })

    summary_fields = [
        "consensus_id", "ercc_id", "ref_len", "consensus_len", "consensus_len_pct_ref",
        "aligned_ref_bases", "aligned_ref_pct", "aligned_query_bases", "query_aligned_pct",
        "matches", "block_len", "mapq", "strand", "target_start_1based", "target_end_1based",
        "nm", "snp_count", "insertion_count", "deletion_count", "has_snp", "has_indel",
    ]
    variant_fields = [
        "consensus_id", "ercc_id", "variant_type", "ref_pos_1based",
        "ref_base", "consensus_base", "strand",
    ]
    write_table(args.summary, summaries, summary_fields)
    write_table(args.variants, variants, variant_fields)

    if missing_cs:
        sys.stderr.write(
            "WARNING: %s PAF records do not have cs tags; exact SNP/indel positions are unknown. "
            "Regenerate PAF with minimap2 --cs=long to enable variant calls.\n" % missing_cs
        )
    sys.stderr.write("Wrote %s ERCC consensus summaries and %s variants.\n" % (
        len(summaries), len(variants)
    ))


if __name__ == "__main__":
    main()
