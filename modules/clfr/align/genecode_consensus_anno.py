#!/usr/bin/env python
"""Composition of consensus.fasta from a minimap2 PAF against a transcript reference.

Best hit per query = max matching bases (col 10). Target name formats:
  GENCODE : ENST|ENSG|OTTHUMG|OTTHUMT|tx_name|gene|len|biotype|   -> gene + biotype from header
  RefSeq  : NM_/XM_ (mRNA), NR_/XR_ (ncRNA)                       -> coarse class from prefix;
            give --map (TSV: accession<TAB>gene<TAB>biotype) for gene names / real biotypes.
Queries with best-hit query coverage < --min_cov are 'low_cov'; queries absent from the PAF
are 'unmapped' (needs --fasta for the total).

usage: genecode_consensus_anno.py in.paf out_prefix [--fasta cons.fasta] [--map map.tsv] [--min_cov 0.8]
outputs: <prefix>.per_read.tsv  <prefix>.composition.tsv  <prefix>.genes.tsv  <prefix>.png
"""
import argparse, re, collections
import pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

GLOBIN = re.compile(r"^HB(A1|A2|B|D|E1|G1|G2|M|Q1|Z)$")
REFSEQ = {"NM_": "mRNA (NM_)", "XM_": "mRNA predicted (XM_)",
          "NR_": "ncRNA (NR_)", "XR_": "ncRNA predicted (XR_)"}

def group(biotype):
    if biotype == "protein_coding": return "protein_coding"
    if biotype == "lncRNA": return "lncRNA"
    if "rRNA" in biotype: return "rRNA"
    if "pseudogene" in biotype: return "pseudogene"
    if biotype.startswith(("IG_", "TR_")): return "IG/TR"
    if biotype in ("retained_intron", "nonsense_mediated_decay", "processed_transcript"):
        return "protein_coding-gene non-coding isoform"
    return biotype or "other"

def parse_target(t, amap):
    if "|" in t:
        a = t.split("|")
        return a[5], a[7]
    acc = t.split(".")[0]
    if acc in amap: return amap[acc]
    return "", REFSEQ.get(t[:3], "other")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("paf"); ap.add_argument("prefix")
    ap.add_argument("--fasta"); ap.add_argument("--map")
    ap.add_argument("--min_cov", type=float, default=0.8)
    ap.add_argument("--no_genes", action="store_true", help="do not plot top genes")
    a = ap.parse_args()

    amap = {}
    if a.map:
        for l in open(a.map):
            f = l.rstrip("\n").split("\t")
            if len(f) >= 3: amap[f[0].split(".")[0]] = (f[1], f[2])

    best = {}   # qname -> (matches, row)
    for l in open(a.paf):
        f = l.split("\t", 12)
        q, m = f[0], int(f[9])
        if q not in best or m > best[q][0]:
            best[q] = (m, f)
    rows = []
    for q, (m, f) in best.items():
        qlen, qs, qe = int(f[1]), int(f[2]), int(f[3])
        gene, bt = parse_target(f[5], amap)
        cov = (qe - qs) / qlen if qlen else 0
        cls = "globin" if GLOBIN.match(gene) else group(bt)
        if cls not in ("globin", "rRNA") and cov < a.min_cov: cls = "low_cov"   # globin/rRNA kept regardless of coverage
        rows.append((q, qlen, round(cov, 3), f[5], gene, bt, cls))
    d = pd.DataFrame(rows, columns=["id", "qlen", "qcov", "target", "gene", "biotype", "class"])

    total = len(d)
    if a.fasta:
        total = sum(1 for l in open(a.fasta) if l.startswith(">"))
    comp = d["class"].value_counts()
    if total > len(d): comp["unmapped"] = total - len(d)
    comp = comp.sort_values(ascending=False)
    out = pd.DataFrame({"n": comp, "pct": (100 * comp / total).round(2)})
    d.to_csv(a.prefix + ".per_read.tsv", sep="\t", index=False)
    out.to_csv(a.prefix + ".composition.tsv", sep="\t")
    genes = d[d.gene != ""].gene.value_counts()
    genes.to_csv(a.prefix + ".genes.tsv", sep="\t", header=["n"])
    print(f"total={total} with_hit={len(d)}\n{out.to_string()}")

    has_genes = len(genes) > 0 and not a.no_genes
    fig, ax = plt.subplots(1, 2 if has_genes else 1, figsize=(11 if has_genes else 6, 4.5), squeeze=False)
    s = out["n"].sort_values()
    bars = ax[0][0].barh(s.index, 100 * s / total, color="#4c78a8")
    for b, v in zip(bars, s):
        ax[0][0].text(b.get_width() + 0.5, b.get_y() + b.get_height() / 2,
                      f"{v} ({100*v/total:.1f}%)", va="center", fontsize=8)
    ax[0][0].set_title(f"Composition (n={total})"); ax[0][0].set_xlabel("% of consensus")
    ax[0][0].set_xlim(0, 120)
    if has_genes:
        g = genes.head(15)[::-1]
        ax[0][1].barh(g.index, g.values, color="#72b7b2")
        ax[0][1].set_title("Top genes"); ax[0][1].set_xlabel("consensus count")
    for x in ax.ravel(): x.spines[["top", "right"]].set_visible(False)
    plt.tight_layout(); plt.savefig(a.prefix + ".png", dpi=150)

if __name__ == "__main__":
    main()
