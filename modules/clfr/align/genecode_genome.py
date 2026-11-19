#!/usr/bin/env python3
"""Where do the sampled reads land in the GENOME? Cross it with the transcript class from genecode_classify.py.

The competitive reference (transcripts + rRNA + genome decoy, linear map-ont) cannot say what an 'unmapped' / 'low_cov' /
'genome_decoy' read is. Here the same reads are mapped spliced to the genome (minimap2 -x splice) and classified by the
GENCODE gtf: exon_same_strand / exon_antisense (overlaps an annotated exon), intronic (inside a gene, no exon overlap),
intergenic, genome_low_cov (best genome hit covers < --min_cov of the read), genome_no_hit. Overlap is by interval only.

usage: genecode_genome.py reads.paf gtf.gz reads_class.tsv out_prefix [--min_cov 0.5] [--fq reads.fq.gz]
  --fq = the mapped reads: adds n_bases / pct_bases_of_transcript_class to _by_class.tsv (reads are very uneven in length,
  so a read count over-weights the short, low-quality ones)
  reads_class.tsv = genecode_classify.py --out_reads (defines the reads and the denominator)
outputs: <out_prefix>.tsv (genome class: n, pct)   <out_prefix>_by_class.tsv (transcript class x genome class, mean exon_frac)
         <out_prefix>_reads.tsv (per read: genome class, location, exon_frac = share of the alignment covered by same-strand
         exons, 0 = entirely intron)
"""
import argparse, collections, gzip, os, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from genecode_consensus_genome import exon_frac, load_gtf   # same interval index / gtf parser as the consensus version


def read_lengths(fq):
    """read name (first word of the header) -> length."""
    out = {}
    with (gzip.open(fq, "rt") if fq.endswith(".gz") else open(fq)) as fh:
        for i, l in enumerate(fh):
            if i % 4 == 0: name = l[1:].split()[0]
            elif i % 4 == 1: out[name] = len(l.strip())
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("paf"); ap.add_argument("gtf"); ap.add_argument("reads_class"); ap.add_argument("prefix")
    ap.add_argument("--min_cov", type=float, default=0.5)
    ap.add_argument("--fq", help="the mapped reads; adds base counts to _by_class.tsv")
    a = ap.parse_args()

    best = {}   # read -> (aligned read bases, PAF fields)
    for l in open(a.paf):
        f = l.split("\t", 12)
        span = int(f[3]) - int(f[2])
        if f[0] not in best or span > best[f[0]][0]: best[f[0]] = (span, f)
    genes, exons = load_gtf(a.gtf)
    lens = read_lengths(a.fq) if a.fq else {}
    nbase = collections.Counter()

    cross, total, fsum = collections.Counter(), collections.Counter(), collections.Counter()
    per_read = open(a.prefix + "_reads.tsv", "w")
    per_read.write("read\ttranscript_class\tgenome_class\tchr\tstart\tend\tstrand\tqcov\texon_frac\n")
    with open(a.reads_class) as fh:
        next(fh)
        for l in fh:
            r, tcls = l.rstrip("\n").split("\t")
            ef, loc, qc = 0.0, ("", "", "", ""), ""
            if r not in best: g = "genome_no_hit"
            else:
                span, f = best[r]
                strand, chrom, ts, te = f[4], f[5], int(f[7]), int(f[8])
                x = exons.overlap(chrom, ts, te)
                ef, loc, qc = exon_frac(x, strand, ts, te), (chrom, ts, te, strand), f"{span / int(f[1]):.3f}"
                if span / int(f[1]) < a.min_cov: g = "genome_low_cov"
                elif x: g = "exon_same_strand" if any(e[2] == strand for e in x) else "exon_antisense"
                elif genes.overlap(chrom, ts, te): g = "intronic"
                else: g = "intergenic"
            cross[(tcls, g)] += 1; total[g] += 1; fsum[(tcls, g)] += ef; nbase[(tcls, g)] += lens.get(r, 0)
            per_read.write("\t".join(map(str, (r, tcls, g, *loc, qc, f"{ef:.3f}"))) + "\n")
    per_read.close()
    n = sum(total.values())
    with open(a.prefix + ".tsv", "w") as o:
        o.write("genome_class\tn_reads\tpct_reads\n")
        for g, v in total.most_common(): o.write(f"{g}\t{v}\t{100*v/n:.3f}\n")
        o.write(f"total\t{n}\t100.000\n")
    tn, tb = collections.Counter(), collections.Counter()
    for (t, g), v in cross.items(): tn[t] += v; tb[t] += nbase[(t, g)]
    with open(a.prefix + "_by_class.tsv", "w") as o:
        o.write("transcript_class\tgenome_class\tn_reads\tpct_of_transcript_class\tmean_exon_frac"
                + ("\tn_bases\tpct_bases_of_transcript_class" if a.fq else "") + "\n")
        for (t, g), v in sorted(cross.items(), key=lambda kv: (-tn[kv[0][0]], -kv[1])):
            o.write(f"{t}\t{g}\t{v}\t{100*v/tn[t]:.3f}\t{fsum[(t, g)]/v:.3f}")
            if a.fq: o.write(f"\t{nbase[(t, g)]}\t{100*nbase[(t, g)]/max(1, tb[t]):.3f}")
            o.write("\n")
    for g, v in total.most_common(): print(f"{g}\t{v}\t{100*v/n:.1f}%", file=sys.stderr)


if __name__ == "__main__":
    main()
