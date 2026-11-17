#!/usr/bin/env python
"""What are the consensus that did not map to the GENCODE transcripts? Map them to the genome and classify.

The transcript reference has no intron / intergenic sequence, so 'unmapped' there is a mix of intronic (pre-mRNA /
retained intron), exon-boundary, antisense and intergenic molecules. This splits it up by where it lands in the genome.

  extract : consensus fasta + <prefix>.per_read.tsv (genecode_consensus_anno.py) -> fasta of the unmapped consensus
  classify: minimap2 -x splice PAF of that fasta vs the genome + GENCODE gtf -> per-consensus class + summary

classes: exon_same_strand / exon_antisense (alignment overlaps an annotated exon; same = any overlapped exon is on the
alignment strand), intronic (inside a gene, no exon overlap), intergenic, genome_low_cov (best genome hit covers
< --min_cov of the consensus), genome_no_hit. Overlap is by interval only; a spliced alignment is not checked exon by exon.

usage: genecode_consensus_genome.py extract in.fasta per_read.tsv out.fasta [--n N --seed S]
       genecode_consensus_genome.py classify in.paf gtf.gz out_prefix mapped.fasta --total_consensus T [--n_unmapped U --min_cov 0.8]
exon_frac = share of the alignment covered by exons on the alignment strand (0 = entirely intron; 'exon_same_strand'
only means ANY overlap, often a few bp at an exon edge, so read exon_frac together with the class).
outputs: <out_prefix>.tsv (per consensus) and <out_prefix>.summary.tsv (n, pct of unmapped, pct of all consensus)
"""
import argparse, collections, gzip, random, re, sys

BIN = 100_000


def read_fasta(path):
    name, seq = None, []
    for l in open(path):
        if l[0] == ">":
            if name is not None: yield name, "".join(seq)
            name, seq = l[1:].split()[0], []
        else:
            seq.append(l.strip())
    if name is not None: yield name, "".join(seq)


def extract(a):
    mapped = set()
    with open(a.per_read) as fh:
        next(fh)
        for l in fh: mapped.add(l.split("\t", 1)[0])
    un = [(n, s) for n, s in read_fasta(a.fasta) if n not in mapped]
    print(f"unmapped consensus: {len(un)}", file=sys.stderr)
    if a.n and a.n < len(un):
        random.seed(a.seed)
        un = random.sample(un, a.n)
        print(f"sampled {len(un)} (seed {a.seed})", file=sys.stderr)
    with open(a.out, "w") as f:
        for n, s in un: f.write(f">{n}\n{s}\n")


class Intervals:
    """chr -> 100 kb bins -> (start, end, strand, gene, gene_type); overlap query without scanning the whole chromosome."""
    def __init__(self): self.b = collections.defaultdict(lambda: collections.defaultdict(list))

    def add(self, c, s, e, rec):
        for i in range(s // BIN, (e - 1) // BIN + 1): self.b[c][i].append(rec)

    def overlap(self, c, s, e):
        seen, out = set(), []
        for i in range(s // BIN, (e - 1) // BIN + 1):
            for r in self.b[c].get(i, ()):
                if r[0] < e and r[1] > s and id(r) not in seen:
                    seen.add(id(r)); out.append(r)
        return out


def exon_frac(x, strand, s, e):
    """Fraction of the alignment [s, e) covered by the union of the overlapping exons on the alignment strand (0 = all intron)."""
    iv = sorted((max(r[0], s), min(r[1], e)) for r in x if r[2] == strand)
    if not iv or e <= s: return 0.0
    cov, (cs, ce) = 0, iv[0]
    for a, b in iv[1:]:
        if a > ce: cov += ce - cs; cs, ce = a, b
        else: ce = max(ce, b)
    return (cov + ce - cs) / (e - s)


def load_gtf(path):
    genes, exons = Intervals(), Intervals()
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as fh:
        for l in fh:
            if l[0] == "#": continue
            f = l.split("\t")
            if f[2] not in ("gene", "exon"): continue
            gn = re.search('gene_name "([^"]+)"', f[8]); bt = re.search('gene_type "([^"]+)"', f[8])
            rec = (int(f[3]) - 1, int(f[4]), f[6], gn.group(1) if gn else "", bt.group(1) if bt else "")
            (genes if f[2] == "gene" else exons).add(f[0], rec[0], rec[1], rec)
    return genes, exons


def classify(a):
    best = {}   # query -> (aligned query bases, PAF fields); one best genome hit per consensus
    for l in open(a.paf):
        f = l.rstrip("\n").split("\t")
        span = int(f[3]) - int(f[2])
        if f[0] not in best or span > best[f[0]][0]: best[f[0]] = (span, f)
    genes, exons = load_gtf(a.gtf)

    names = [l[1:].split()[0] for l in open(a.fasta) if l[0] == ">"]
    rows, cnt = [], collections.Counter()
    for q in names:
        if q not in best:
            cnt["genome_no_hit"] += 1; rows.append((q, "", "*", "", "", "", "", "", "genome_no_hit", "", "", "")); continue
        span, f = best[q]
        qlen, strand, chrom, ts, te = int(f[1]), f[4], f[5], int(f[7]), int(f[8])
        cov = span / qlen
        x, g = exons.overlap(chrom, ts, te), genes.overlap(chrom, ts, te)
        if cov < a.min_cov: k = "genome_low_cov"
        elif x: k = "exon_same_strand" if any(e[2] == strand for e in x) else "exon_antisense"
        elif g: k = "intronic"
        else: k = "intergenic"
        cnt[k] += 1
        gt = collections.Counter(e[4] for e in (x or g)).most_common(1)
        rows.append((q, qlen, chrom, ts, te, strand, f"{cov:.3f}", f[11], k,
                     ",".join(sorted({e[3] for e in (x or g)})[:3]), gt[0][0] if gt else "",
                     f"{exon_frac(x, strand, ts, te):.3f}"))
    with open(a.prefix + ".tsv", "w") as o:
        o.write("id\tqlen\tchr\tstart\tend\tstrand\tqcov\tmapq\tclass\tgene\tgene_type\texon_frac\n")
        for r in rows: o.write("\t".join(map(str, r)) + "\n")
    n = len(names)
    n_un = a.n_unmapped or n   # all unmapped consensus, also when only a --n subset was mapped
    with open(a.prefix + ".summary.tsv", "w") as o:
        o.write("class\tn\tpct_of_unmapped\tpct_of_all_consensus\n")
        for k, v in cnt.most_common():
            o.write(f"{k}\t{v}\t{100*v/n:.2f}\t{100*(v/n)*n_un/a.total_consensus:.2f}\n")
            print(f"{k}\t{v}\t{100*v/n:.1f}%", file=sys.stderr)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); sp = ap.add_subparsers(dest="cmd", required=True)
    e = sp.add_parser("extract"); e.add_argument("fasta"); e.add_argument("per_read"); e.add_argument("out")
    e.add_argument("--n", type=int, default=0, help="0 = all unmapped"); e.add_argument("--seed", type=int, default=42)
    c = sp.add_parser("classify"); c.add_argument("paf"); c.add_argument("gtf"); c.add_argument("prefix")
    c.add_argument("fasta", help="the fasta that was mapped (defines the denominator)")
    c.add_argument("--total_consensus", type=int, required=True); c.add_argument("--n_unmapped", type=int, default=0, help="unmapped before sampling; 0 = the mapped fasta"); c.add_argument("--min_cov", type=float, default=0.8)
    a = ap.parse_args()
    extract(a) if a.cmd == "extract" else classify(a)
