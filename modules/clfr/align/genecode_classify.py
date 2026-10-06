#!/usr/bin/env python3
"""Classify reads of a minimap2 BAM (reads mapped competitively to GENCODE transcripts
+ rRNA + genome decoy) into RNA classes and report the library composition.

The BAM must be grouped by read (minimap2 output order; do NOT coordinate-sort) and
should contain unmapped records, so unmapped reads are part of the denominator.

Reference names are the combined-fasta headers: ENST|ENSG|OTTHUMG|OTTHUMT|tx|gene|len|biotype|
Class comes from --gtf (transcript_id -> gene_type, gene_name) when given, else from the
header biotype (transcript_type). Custom rRNA / genome decoy records carry the biotype in the header.
"""
import argparse
import gzip
import re
import sys
from collections import Counter, defaultdict

import pysam

GLOBIN = {'HBA1', 'HBA2', 'HBB', 'HBD', 'HBE1', 'HBG1', 'HBG2', 'HBM', 'HBQ1', 'HBZ'}
OTHER_NCRNA = {'snRNA', 'snoRNA', 'scaRNA', 'miRNA', 'misc_RNA', 'ribozyme', 'scRNA',
               'sRNA', 'vault_RNA', 'Y_RNA', 'vaultRNA', 'piRNA', 'lincRNA'}
ORDER = ['protein_coding', 'globin', 'lncRNA', 'rRNA', 'Mt_rRNA', 'tRNA', 'other_ncRNA',
         'pseudogene', 'IG_TR', 'TEC', 'other', 'genome_decoy', 'ambiguous', 'low_cov', 'unmapped']
CIGAR_QUERY = {0, 1, 7, 8}      # M I = X consume query and count as aligned
CIGAR_CLIP = {4, 5}             # S H


def to_class(biotype, gene_name):
    if gene_name in GLOBIN:
        return 'globin'
    if biotype == 'protein_coding':
        return 'protein_coding'
    if biotype == 'lncRNA':
        return 'lncRNA'
    if biotype == 'rRNA':
        return 'rRNA'
    if biotype == 'Mt_rRNA':
        return 'Mt_rRNA'
    if biotype in ('tRNA', 'Mt_tRNA'):
        return 'tRNA'
    if biotype in OTHER_NCRNA:
        return 'other_ncRNA'
    if biotype.endswith('pseudogene'):
        return 'pseudogene'
    if biotype.startswith('IG_') or biotype.startswith('TR_'):
        return 'IG_TR'
    if biotype == 'TEC':
        return 'TEC'
    if biotype == 'genome_decoy':
        return 'genome_decoy'
    return 'other'


def load_gtf(path):
    """transcript_id -> (gene_type, gene_name) from 'transcript' lines."""
    tid_re = re.compile(r'transcript_id "([^"]+)"')
    gt_re = re.compile(r'gene_type "([^"]+)"')
    gn_re = re.compile(r'gene_name "([^"]+)"')
    out = {}
    opener = gzip.open if path.endswith('.gz') else open
    with opener(path, 'rt') as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) < 9 or f[2] != 'transcript':
                continue
            t = tid_re.search(f[8])
            g = gt_re.search(f[8])
            if t and g:
                n = gn_re.search(f[8])
                out[t.group(1)] = (g.group(1), n.group(1) if n else '')
    return out


def ref_class(rname, gtf, cache):
    c = cache.get(rname)
    if c is None:
        p = rname.split('|')
        tid = p[0]
        if tid in gtf:
            biotype, gene = gtf[tid]
        else:   # custom records (rRNA, genome_decoy) or no --gtf: fall back to header
            biotype = p[7] if len(p) >= 8 else 'other'
            gene = p[5] if len(p) >= 6 else ''
        c = cache[rname] = to_class(biotype, gene)
    return c


def umi_of(qname):
    """UMI from the read name written by split_barcode_LFR.py: <id>#<barcode>#[extra] (minimap2
    strips the trailing /1 /2). The barcode is the field right after the first '#'."""
    parts = qname.split('#')
    if len(parts) < 2 or not parts[1]:
        return None
    return parts[1]


def aligned_and_total(rec):
    aligned = clip = 0
    for op, n in rec.cigartuples or []:
        if op in CIGAR_QUERY:
            aligned += n
        elif op in CIGAR_CLIP:
            clip += n
    return aligned, aligned + clip


def classify_read(recs, gtf, cache, min_cov, tie_frac):
    """recs: all records of one read. Returns (class, ambiguous_combo or None)."""
    hits = [r for r in recs if not r.is_unmapped and not r.is_supplementary]
    if not hits:
        return 'unmapped', None
    scored = []
    for r in hits:
        a, total = aligned_and_total(r)
        score = r.get_tag('AS') if r.has_tag('AS') else a
        scored.append((score, a, total, r))
    best = max(s[0] for s in scored)
    top = [s for s in scored if s[0] >= tie_frac * best]
    # read length: the primary record keeps SEQ, so its clipped length is the read length
    read_len = max(s[2] for s in scored)
    if max(s[1] for s in top) < min_cov * read_len:
        return 'low_cov', None
    classes = {ref_class(s[3].reference_name, gtf, cache) for s in top}
    # An exonic read aligns equally well to its transcript and to the genome decoy; the decoy is only
    # meant to absorb reads that fit the genome better, so a tie goes to the transcript class.
    if len(classes) > 1:
        classes.discard('genome_decoy')
    if len(classes) == 1:
        return classes.pop(), None
    return 'ambiguous', '/'.join(sorted(classes))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--bam', required=True)
    ap.add_argument('--out', required=True, help='composition TSV')
    ap.add_argument('--out_ambiguous', required=True, help='ambiguous class-combination TSV')
    ap.add_argument('--gtf', default='', help='GENCODE annotation GTF(.gz) for gene_type; optional')
    ap.add_argument('--min_cov', type=float, default=0.8,
                    help='best hit must cover >= this fraction of the read (default 0.8)')
    ap.add_argument('--tie_frac', type=float, default=0.99,
                    help='hits with AS >= tie_frac*best AS count as best hits (default 0.99)')
    a = ap.parse_args()

    gtf = load_gtf(a.gtf) if a.gtf else {}
    cache = {}
    reads = Counter()
    umis = defaultdict(set)
    all_umis = set()
    ambiguous = Counter()

    def flush(recs):
        if not recs:
            return
        cls, combo = classify_read(recs, gtf, cache, a.min_cov, a.tie_frac)
        reads[cls] += 1
        if combo:
            ambiguous[combo] += 1
        u = umi_of(recs[0].query_name)
        if u is not None:
            umis[cls].add(u)
            all_umis.add(u)

    cur, name = [], None
    with pysam.AlignmentFile(a.bam, 'rb', check_sq=False) as bam:
        for rec in bam:
            if rec.query_name != name:
                flush(cur)
                cur, name = [], rec.query_name
            cur.append(rec)
        flush(cur)

    total = sum(reads.values())
    if total == 0:
        sys.exit('no reads in {}'.format(a.bam))
    n_umi = len(all_umis)
    classes = [c for c in ORDER if c in reads] + sorted(c for c in reads if c not in ORDER)
    with open(a.out, 'w') as o:
        o.write('class\tn_reads\tpct_reads\tn_umis\tpct_umis\n')
        for c in classes:
            nu = len(umis[c])
            o.write('{}\t{}\t{:.3f}\t{}\t{:.3f}\n'.format(
                c, reads[c], 100.0 * reads[c] / total, nu, 100.0 * nu / n_umi if n_umi else 0.0))
        o.write('total\t{}\t100.000\t{}\t100.000\n'.format(total, n_umi))
    with open(a.out_ambiguous, 'w') as o:
        o.write('classes\tn_reads\tpct_reads\n')
        for combo, n in ambiguous.most_common():
            o.write('{}\t{}\t{:.3f}\n'.format(combo, n, 100.0 * n / total))


if __name__ == '__main__':
    main()
