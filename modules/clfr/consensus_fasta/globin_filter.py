"""
Drop globin consensus fragments from consensus.fixRC.fasta.

A consensus is called globin when its best hit in consensus.paf (most matching
residues) is a transcript whose gene symbol is in --genes. The gene symbol of a
target is read from the headers of the reference the PAF was made against
(GENCODE pipe-delimited or RefSeq "... (HBB), mRNA"), because minimap2's PAF
only carries the accession.

Records with no hit in the PAF are kept (nothing to classify them with).
"""
import argparse
import gzip
import re
import sys


def open_text(path):
    return gzip.open(path, 'rt') if path.endswith('.gz') else open(path, 'r')


def read_fasta(path):
    """Yield (description, id, sequence); id is the first whitespace token."""
    desc, chunks = None, []
    with open_text(path) as f:
        for line in f:
            line = line.rstrip('\n')
            if line.startswith('>'):
                if desc is not None:
                    yield desc, desc.split()[0], ''.join(chunks)
                desc, chunks = line[1:], []
            elif desc is not None:
                chunks.append(line.strip())
    if desc is not None:
        yield desc, desc.split()[0], ''.join(chunks)


def gene_of_header(header):
    """header: FASTA header without '>'. Returns gene symbol or None."""
    if '|' in header:
        # GENCODE transcripts: ENST|ENSG|OTTHUMG|OTTHUMT|HBB-201|HBB|628|protein_coding|
        fields = header.split('|')
        if len(fields) > 5 and fields[5]:
            return fields[5]
    m = re.search(r'gene_symbol:(\S+)|gene_name:(\S+)|\[gene=([^\]]+)\]', header)
    if m:
        return next(g for g in m.groups() if g)
    # RefSeq: "NM_000518.5 Homo sapiens hemoglobin subunit beta (HBB), mRNA"
    # a description can hold several "(...)," groups, e.g. "kinase (cAMP-dependent), ... (PRKACA), ..."; the symbol is the last
    groups = re.findall(r'\(([^()\s]+)\),', header)
    if groups:
        return groups[-1]
    return None


def globin_targets(ref, genes):
    """{target name as minimap2 prints it: gene} for targets whose gene is in genes."""
    targets = {}
    with open_text(ref) as f:
        for line in f:
            if not line.startswith('>'):
                continue
            header = line[1:].rstrip('\n')
            gene = gene_of_header(header)
            if gene in genes:
                targets[header.split()[0]] = gene
    return targets


def globin_queries(paf, targets):
    """Query ids whose best PAF hit (most matching residues) is a globin target."""
    best = {}  # query -> (n_match, target)
    with open(paf) as f:
        for line in f:
            info = line.rstrip('\n').split('\t')
            if len(info) < 12:
                continue
            query, target, n_match = info[0], info[5], int(info[9])
            if query not in best or n_match > best[query][0]:
                best[query] = (n_match, target)
    return {q: targets[t] for q, (_, t) in best.items() if t in targets}, len(best)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--fasta', required=True, help='consensus.fixRC.fasta')
    p.add_argument('--paf', required=True, help='consensus.paf (consensus vs transcriptome)')
    p.add_argument('--ref', required=True, help='transcriptome FASTA the PAF was made against')
    p.add_argument('--genes', required=True, help='comma-separated gene symbols')
    p.add_argument('--out_fasta', required=True)
    p.add_argument('--out_ids', required=True, help='removed ids and their gene')
    p.add_argument('--out_stats', required=True)
    args = p.parse_args()

    genes = {g.strip() for g in args.genes.split(',') if g.strip()}
    targets = globin_targets(args.ref, genes)
    if not targets:
        sys.exit(f'ERROR: no header in {args.ref} resolved to any of {sorted(genes)}; '
                 'is --ref the transcriptome used for consensus.paf (not a genome)?')
    removed, n_in_paf = globin_queries(args.paf, targets)

    n_total = n_kept = n_removed = 0
    len_kept, len_removed = [], []
    with open(args.out_fasta, 'w') as out, open(args.out_ids, 'w') as ids:
        for desc, rid, seq in read_fasta(args.fasta):
            n_total += 1
            if rid in removed:
                n_removed += 1
                len_removed.append(len(seq))
                ids.write(f'{rid}\t{removed[rid]}\n')
            else:
                n_kept += 1
                len_kept.append(len(seq))
                out.write(f'>{desc}\n{seq}\n')

    def mean(x):
        return sum(x) / len(x) if x else 0

    with open(args.out_stats, 'w') as f:
        f.write(f'genes\t{",".join(sorted(genes))}\n')
        f.write(f'globin_reference_targets\t{len(targets)}\n')
        f.write(f'consensus_in_paf\t{n_in_paf}\n')
        f.write(f'consensus_total\t{n_total}\n')
        f.write(f'removed\t{n_removed}\t{n_removed / n_total:.4f}\n' if n_total else 'removed\t0\t0\n')
        f.write(f'kept\t{n_kept}\n')
        f.write(f'mean_len_removed\t{mean(len_removed):.1f}\n')
        f.write(f'mean_len_kept\t{mean(len_kept):.1f}\n')


if __name__ == '__main__':
    main()
