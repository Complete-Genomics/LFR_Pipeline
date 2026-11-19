#!/usr/bin/env python3
"""Bar plots of the genecode outputs.

composition (default): genecode_classify.py composition TSV, % of reads and % of UMIs per class. Classified RNA
classes share one hue; the non-classifiable bins (ambiguous, low_cov, unmapped) are neutral gray so they read as
"not assigned" rather than as another RNA class.

--kind genome: genecode_genome.py <prefix>_by_class.tsv, one 100% bar per transcript class (unmapped, low_cov,
genome_decoy, ...) split by where its reads land in the genome (exon / intronic / intergenic / no hit).
--metric bases weights by read length instead of counting reads (needs genecode_genome.py --fq).
"""
import argparse
import csv

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

NOT_ASSIGNED = {'ambiguous', 'low_cov', 'unmapped'}
BLUE, GRAY = '#2a6fdb', '#9aa0a6'


def read_table(path):
    rows = []
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            rows.append(r)
    total = next((r for r in rows if r['class'] == 'total'), None)
    return [r for r in rows if r['class'] != 'total'], total


def plot_composition(tsv, out, title=''):
    rows, total = read_table(tsv)
    labels = [r['class'] for r in rows]
    colors = [GRAY if c in NOT_ASSIGNED else BLUE for c in labels]
    fig, axes = plt.subplots(1, 2, figsize=(10, 0.42 * len(rows) + 1.6), sharey=True)
    for ax, (col, ncol, name) in zip(axes, [('pct_reads', 'n_reads', 'reads'), ('pct_umis', 'n_umis', 'UMIs')]):
        vals = [float(r[col]) for r in rows]
        y = range(len(rows))
        ax.barh(y, vals, color=colors, height=0.62)
        for yi, r, v in zip(y, rows, vals):
            ax.text(v + max(vals) * 0.01, yi, '{:.1f}%  ({})'.format(v, r[ncol]), va='center', fontsize=8, color='#444')
        ax.set_yticks(list(y))
        ax.set_yticklabels(labels)
        ax.set_xlim(0, max(vals) * 1.35)
        ax.set_xlabel('% of {}'.format(name) + (' (n={})'.format(total[ncol]) if total else ''))
        ax.xaxis.grid(True, color='#e5e5e5', linewidth=0.6)
        ax.set_axisbelow(True)
        for s in ('top', 'right', 'left'):
            ax.spines[s].set_visible(False)
        ax.tick_params(length=0)
    axes[0].invert_yaxis()   # shared y: once, so the table order reads top-down
    fig.suptitle(title or 'GENCODE composition', x=0.01, ha='left', fontsize=11)
    fig.tight_layout()
    fig.savefig(out, dpi=200)


# stack order = left to right; exon classes blue, intronic orange, intergenic purple, not placed gray
GENOME_COLORS = [('exon_same_strand', '#2a6fdb'), ('exon_antisense', '#8fb4ee'), ('intronic', '#e8833a'),
                 ('intergenic', '#8e5bb5'), ('genome_low_cov', '#b8bdc2'), ('genome_no_hit', '#6b7075')]


def plot_genome_by_class(tsv, out, title='', metric='reads'):
    col = 'n_bases' if metric == 'bases' else 'n_reads'
    data = {}   # transcript_class -> {genome_class: n_reads or n_bases}
    with open(tsv) as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            data.setdefault(r['transcript_class'], {})[r['genome_class']] = int(r[col])
    classes = sorted(data, key=lambda c: -sum(data[c].values()))
    fig, ax = plt.subplots(figsize=(10, 0.5 * len(classes) + 2))
    for yi, c in enumerate(classes):
        tot, left = sum(data[c].values()), 0.0
        for g, col in GENOME_COLORS:
            pct = 100.0 * data[c].get(g, 0) / tot
            ax.barh(yi, pct, left=left, color=col, height=0.66, label=g if yi == 0 else None)
            if pct >= 6:
                ax.text(left + pct / 2, yi, '{:.0f}%'.format(pct), ha='center', va='center', fontsize=8,
                        color='white' if g != 'exon_antisense' else '#222')
            left += pct
        ax.text(101, yi, (('{:.1f} Mb'.format(tot / 1e6) if tot >= 1e6 else '{:.0f} kb'.format(tot / 1e3)) if metric == 'bases' else 'n={:,}'.format(tot)), va='center', fontsize=8, color='#444')
    ax.set_yticks(range(len(classes)))
    ax.set_yticklabels(classes)
    ax.invert_yaxis()
    ax.set_xlim(0, 100)
    ax.set_xlabel('% of {} in the transcript class'.format(metric))
    for s in ('top', 'right', 'left'):
        ax.spines[s].set_visible(False)
    ax.tick_params(length=0)
    ax.legend(ncol=6, loc='lower left', bbox_to_anchor=(0, 1.0), frameon=False, fontsize=8, handlelength=1.2)
    fig.suptitle(title or 'Genome placement by transcript class', x=0.01, ha='left', fontsize=11)
    fig.tight_layout()
    fig.savefig(out, dpi=200, bbox_inches='tight')


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--tsv', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--title', default='')
    ap.add_argument('--kind', choices=['composition', 'genome'], default='composition')
    ap.add_argument('--metric', choices=['reads', 'bases'], default='reads', help='genome kind only')
    a = ap.parse_args()
    if a.kind == 'genome': plot_genome_by_class(a.tsv, a.out, a.title, a.metric)
    else: plot_composition(a.tsv, a.out, a.title)
