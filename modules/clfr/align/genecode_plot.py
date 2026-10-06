#!/usr/bin/env python3
"""Bar plot of the genecode_classify.py composition TSV: % of reads and % of UMIs per class.

Classified RNA classes share one hue; the non-classifiable bins (ambiguous, low_cov, unmapped)
are neutral gray so they read as "not assigned" rather than as another RNA class.
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


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--tsv', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--title', default='')
    a = ap.parse_args()
    plot_composition(a.tsv, a.out, a.title)
