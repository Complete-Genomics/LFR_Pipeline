#!/usr/bin/env python3
"""Per-batch base-quality scan: decide HOW to quality-trim before choosing trim parameters.

Companion of qc_adapter_scan.sh (adapters); this one only looks at base qualities. Standalone diagnostic: it does
not touch trim_reads and is not part of any workflow.

  qc_quality_scan.py reads.fq.gz [-o outdir] [-n 100000] [--stride 1] [--baseline prev/metrics.json]

Reads the first N reads (--stride k keeps every k-th read of the first N*k reads, for a less file-order-biased
sample) and writes outdir/{report.txt, metrics.json, quality_scan.png}:

  1. mean Q by read position (+ 10th percentile, fraction Q<10): smooth decay vs cliff vs mid-read dip
  2. mean Q per read (histogram): is there a population of bad reads (drop the read) vs a bad tail (trim it)
  3. end quality: mean Q of the last 10 bases and the fraction of reads that end below Q6 -- the only thing
     bbduk qtrim=rl trimq=6 (the pipeline default) can act on
  4. trim simulation: mean bases removed / reads touched by (a) end trim (strip trailing bases < q, what bbduk
     qtrim=r does), (b) 4-base sliding window from the 5' end, cut at the first window with mean < q (fastp
     --cut_right style). If (a) removes ~nothing but (b) removes a lot, the low quality sits inside the read
     and interspersed with good bases: end trimming cannot reach it.
  5. cliff: largest drop of the smoothed mean-Q curve over a 20 bp span and where it happens (a fixed position
     points to an adapter/cassette boundary: cut there instead of quality-trimming)

--baseline compares against another batch's metrics.json and flags large deviations (cross-batch alert).
Needs numpy + matplotlib (lfr_pipeline env).
"""
import argparse
import gzip
import json
import os
import sys

import numpy as np

PHRED_OFFSET = 33
WINDOW = 4          # sliding-window size of the simulated window trim
CLIFF_SPAN = 20     # bp over which the drop of the smoothed mean-Q curve is measured
Q_END = 6           # bbduk's default trimq, i.e. what the pipeline uses today
SIM_Q = (6, 10, 15, 20)


def read_quals(path, n, stride):
    """First n reads (every stride-th of the first n*stride reads) -> list of uint8 phred arrays."""
    opener = gzip.open if path.endswith('.gz') else open
    out = []
    with opener(path, 'rt') as fh:
        for i, line in enumerate(fh):
            if i % 4 != 3:
                continue
            k = i // 4
            if k >= n * stride:
                break
            if k % stride == 0:
                out.append(np.frombuffer(line.rstrip('\n').encode(), dtype=np.uint8) - PHRED_OFFSET)
    return out


def end_trim_len(q, thr):
    """bases removed by stripping trailing bases with Q < thr."""
    good = np.nonzero(q >= thr)[0]
    return len(q) - (good[-1] + 1 if len(good) else 0)


def window_trim_len(q, thr, w=WINDOW):
    """bases removed when cutting at the first w-window (5'->3') whose mean Q < thr."""
    if len(q) < w:
        return 0
    c = np.cumsum(np.insert(q.astype(float), 0, 0))
    means = (c[w:] - c[:-w]) / w
    bad = np.nonzero(means < thr)[0]
    return len(q) - bad[0] if len(bad) else 0


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('fastq')
    ap.add_argument('-o', '--outdir', default='qc_quality_scan_out')
    ap.add_argument('-n', '--n_reads', type=int, default=100000)
    ap.add_argument('--stride', type=int, default=1)
    ap.add_argument('--baseline', default='', help='metrics.json of another batch to compare against')
    ap.add_argument('--no_plot', action='store_true')
    a = ap.parse_args()

    quals = read_quals(a.fastq, a.n_reads, a.stride)
    if not quals:
        sys.exit('no reads in {}'.format(a.fastq))
    os.makedirs(a.outdir, exist_ok=True)
    n = len(quals)
    lens = np.array([len(q) for q in quals])
    L = int(lens.max())
    mat = np.full((n, L), np.nan)
    for i, q in enumerate(quals):
        mat[i, :len(q)] = q

    pos_mean = np.nanmean(mat, axis=0)
    pos_p10 = np.nanpercentile(mat, 10, axis=0)
    pos_low = np.nanmean(mat < 10, axis=0)
    pos_n = (~np.isnan(mat)).sum(axis=0)
    ok = pos_n >= 0.5 * n                                  # positions covered by >= half of the reads
    read_mean = np.array([q.mean() for q in quals])
    last10 = np.array([q[-10:].mean() for q in quals if len(q) >= 30])

    sim = {}
    for t in SIM_Q:
        e = np.array([end_trim_len(q, t) for q in quals])
        w = np.array([window_trim_len(q, t) for q in quals])
        sim[str(t)] = {'end_mean_bases': float(e.mean()), 'end_reads_trimmed_pct': float(100 * (e > 0).mean()),
                       'win_mean_bases': float(w.mean()), 'win_reads_trimmed_pct': float(100 * (w > 0).mean())}

    # cliff: biggest drop of the smoothed curve over CLIFF_SPAN bp, among well-covered positions
    sm = np.convolve(np.nan_to_num(pos_mean), np.ones(5) / 5, mode='same')
    cliff_drop, cliff_pos = 0.0, -1
    for p in range(0, int(ok.sum()) - CLIFF_SPAN):
        d = sm[p] - sm[p + CLIFF_SPAN]
        if d > cliff_drop:
            cliff_drop, cliff_pos = float(d), p
    all_q = np.concatenate(quals)
    m = {
        'file': os.path.basename(a.fastq), 'n_reads': n,
        'read_len_median': float(np.median(lens)), 'read_len_max': L,
        'base_q_median': float(np.median(all_q)), 'base_q_mean': float(all_q.mean()),
        'frac_q_lt10': float((all_q < 10).mean()), 'frac_q_ge30': float((all_q >= 30).mean()),
        'read_meanq_median': float(np.median(read_mean)), 'reads_meanq_lt20_pct': float(100 * (read_mean < 20).mean()),
        'last10_meanq': float(last10.mean()), 'reads_end_lt6_pct': float(100 * (last10 < Q_END).mean()),
        'pos_meanq_first50': float(np.nanmean(pos_mean[:50])),
        'pos_meanq_last100': float(np.nanmean(pos_mean[:int(ok.sum())][-100:])),
        'cliff_drop_q': cliff_drop, 'cliff_start_pos': cliff_pos, 'trim_sim': sim,
    }
    json.dump(m, open(os.path.join(a.outdir, 'metrics.json'), 'w'), indent=1)

    lines = ['file {}  reads {}  (first {} x stride {})  read length median {:.0f} max {}'.format(
        m['file'], n, a.n_reads, a.stride, m['read_len_median'], L),
        'base Q: median {:.0f}  mean {:.1f}  Q<10 {:.1%}  Q>=30 {:.1%}'.format(
            m['base_q_median'], m['base_q_mean'], m['frac_q_lt10'], m['frac_q_ge30']),
        'read mean Q: median {:.1f}  reads with mean Q<20: {:.1f}%'.format(
            m['read_meanq_median'], m['reads_meanq_lt20_pct']),
        'mean Q by position: first 50 bp {:.1f} -> last 100 bp {:.1f}'.format(
            m['pos_meanq_first50'], m['pos_meanq_last100']),
        'end quality: last-10-base mean Q {:.1f}; reads ending below Q{}: {:.2f}%'.format(
            m['last10_meanq'], Q_END, m['reads_end_lt6_pct']),
        'cliff: curve drops {:.1f} Q over {} bp starting at position {}'.format(cliff_drop, CLIFF_SPAN, cliff_pos + 1),
        '', 'trim simulation (mean bases removed per read / % reads touched):',
        '   Q    end-trim (bbduk qtrim=r)    {}-bp sliding window (fastp cut_right)'.format(WINDOW)]
    for t in SIM_Q:
        s = sim[str(t)]
        lines.append('  {:>2}    {:6.1f} / {:5.1f}%           {:6.1f} / {:5.1f}%'.format(
            t, s['end_mean_bases'], s['end_reads_trimmed_pct'], s['win_mean_bases'], s['win_reads_trimmed_pct']))

    # rule-of-thumb reading (a hint, not a decision)
    s6, s20 = sim['6'], sim['20']
    hints = []
    med = max(m['read_len_median'], 1.0)
    if s6['end_mean_bases'] < 1.0 and s20['win_mean_bases'] > 10 * max(s6['end_mean_bases'], 0.1):
        hints.append('end trim at Q6 removes ~nothing ({:.1f} bp) while a Q20 window would remove {:.0f} bp/read '
                     '({:.0f}% of the median read): the low quality is spread inside the read, not at its end. '
                     'Quality-trimming that hard is not free; only adopt a window trim if the aligned data '
                     '(low_cov, clip length) improves more than it loses.'
                     .format(s6['end_mean_bases'], s20['win_mean_bases'], 100 * s20['win_mean_bases'] / med))
    if cliff_drop >= 8:
        hints.append('sharp quality drop ({:.0f} Q within {} bp, position ~{}): a fixed position suggests a '
                     'cassette/adapter boundary or a cycle problem; cut there (and run qc_adapter_scan.sh) instead '
                     'of quality-trimming.'.format(cliff_drop, CLIFF_SPAN, cliff_pos + 1))
    if m['reads_meanq_lt20_pct'] > 10:
        hints.append('{:.0f}% of reads have mean Q<20: a population of bad reads, consider dropping whole reads '
                     '(mean-Q filter) rather than trimming.'.format(m['reads_meanq_lt20_pct']))
    if m['reads_end_lt6_pct'] > 5:
        hints.append('{:.0f}% of reads end below Q6: end trimming (qtrim=r) is worthwhile.'.format(m['reads_end_lt6_pct']))
    if not hints:
        hints.append('no strong quality signal: default end trim is enough; look at adapters instead (qc_adapter_scan.sh).')
    lines += ['', 'reading:'] + ['  - ' + h for h in hints]

    if a.baseline:
        b = json.load(open(a.baseline))
        lines += ['', 'vs baseline {}:'.format(a.baseline)]
        flag = False
        for k, tol_abs, tol_rel in (('base_q_mean', 2.0, None), ('read_meanq_median', 3.0, None),
                                    ('last10_meanq', 3.0, None), ('reads_meanq_lt20_pct', 5.0, None),
                                    ('reads_end_lt6_pct', 2.0, None), ('frac_q_lt10', 0.03, None)):
            d = m[k] - b[k]
            bad = abs(d) > tol_abs
            flag |= bad
            lines.append('  {:<22} {:>8.3f} vs {:>8.3f}  {}'.format(k, m[k], b[k], '<-- DEVIATES' if bad else ''))
        lines.append('  => ' + ('this batch differs from the baseline: re-decide the quality trim.'
                                if flag else 'within tolerance: reuse the baseline trim settings.'))

    rep = '\n'.join(lines)
    open(os.path.join(a.outdir, 'report.txt'), 'w').write(rep + '\n')
    print(rep)

    if not a.no_plot:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))
        x = np.arange(1, ok.sum() + 1)
        ax[0].plot(x, pos_mean[:ok.sum()], color='#2b6fdb', label='mean')
        ax[0].plot(x, pos_p10[:ok.sum()], color='#9aa0a6', label='10th percentile')
        ax[0].set_ylim(0, 42)
        ax[0].set_xlabel('position in read (bp)')
        ax[0].set_ylabel('Phred Q')
        ax[0].set_title('quality by position')
        ax[0].legend(frameon=False)
        ax[1].hist(read_mean, bins=np.arange(0, 43, 1), color='#2b6fdb')
        ax[1].axvline(20, color='#9aa0a6', ls='--')
        ax[1].set_xlabel('mean Q of read')
        ax[1].set_ylabel('reads')
        ax[1].set_title('per-read mean quality')
        xs = np.arange(len(SIM_Q))
        ax[2].bar(xs - 0.2, [sim[str(t)]['end_mean_bases'] for t in SIM_Q], 0.4, color='#9aa0a6', label='end trim')
        ax[2].bar(xs + 0.2, [sim[str(t)]['win_mean_bases'] for t in SIM_Q], 0.4, color='#2b6fdb',
                  label='{}-bp window'.format(WINDOW))
        ax[2].set_xticks(xs)
        ax[2].set_xticklabels(['Q{}'.format(t) for t in SIM_Q])
        ax[2].set_ylabel('mean bases removed / read')
        ax[2].set_title('trim simulation')
        ax[2].legend(frameon=False)
        for s in ax:
            s.spines[['top', 'right']].set_visible(False)
        fig.suptitle('{}  (n={})'.format(m['file'], n), x=0.01, ha='left')
        fig.tight_layout()
        fig.savefig(os.path.join(a.outdir, 'quality_scan.png'), dpi=130)


if __name__ == '__main__':
    main()
