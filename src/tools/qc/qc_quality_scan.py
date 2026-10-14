#!/usr/bin/env python3
"""
usage:
qc_quality_scan.py raw.fq.gz --raw --config config.yaml [--raw_r2 R2.fq.gz] [--read 1|2]

Per-batch quality-trim scan: measure quality BEFORE and AFTER candidate trims and print the config.yaml values.

Companion of qc_adapter_scan.sh (adapters). Standalone diagnostic: it does not touch trim_reads and is not part of
any workflow; its output is a config snippet you paste into the run's config.yaml.

  qc_quality_scan.py reads.fq.gz [-o outdir] [-n 100000] [--config config.yaml] [--trimmed trimmed.fq.gz]
                     [--baseline other_batch/metrics.json]

reads.fq.gz = the reads that go INTO trim_reads (data/split_read.2.fq.gz; for PE run it once per read file).
Takes the first N reads (default 100000; --stride k keeps every k-th of the first N*k reads), then:

  1. BEFORE: quality of the sample (by position, per read, end quality, Q<10 fraction, quality cliff).
  2. Runs the real bbduk quality pass on the sample for a grid of candidates (qtrim=rl trimq=Q [maq=M]; Q:M pairs
     in --candidates, default 6:0 = what the pipeline did so far, then 8:0 10:0 15:0 6:10 6:15 10:15) and measures
     each: reads kept, bases kept, mean Q, reads with mean Q<20, reads ending below Q6, Q<10 fraction.
     NB bbduk maq is an expected-error average (a few Q2 bases pull a read far below its arithmetic mean Q), so
     maq 10-15 is already demanding; maq 20 discards most reads of a batch with decaying 3' quality.
  3. Recommends the candidate that minimises the Q<10 base fraction among those that lose at most
     --max_extra_base_loss percentage points of bases (default 5) and --max_extra_read_loss points of reads
     (default 2) MORE than the current setting (first candidate), unless it improves Q<10 by less than 10%
     (then the current setting is kept). Relative to the current setting because the loss of the baseline trim
     itself depends on the batch.
  4. If --trimmed is given (the pipeline's actual trim_reads output for the same reads) it is measured as well.
  5. LOW-QUALITY BASES per read: for every read the % of its bases with Q < --low_q (default 10), then the distribution
     over all reads (mean/median/p90/p99 and the share of reads in 0%, 0-5%, 5-10%, 10-20%, 20-50%, >50%), before and
     after trim; plot low_quality_distribution.png.
  6. Writes suggested_config.yaml with trim_q / trim_maq (and the adapter keys, to be set from qc_adapter_scan.sh).

RAW MODE (--raw --config config.yaml): the input is the RAW fastq (SE: the one read file; PE: R1 plus --raw_r2).
The script copies the first N raw reads into a scratch run directory and lets the pipeline itself do the work with
the run's config.yaml: split_reads (its own split settings), then trim_reads for every candidate with
trim_q/trim_maq overridden (all other trim keys exactly as in the config). So BEFORE = the split reads that go
into trim_reads, and every AFTER number is the real pipeline trim (adapters + quality + polyA as configured), not
a quality-only approximation. The first candidate is the config's current trim_q:trim_maq. When you are happy,
copy trim_q/trim_maq from suggested_config.yaml into config.yaml and run the pipeline.

NAMING: the scratch pipeline_run/ uses the pipeline's own names (data/split_read.*), because the real rules are run.
The QC products you keep are written with a test_ prefix, so they are not confused with the full-data
split_read.* files: outdir/test_split_read.{1,2}.fq.gz (the split sample) and
outdir/test_split_read_<read>_trimmed.fastq.gz (the sample trimmed with the config's current trim settings).
--split_only stops right after the split.

Outputs in outdir: report.txt, metrics.json, quality_scan.png, low_quality_distribution.png, suggested_config.yaml, sample.fq.gz.
Quality is only one side of the trade-off: this cannot see what the trim does downstream. Validate a changed setting on
the aligned data (low_cov fraction, soft-clip length) before adopting it. --baseline compares the BEFORE metrics with
another batch and flags large deviations (re-decide the trim when the batch differs).
Needs numpy + matplotlib (+ pyyaml for --config) and bbduk (--bbduk, default 'bbduk.sh -Xmx2g').
"""
import argparse
import gzip
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np

PHRED_OFFSET = 33
CLIFF_SPAN = 20
Q_END = 6


def read_sample(path, n, stride, out_fq):
    """First n reads (every stride-th of the first n*stride) -> out_fq (gz); returns phred arrays."""
    opener = gzip.open if path.endswith('.gz') else open
    quals = []
    with opener(path, 'rt') as fh, gzip.open(out_fq, 'wt') as out:
        rec = []
        for i, line in enumerate(fh):
            rec.append(line)
            if i % 4 != 3:
                continue
            k = i // 4
            if k >= n * stride:
                break
            if k % stride == 0:
                out.write(''.join(rec))
                quals.append(np.frombuffer(rec[3].rstrip('\n').encode(), dtype=np.uint8) - PHRED_OFFSET)
            rec = []
    return quals


def read_fastq(path):
    """(phred arrays, sequences) of a whole fastq(.gz)."""
    path = str(path)
    opener = gzip.open if path.endswith('.gz') else open
    quals, seqs = [], []
    with opener(path, 'rt') as fh:
        for i, line in enumerate(fh):
            if i % 4 == 1:
                seqs.append(line.rstrip('\n'))
            elif i % 4 == 3:
                quals.append(np.frombuffer(line.rstrip('\n').encode(), dtype=np.uint8) - PHRED_OFFSET)
    return quals, seqs


LOW_Q = 10          # a base is "low quality" below this Phred value (--low_q)
LOWQ_BINS = (('0%', None, 0), ('0-5%', 0, 5), ('5-10%', 5, 10), ('10-20%', 10, 20), ('20-50%', 20, 50), ('>50%', 50, 100))


def lowq_pct(quals):
    """Per read: low-quality bases (Q < LOW_Q) as a percentage of the read length."""
    return np.array([100.0 * float((q < LOW_Q).mean()) for q in quals])


def lowq_summary(lq):
    """Distribution of the per-read low-quality percentage over all reads."""
    bins = {}
    for label, lo, hi in LOWQ_BINS:
        sel = (lq == 0) if lo is None else ((lq > lo) & (lq <= hi))
        bins[label] = float(100.0 * sel.mean())
    return {'lowq_q': LOW_Q, 'lowq_pct_mean': float(lq.mean()), 'lowq_pct_median': float(np.median(lq)),
            'lowq_pct_p90': float(np.percentile(lq, 90)), 'lowq_pct_p99': float(np.percentile(lq, 99)),
            'lowq_bins_reads_pct': bins}


def metrics(quals, bases_ref, reads_ref):
    m = _metrics_base(quals, bases_ref, reads_ref)
    if quals:
        m.update(lowq_summary(lowq_pct(quals)))
    return m


def _metrics_base(quals, bases_ref, reads_ref):
    if not quals:
        return {'reads': 0, 'reads_kept_pct': 0.0, 'bases_kept_pct': 0.0}
    n = len(quals)
    lens = np.array([len(q) for q in quals])
    allq = np.concatenate(quals)
    rmean = np.array([q.mean() for q in quals])
    last10 = np.array([q[-10:].mean() for q in quals if len(q) >= 30]) if (lens >= 30).any() else np.array([np.nan])
    return {
        'reads': n, 'reads_kept_pct': 100.0 * n / reads_ref, 'bases_kept_pct': 100.0 * int(lens.sum()) / bases_ref,
        'read_len_median': float(np.median(lens)), 'read_len_mean': float(lens.mean()),
        'base_q_mean': float(allq.mean()),
        'base_q_median': float(np.median(allq)), 'frac_q_lt10': float((allq < 10).mean()),
        'frac_q_ge30': float((allq >= 30).mean()), 'read_meanq_median': float(np.median(rmean)),
        'reads_meanq_lt20_pct': float(100 * (rmean < 20).mean()), 'last10_meanq': float(np.nanmean(last10)),
        'reads_end_lt6_pct': float(100 * (np.nan_to_num(last10, nan=99) < Q_END).mean()),
    }


def pos_curves(quals):
    L = max(len(q) for q in quals)
    mat = np.full((len(quals), L), np.nan)
    for i, q in enumerate(quals):
        mat[i, :len(q)] = q
    cover = (~np.isnan(mat)).sum(axis=0) >= 0.5 * len(quals)
    k = int(cover.sum())
    return np.nanmean(mat, axis=0)[:k], np.nanpercentile(mat, 10, axis=0)[:k]


def cliff(mean_curve):
    sm = np.convolve(mean_curve, np.ones(5) / 5, mode='same')
    best, pos = 0.0, -1
    for p in range(0, len(sm) - CLIFF_SPAN):
        d = sm[p] - sm[p + CLIFF_SPAN]
        if d > best:
            best, pos = float(d), p
    return best, pos + 1


def run_bbduk(bbduk, sample, out, q, maq):
    cmd = '{} in={} out={} qtrim=rl trimq={}{}'.format(bbduk, sample, out, q, ' maq=%d' % maq if maq > 0 else '')
    r = subprocess.run(cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, universal_newlines=True)
    if r.returncode != 0 or not os.path.exists(out):
        raise RuntimeError('bbduk failed: {}\n{}'.format(cmd, r.stderr[-400:]))


def head_reads(src, dst, n):
    """First n reads of a fastq(.gz) -> gz copy."""
    opener = gzip.open if str(src).endswith('.gz') else open
    with opener(src, 'rt') as fh, gzip.open(dst, 'wt') as out:
        for i, line in enumerate(fh):
            if i >= 4 * n:
                break
            out.write(line)


def run_snakemake(a, workdir, target, extra=()):
    cmd = [a.snakemake, '-s', a.smk, target, '-d', str(workdir), '--cores', str(a.cores), '--quiet'] + list(extra)
    r = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True)
    if r.returncode != 0 or not (Path(workdir) / target).exists():
        raise RuntimeError('snakemake failed ({}):\n{}'.format(' '.join(cmd), r.stdout[-1500:]))


def raw_setup(a):
    """Scratch pipeline run from the raw reads; returns (workdir, base config dict, split fastq path)."""
    import yaml
    cfg = yaml.safe_load(open(a.config))
    cfg.setdefault('samples', {})
    work = Path(a.outdir).resolve() / 'pipeline_run'
    if work.exists():
        shutil.rmtree(work)
    (work / 'fastq' / 'scan').mkdir(parents=True)
    se = str(cfg['params'].get('sequence_type', 'se')).lower() == 'se'
    head_reads(a.fastq, work / 'fastq' / 'scan' / 'scan_read_1.fq.gz', a.n_reads)
    if not se:
        if not a.raw_r2:
            sys.exit('PE config: give the second raw file with --raw_r2')
        head_reads(a.raw_r2, work / 'fastq' / 'scan' / 'scan_read_2.fq.gz', a.n_reads)
    # fq_path must be relative to the run dir: cat_read_two links it as ../<path>
    cfg['samples'].update({'fq_path': 'fastq', 'lanes': ['scan']})
    # only split_reads / trim_reads are used: the genecode composition mode would insist on its own inputs
    if str(cfg['params'].get('mrna_mapper', '')).lower() == 'minimap_genecode':
        cfg['params']['mrna_mapper'] = 'minimap2'
    # a config copied from a server may point src_dir at a path that does not exist here: use this repo (what
    # bootstrap.smk does for placeholder values)
    repo = Path(a.smk).resolve().parent.parent
    if not Path(str(cfg['params'].get('src_dir', ''))).is_dir():
        cfg['params']['src_dir'] = str(repo) + '/'
    # keys the config lacks (older configs) are filled from the repo's default config, as the workflow needs them to load
    default_cfg = repo / 'config' / 'config.yaml'
    if default_cfg.exists():
        def fill(dst, src):
            for k, v in src.items():
                if k not in dst:
                    dst[k] = v
                elif isinstance(v, dict) and isinstance(dst[k], dict):
                    fill(dst[k], v)
        fill(cfg, yaml.safe_load(open(default_cfg)) or {})
    # tool paths copied from a server config (e.g. /home/.../anaconda3/bin/python) that do not exist here: use what
    # is available locally, only for the two tools split_reads / trim_reads run
    for key, fallback in (('general_python', sys.executable), ('bbduk', 'bbduk.sh -Xmx2g')):
        val = str(cfg['params'].get(key, '') or '')
        tok = val.split()[0] if val.split() else ''
        if tok.startswith('/') and not os.path.exists(tok):
            print('NOTE: params.{} = {} does not exist here; using {} for this scan'.format(key, val, fallback),
                  file=sys.stderr)
            cfg['params'][key] = fallback
    # bbduk on macOS fails to start without an explicit heap size (its memory autodetect gives Xms > Xmx)
    bb = str(cfg['params'].get('bbduk', '') or '')
    if sys.platform == 'darwin' and bb and '-Xmx' not in bb:
        print('NOTE: params.bbduk = {} has no -Xmx; adding -Xmx2g for this scan (macOS)'.format(bb), file=sys.stderr)
        cfg['params']['bbduk'] = bb + ' -Xmx2g'
    # only split_reads / trim_reads are used: never the genecode composition module (which wants its own inputs)
    cfg.setdefault('modules', {})['genecode'] = False
    yaml.safe_dump(cfg, open(work / 'config.yaml', 'w'))
    return work, cfg, se


def pipeline_trim(a, base_cfg, split_dir, q, m, extra=None):
    """Run the pipeline's trim_reads on the split reads with trim_q/trim_maq overridden; returns the output fastq."""
    import yaml
    d = Path(a.outdir).resolve() / ('cand_q{}_maq{}'.format(q, m) + ''.join('_{}{}'.format(k, v) for k, v in (extra or {}).items()))
    if d.exists():
        shutil.rmtree(d)
    (d / 'data').mkdir(parents=True)
    for f in ('split_read.1.fq.gz', 'split_read.2.fq.gz'):
        os.symlink(split_dir / 'data' / f, d / 'data' / f)
    cfg = json.loads(json.dumps(base_cfg))
    cfg['params']['trim_q'], cfg['params']['trim_maq'] = q, m
    cfg['params'].update(extra or {})
    yaml.safe_dump(cfg, open(d / 'config.yaml', 'w'))
    out = 'data/split_read_{}_trimmed.fastq.gz'.format(a.read)
    run_snakemake(a, d, out, ['--allowed-rules', 'trim_reads'])
    return d, d / out


def default_smk():
    """workflows/clfr.smk of the repo this script lives in (found by walking up, so the tool can be moved)."""
    for p in Path(__file__).resolve().parents:
        if (p / 'workflows' / 'clfr.smk').exists():
            return str(p / 'workflows' / 'clfr.smk')
    return 'clfr.smk'


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('fastq')
    ap.add_argument('-o', '--outdir', default='qc_quality_scan_out')
    ap.add_argument('-n', '--n_reads', type=int, default=100000)
    ap.add_argument('--stride', type=int, default=1)
    ap.add_argument('--bbduk', default='bbduk.sh -Xmx2g')
    ap.add_argument('--candidates', default='6:0,8:0,10:0,15:0,6:10,6:15,10:15',
                    help='trimq:maq pairs; the first one is treated as the current/default setting')
    ap.add_argument('--max_extra_base_loss', type=float, default=5.0,
                    help='max extra %% points of bases lost vs the current setting (default 5)')
    ap.add_argument('--max_extra_read_loss', type=float, default=2.0,
                    help='max extra %% points of reads lost vs the current setting (default 2)')
    ap.add_argument('--config', default='', help='run config.yaml: show its current trim keys next to the suggestion')
    ap.add_argument('--trimmed', default='', help="the pipeline's trim_reads output for the same reads (optional)")
    ap.add_argument('--baseline', default='')
    ap.add_argument('--no_plot', action='store_true')
    ap.add_argument('--low_q', type=int, default=10, help='a base with Phred Q below this is low quality (default 10)')
    ap.add_argument('--raw', action='store_true', help='input is RAW fastq: split + trim through the pipeline (needs --config)')
    ap.add_argument('--raw_r2', default='', help='PE raw mode: the second raw read file')
    ap.add_argument('--read', default='2', choices=['1', '2'], help='raw mode: which split read to scan (default 2)')
    ap.add_argument('--smk', default=default_smk())
    ap.add_argument('--snakemake', default='snakemake')
    ap.add_argument('--cores', type=int, default=4)
    ap.add_argument('--split_only', action='store_true', help='raw mode: stop after writing test_split_read.*')
    a = ap.parse_args()

    globals()['LOW_Q'] = a.low_q
    os.makedirs(a.outdir, exist_ok=True)
    sample = os.path.join(a.outdir, 'sample.fq.gz')
    src_fastq = a.fastq
    work = base_cfg = None
    if a.raw:
        if not a.config:
            sys.exit('--raw needs --config config.yaml')
        work, base_cfg, _ = raw_setup(a)
        try:
            run_snakemake(a, work, 'data/split_read.{}.fq.gz'.format(a.read))
        except RuntimeError as e:
            sys.exit(str(e))
        src_fastq = str(work / 'data' / 'split_read.{}.fq.gz'.format(a.read))
        for i in ('1', '2'):
            f = work / 'data' / 'split_read.{}.fq.gz'.format(i)
            if f.exists():
                shutil.copy(f, os.path.join(a.outdir, 'test_split_read.{}.fq.gz'.format(i)))
        if a.split_only:
            print('split sample written: {}'.format(os.path.join(a.outdir, 'test_split_read.{}.fq.gz'.format(a.read))))
            return
    quals = read_sample(src_fastq, a.n_reads if not a.raw else 10 ** 9, a.stride, sample)
    if not quals:
        sys.exit('no reads in {}'.format(a.fastq))
    n = len(quals)
    bases = int(sum(len(q) for q in quals))
    before = metrics(quals, bases, n)
    mean0, p10_0 = pos_curves(quals)
    cdrop, cpos = cliff(mean0)
    before.update({'file': os.path.basename(a.fastq), 'cliff_drop_q': cdrop, 'cliff_start_pos': cpos,
                   'pos_meanq_first50': float(mean0[:50].mean()), 'pos_meanq_last100': float(mean0[-100:].mean())})

    cands = []
    if a.raw:   # the config's current setting first, then the grid
        cands.append((int(base_cfg['params'].get('trim_q', 6)), int(base_cfg['params'].get('trim_maq', 0))))
    for tok in a.candidates.split(','):
        q, m = (int(x) for x in tok.split(':'))
        if (q, m) not in cands:
            cands.append((q, m))
    results, curves, notes, rmeans, lqs = {}, {}, [], {}, {}
    for q, m in cands:
        out = os.path.join(a.outdir, 'cand_q{}_maq{}.fq.gz'.format(q, m))
        cdir = None
        try:
            if a.raw:
                cdir, out = pipeline_trim(a, base_cfg, work, q, m)
            else:
                run_bbduk(a.bbduk, sample, out, q, m)
        except (RuntimeError, OSError) as e:
            notes.append(str(e))
            break
        cq, cseq = read_fastq(out)
        if a.raw and (q, m) == cands[0]:
            cur_seqs = cseq
            shutil.copy(str(out), os.path.join(a.outdir, 'test_split_read_{}_trimmed.fastq.gz'.format(a.read)))
        results['{}:{}'.format(q, m)] = metrics(cq, bases, n)
        curves['{}:{}'.format(q, m)] = pos_curves(cq) if cq else None
        rmeans['{}:{}'.format(q, m)] = [x.mean() for x in cq]
        lqs['{}:{}'.format(q, m)] = lowq_pct(cq)
        if cdir is not None:
            shutil.rmtree(cdir)
        else:
            os.remove(out)

    if a.raw and not results:      # the trim never ran: fail loudly instead of writing an empty report
        sys.exit('ERROR: the trim_reads run failed, no trimmed sample was produced:\n' + '\n'.join(notes))
    trimmed = None
    if a.raw and results:
        k0 = '{}:{}'.format(*cands[0])
        trimmed = dict(results[k0])
        trimmed['polyA_tail_pct'] = float(100 * np.mean([x.endswith('A' * 10) or x.startswith('T' * 10) for x in cur_seqs]))
        curves['trimmed'] = curves.get(k0)
        rmeans['trimmed'] = rmeans.get(k0)
        lqs['trimmed'] = lqs.get(k0)
    if a.trimmed:
        tq, tseq = read_fastq(a.trimmed)
        trimmed = metrics(tq, bases, n)
        trimmed['reads_kept_pct'] = 100.0 * len(tq) / n          # only meaningful when --trimmed covers the same reads
        trimmed['polyA_tail_pct'] = float(100 * np.mean([s.endswith('A' * 10) or s.startswith('T' * 10) for s in tseq]))
        curves['trimmed'] = pos_curves(tq) if tq else None
        rmeans['trimmed'] = [x.mean() for x in tq]
        lqs['trimmed'] = lowq_pct(tq)

    # trim_qtrim on vs off at the current trim_q:trim_maq (the other trim keys as in the config): raw mode only
    qtrim_cmp = None
    cur_key = '{}:{}'.format(*cands[0])
    if a.raw and results:
        cfg_on = str(base_cfg['params'].get('trim_qtrim', True)).strip().lower() in ('true', '1', 'yes')
        try:
            cdir, out = pipeline_trim(a, base_cfg, work, cands[0][0], cands[0][1], {'trim_qtrim': not cfg_on})
            oq, _ = read_fastq(out)
            shutil.rmtree(cdir)
            other = metrics(oq, bases, n)
            qtrim_cmp = {'on': results[cur_key] if cfg_on else other, 'off': other if cfg_on else results[cur_key]}
        except (RuntimeError, OSError) as e:
            notes.append('trim_qtrim comparison skipped: ' + str(e)[:300])

    # recommendation
    cur_key = '{}:{}'.format(*cands[0])
    rec_key, why = cur_key, 'no candidate could be run' if not results else 'current setting kept'
    if results:
        cur = results[cur_key]
        ok = [k for k, r in results.items()
              if r['bases_kept_pct'] >= cur['bases_kept_pct'] - a.max_extra_base_loss
              and r['reads_kept_pct'] >= cur['reads_kept_pct'] - a.max_extra_read_loss]
        best = min(ok, key=lambda k: results[k]['frac_q_lt10'])
        gain = 1 - results[best]['frac_q_lt10'] / max(cur['frac_q_lt10'], 1e-9)
        if best != cur_key and gain >= 0.10:
            rec_key = best
            why = ('lowest Q<10 fraction among candidates losing <= {:g} extra points of bases and <= {:g} of reads '
                   'vs {}: {:.0%} fewer Q<10 bases'.format(a.max_extra_base_loss, a.max_extra_read_loss, cur_key, gain))
        else:
            why = ('no candidate within +{:g} points of base loss / +{:g} of read loss improves Q<10 by >= 10%: '
                   'keep {}'.format(a.max_extra_base_loss, a.max_extra_read_loss, cur_key))
    rq, rm = (int(x) for x in rec_key.split(':'))

    cfg_now = {}
    if a.config:
        import yaml
        p = (yaml.safe_load(open(a.config)) or {}).get('params', {})
        cfg_now = {k: p.get(k, '(not set)') for k in
                   ('trim_q', 'trim_maq', 'trim_qtrim', 'trim_mgi_adapter', 'adapter_ref_5p', 'trim_5p_restrict_left',
                    'adapter_ref_3p', 'adapter_ref_k', 'trim_polya', 'trim_minlen', 'trim_maxns')}

    snippet = ['    # quality trim of trim_reads, from qc_quality_scan.py on {} (n={})'.format(before['file'], n),
               '    trim_q: {}'.format(rq), '    trim_maq: {}'.format(rm),
               '    # adapter keys are NOT decided by quality; set them from qc_adapter_scan.sh / the trimmed-read check:',
               '    #   adapter_ref_5p: true/false   (5\' fixed contaminant, e.g. adapter183)',
               '    #   adapter_ref_3p: true/false   (3\' linker read-through, SE600)',
               '    #   adapter_ref_k: 18-23          (18 for SE600)',
               '    #   trim_polya: 10                (0 = off; use 10 when many trimmed reads still end in polyA)',
               '    #   trim_mgi_adapter / trim_qtrim / trim_5p_restrict_left / trim_minlen / trim_maxns: per-batch switches',
               '    #   (they replace the hand-written trim.sh; unchanged here, shown above under "current config.yaml")']
    if trimmed is not None and trimmed['polyA_tail_pct'] >= 1.0:
        snippet.append('    # NOTE: {:.1f}% of the --trimmed reads still end in >=10 A (or start with >=10 T): '
                       'set trim_polya: 10'.format(trimmed['polyA_tail_pct']))
    open(os.path.join(a.outdir, 'suggested_config.yaml'), 'w').write('params:\n' + '\n'.join(snippet) + '\n')

    json.dump({'before': before, 'candidates': results, 'recommended': rec_key, 'trimmed': trimmed, 'trim_qtrim': qtrim_cmp},
              open(os.path.join(a.outdir, 'metrics.json'), 'w'), indent=1)

    L = ['file {}  reads {} (first {} x stride {})  bases {}  read length median {:.0f}'.format(
        before['file'], n, a.n_reads, a.stride, bases, before['read_len_median']),
        '', 'BEFORE trim: base Q mean {:.1f} median {:.0f} | Q<10 {:.1%} | Q>=30 {:.1%}'.format(
            before['base_q_mean'], before['base_q_median'], before['frac_q_lt10'], before['frac_q_ge30']),
        '  read mean Q median {:.1f}; reads with mean Q<20: {:.1f}%'.format(
            before['read_meanq_median'], before['reads_meanq_lt20_pct']),
        '  mean Q by position: first 50 bp {:.1f} -> last 100 bp {:.1f}; end (last 10 bases) mean Q {:.1f}, '
        'reads ending below Q{}: {:.2f}%'.format(before['pos_meanq_first50'], before['pos_meanq_last100'],
                                                 before['last10_meanq'], Q_END, before['reads_end_lt6_pct']),
        '  biggest drop: {:.1f} Q within {} bp starting at position {}{}'.format(
            cdrop, CLIFF_SPAN, cpos, '  (sharp: look for a cassette/adapter boundary or a cycle problem)' if cdrop >= 8 else '')]
    if results:
        L += ['', 'AFTER trim ({}):'.format('real pipeline trim_reads with your config, only trim_q:trim_maq varied; '
                                   'bases/reads kept are vs the split reads' if a.raw else
                                   'bbduk qtrim=rl trimq=Q maq=M, quality only'),
              '  Q:M    reads kept  bases kept  mean Q  Q<10   reads meanQ<20  end<Q6   read len (median)']
        for k, r in results.items():
            L.append('  {:<6} {:>9.1f}%  {:>9.1f}%  {:>6.1f}  {:>5.1%}  {:>13.1f}%  {:>5.2f}%  {:>7.0f}{}'.format(
                k, r['reads_kept_pct'], r['bases_kept_pct'], r['base_q_mean'], r['frac_q_lt10'],
                r['reads_meanq_lt20_pct'], r['reads_end_lt6_pct'], r['read_len_median'],
                '   <- current/default' if k == cur_key else ('   <- recommended' if k == rec_key else '')))
    if trimmed is not None and not a.raw:
        L += ['', 'pipeline output (--trimmed, adapters + quality + polyA as configured): reads kept {:.1f}%  bases kept {:.1f}%  '
              'mean Q {:.1f}  Q<10 {:.1%}  reads meanQ<20 {:.1f}%  end<Q6 {:.2f}%  still ending in polyA/T: {:.1f}%'.format(
                  trimmed['reads_kept_pct'], trimmed['bases_kept_pct'], trimmed['base_q_mean'], trimmed['frac_q_lt10'],
                  trimmed['reads_meanq_lt20_pct'], trimmed['reads_end_lt6_pct'], trimmed['polyA_tail_pct'])]
    rows = [('before trim', before)]
    if results:
        rows.append(('after {} (current)'.format(cur_key), results[cur_key]))
        if rec_key != cur_key:
            rows.append(('after {} (recommended)'.format(rec_key), results[rec_key]))
    if trimmed is not None and not a.raw:
        rows.append(('pipeline output (--trimmed)', trimmed))
    L += ['', 'LOW-QUALITY BASES per read: bases with Q<{} as % of the read length, distribution over reads'.format(LOW_Q),
          '  {:<30} {:>6} {:>7} {:>6} {:>6}   {}'.format('', 'mean', 'median', 'p90', 'p99',
                                                        '  '.join('{:>7}'.format(b[0]) for b in LOWQ_BINS) + '   (% of reads in each bin)')]
    for name, r in rows:
        L.append('  {:<30} {:>5.1f}% {:>6.1f}% {:>5.1f}% {:>5.1f}%   {}'.format(
            name, r['lowq_pct_mean'], r['lowq_pct_median'], r['lowq_pct_p90'], r['lowq_pct_p99'],
            '  '.join('{:>6.1f}%'.format(r['lowq_bins_reads_pct'][b[0]]) for b in LOWQ_BINS)))
    if qtrim_cmp is not None:
        L += ['', 'trim_qtrim on vs off (trim_q:trim_maq = {}, other trim keys as in the config; bases/reads kept vs the split reads):'.format(cur_key),
              '  {:<12} {:>11} {:>9} {:>11} {:>11} {:>7} {:>7} {:>13}'.format(
                  '', 'reads kept', 'mean len', 'median len', 'bases kept', 'mean Q', 'Q<10', 'lowq/read')]
        for name in ('on', 'off'):
            r = qtrim_cmp[name]
            L.append('  {:<12} {:>10.2f}% {:>9.0f} {:>11.0f} {:>10.1f}% {:>7.1f} {:>6.1%} {:>12.1f}%'.format(
                'qtrim ' + name, r['reads_kept_pct'], r['read_len_mean'], r['read_len_median'], r['bases_kept_pct'],
                r['base_q_mean'], r['frac_q_lt10'], r['lowq_pct_mean']))
        on, off = qtrim_cmp['on'], qtrim_cmp['off']
        L.append('  => qtrim off keeps {:+.0f} bp mean length per read, {:+.2f} points of reads and {:+.1f} points of bases; '
                 'mean Q {:+.1f}, Q<10 {:+.1f} points'.format(
                     off['read_len_mean'] - on['read_len_mean'], off['reads_kept_pct'] - on['reads_kept_pct'],
                     off['bases_kept_pct'] - on['bases_kept_pct'], off['base_q_mean'] - on['base_q_mean'],
                     100 * (off['frac_q_lt10'] - on['frac_q_lt10'])))
    L += ['', 'recommendation: trim_q: {}  trim_maq: {}   ({})'.format(rq, rm, why)] + ['  ! ' + s for s in notes]
    if cfg_now:
        L += ['', 'current config.yaml: ' + ', '.join('{}={}'.format(k, v) for k, v in cfg_now.items())]
    L += ['', 'config.yaml (also in suggested_config.yaml):'] + snippet
    L += ['', 'caution: quality is one side only. A changed trim_q/trim_maq must be checked on the aligned data '
              '(low_cov fraction, soft-clip length) before adopting.']
    if a.baseline:
        b = json.load(open(a.baseline))['before']
        L += ['', 'vs baseline {}:'.format(a.baseline)]
        flag = False
        for k, tol in (('base_q_mean', 2.0), ('read_meanq_median', 3.0), ('last10_meanq', 3.0),
                       ('reads_meanq_lt20_pct', 5.0), ('reads_end_lt6_pct', 2.0), ('frac_q_lt10', 0.03)):
            d = before[k] - b[k]
            flag |= abs(d) > tol
            L.append('  {:<22} {:>8.3f} vs {:>8.3f}  {}'.format(k, before[k], b[k], '<-- DEVIATES' if abs(d) > tol else ''))
        L.append('  => ' + ('this batch differs from the baseline: re-decide the quality trim.'
                            if flag else 'within tolerance: the baseline trim settings can be reused.'))
    rep = '\n'.join(L)
    open(os.path.join(a.outdir, 'report.txt'), 'w').write(rep + '\n')
    print(rep)

    if not a.no_plot and results:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 3, figsize=(16, 4.4))
        ax[0].plot(mean0, color='#9aa0a6', label='before (mean)')
        ax[0].plot(p10_0, color='#9aa0a6', ls=':', label='before (10th pct)')
        for key, color, lab in ((rec_key, '#2b6fdb', 'after, recommended ' + rec_key),):
            c = curves.get(key)
            if c is not None:
                ax[0].plot(c[0], color=color, label=lab + ' (mean)')
                ax[0].plot(c[1], color=color, ls=':', label=lab + ' (10th pct)')
        if trimmed is not None and curves.get('trimmed') is not None:
            ax[0].plot(curves['trimmed'][0], color='#d9822b', label='pipeline output (mean)')
        ax[0].set_ylim(0, 42)
        ax[0].set_xlabel('position in read (bp)')
        ax[0].set_ylabel('Phred Q')
        ax[0].set_title('quality by position')
        ax[0].legend(frameon=False, fontsize=7)
        rmean0 = [q.mean() for q in quals]
        ax[1].hist(rmean0, bins=np.arange(0, 43, 1), color='#9aa0a6', alpha=0.9, label='before')
        ax[1].hist(rmeans[rec_key], bins=np.arange(0, 43, 1), histtype='step', color='#2b6fdb', lw=1.6,
                   label='after ' + rec_key)
        if trimmed is not None:
            ax[1].hist(rmeans['trimmed'], bins=np.arange(0, 43, 1), histtype='step', color='#d9822b', lw=1.6,
                       label='pipeline output')
        ax[1].axvline(20, color='#444', ls='--', lw=0.8)
        ax[1].set_xlabel('mean Q of read')
        ax[1].set_ylabel('reads')
        ax[1].set_title('per-read mean quality')
        ax[1].legend(frameon=False)
        keys = list(results)
        ax[2].scatter([results[k]['bases_kept_pct'] for k in keys], [100 * results[k]['frac_q_lt10'] for k in keys],
                      color=['#2b6fdb' if k == rec_key else '#9aa0a6' for k in keys], zorder=3)
        for k in keys:
            ax[2].annotate(k, (results[k]['bases_kept_pct'], 100 * results[k]['frac_q_lt10']), fontsize=8,
                           xytext=(4, 3), textcoords='offset points')
        ax[2].axvline(results[cur_key]['bases_kept_pct'] - a.max_extra_base_loss, color='#444', ls='--', lw=0.8)
        ax[2].set_xlabel('bases kept (%)')
        ax[2].set_ylabel('bases with Q<10 (%)')
        ax[2].set_title('trim trade-off (labels = trimq:maq)')
        for s in ax:
            for side in ('top', 'right'):
                s.spines[side].set_visible(False)
        fig.suptitle('{}  (n={})'.format(before['file'], n), x=0.01, ha='left')
        fig.tight_layout()
        fig.savefig(os.path.join(a.outdir, 'quality_scan.png'), dpi=130)

        # distribution over reads of the low-quality percentage
        fig2, bx = plt.subplots(1, 2, figsize=(11, 4.2))
        edges = np.arange(0, 105, 5)
        series = [('before trim', lowq_pct(quals), '#9aa0a6'), ('after ' + cur_key + ' (current)', lqs.get(cur_key), '#d9822b')]
        if rec_key != cur_key:
            series.append(('after ' + rec_key + ' (recommended)', lqs.get(rec_key), '#2b6fdb'))
        for name, arr, color in series:
            if arr is None or len(arr) == 0:
                continue
            bx[0].hist(arr, bins=edges, weights=np.full(len(arr), 100.0 / len(arr)), histtype='step', lw=1.8, color=color, label=name)
            xs = np.sort(arr)
            bx[1].plot(xs, 100.0 * np.arange(1, len(xs) + 1) / len(xs), color=color, lw=1.8, label=name)
        bx[0].set_yscale('log')
        bx[0].set_xlabel('low-quality bases (Q<{}) in the read, % of read length'.format(LOW_Q))
        bx[0].set_ylabel('% of reads (log)')
        bx[0].set_title('per-read low-quality percentage')
        bx[1].set_xlabel('low-quality bases (Q<{}) in the read, % of read length'.format(LOW_Q))
        bx[1].set_ylabel('cumulative % of reads')
        bx[1].set_title('cumulative distribution')
        bx[0].legend(frameon=False, fontsize=8)
        for s_ in bx:
            for side in ('top', 'right'):
                s_.spines[side].set_visible(False)
        fig2.suptitle('{}  (n={})'.format(before['file'], n), x=0.01, ha='left')
        fig2.tight_layout()
        fig2.savefig(os.path.join(a.outdir, 'low_quality_distribution.png'), dpi=130)


if __name__ == '__main__':
    main()
