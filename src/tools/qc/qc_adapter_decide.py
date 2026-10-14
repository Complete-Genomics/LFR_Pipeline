#!/usr/bin/env python3
"""Evidence -> adapter trim keys, for trim.sh. Three small subcommands (all print their evidence):

  qc_adapter_decide.py split   --reads test_split_read.2.fq.gz            # adapter_ref_5p, trim_5p_restrict_left, adapter_ref_3p
  qc_adapter_decide.py mgi     --reads test_split_read_2_trimmed.fastq.gz  # trimmed with 5p/3p only -> trim_mgi_adapter, adapter_ref_k
  qc_adapter_decide.py polya   --reads test_split_read_2_trimmed.fastq.gz  # trimmed with the final adapter keys -> trim_polya
  qc_adapter_decide.py setcfg  --config in.yaml --out out.yaml key=value ...   # copy a config with params keys set
  qc_adapter_decide.py params  --config c.yaml --kv a.kv b.kv ... --trim_q Q --trim_maq M --out f   # the full key block

Rules (the thresholds are arguments; they are rules of thumb, the evidence lines let you overrule them):
  split : adapter183 (GAGACGTTCTCGACTC) anchored at the read start (hit starting before --anchor_max, default 100)
          in >= --min_pct (3%) of the reads -> adapter_ref_5p true, trim_5p_restrict_left = where it ends in the 95th
          percentile read + k(17) + 5, rounded up to 10 -- deliberately the SMALL end: the pass only needs one full
          17-mer of the adapter inside the first restrict_left bases (measured on the 694 split reads: 20 does
          nothing, >= 25 works, 25-40 equivalent), and a larger window can only cut more real sequence, while the
          goal is long trimmed reads (split reads: start 5 -> 30; raw reads: start 30 -> 60);
          linker GTATCTGAGTCC anywhere in >= 3% -> adapter_ref_3p true. Reads here are SPLIT reads (the ones entering
          trim_reads): positions in raw reads differ by the barcode offset.
  mgi   : after the 5p/3p passes, MGI forward / reverse-RC fragments still in >= 3% of the reads -> trim_mgi_adapter
          true with adapter_ref_k 18 (the k that matches the 31 bp SE600 variant); else false.
  polya : trimmed reads ending in >= 10 A (or starting with >= 10 T) in >= 1% -> trim_polya 10, else 0.
"""
import argparse
import gzip
import math
import sys

A183 = 'GAGACGTTCTCGACTC'      # adapter183 (25 bp GAGACGTTCTCGACTCAGCAGAGGG), 16-mer is enough to find it
A183_LEN = 25
KMER = 17                # k of the 5' pass (bbduk k=17): one full k-mer must lie inside restrict_left
MARGIN = 5
LINKER = 'GTATCTGAGTCC'
MGI_SIGS = {'MGI fwd mid': 'CCAAGCGGTC', 'MGI fwd GGAAGACAA': 'GGAAGACAA', 'rev-RC junction': 'CAACTCCTTGGCTCAC'}


def seqs(path, n):
    opener = gzip.open if str(path).endswith('.gz') else open
    out = []
    with opener(path, 'rt') as fh:
        for i, line in enumerate(fh):
            if i % 4 == 1:
                out.append(line.strip())
                if len(out) >= n:
                    break
    return out


def pct(x, n):
    return 100.0 * x / n if n else 0.0


def cmd_split(a):
    r = seqs(a.reads, a.n)
    n = len(r)
    starts = [x.find(A183) for x in r]
    anchored = sorted(p for p in starts if 0 <= p < a.anchor_max)
    anywhere = sum(p >= 0 for p in starts)
    p5 = pct(len(anchored), n)
    link = sum(LINKER in x for x in r)
    out = {}
    print('evidence ({} split reads):'.format(n))
    print('  adapter183: {:.1f}% anchored within the first {} bp ({:.1f}% anywhere)'.format(p5, a.anchor_max, pct(anywhere, n)))
    if p5 >= a.min_pct:
        p95 = anchored[int(0.95 * (len(anchored) - 1))]
        restrict = int(math.ceil((p95 + KMER + MARGIN) / 10.0) * 10)
        print('    start (0-based) median {} p95 {} -> needs start+{} = {} (+{} margin) -> trim_5p_restrict_left {}'.format(
            anchored[len(anchored) // 2], p95, KMER, p95 + KMER, MARGIN, restrict))
        out.update({'adapter_ref_5p': 'true', 'trim_5p_restrict_left': restrict})
    else:
        print('    below {:.0f}%: 5\' pass off'.format(a.min_pct))
        out.update({'adapter_ref_5p': 'false'})
    print('  3\' linker GTATCTGAGTCC: {:.1f}% of reads'.format(pct(link, n)))
    out['adapter_ref_3p'] = 'true' if pct(link, n) >= a.min_pct else 'false'
    return out


def cmd_mgi(a):
    r = seqs(a.reads, a.n)
    n = len(r)
    print('evidence ({} reads after the 5p/3p passes):'.format(n))
    any_hit = 0
    for name, k in MGI_SIGS.items():
        c = sum(k in x for x in r)
        print('  {:<18} {:.1f}%'.format(name, pct(c, n)))
    any_hit = sum(any(k in x for k in MGI_SIGS.values()) for x in r)
    print('  any MGI fragment   {:.1f}%'.format(pct(any_hit, n)))
    on = pct(any_hit, n) >= a.min_pct
    return {'trim_mgi_adapter': 'true' if on else 'false', 'adapter_ref_k': 18} if on else {'trim_mgi_adapter': 'false'}


def cmd_polya(a):
    r = seqs(a.reads, a.n)
    n = len(r)
    tails = sum(x.endswith('A' * 10) or x.startswith('T' * 10) for x in r)
    print('evidence ({} reads after the adapter passes): ending in >=10 A / starting with >=10 T: {:.2f}%'.format(n, pct(tails, n)))
    return {'trim_polya': 10 if pct(tails, n) >= a.min_pct else 0}


def cmd_params(a):
    """All trim keys, decided ones from the kv files, the rest from the config or the pipeline default."""
    import yaml
    cfg = (yaml.safe_load(open(a.config)) or {}).get('params', {})
    se = str(cfg.get('sequence_type', 'se')).lower() == 'se'
    decided = {}
    for f in a.kv:
        for line in open(f):
            if '=' in line:
                k, v = line.strip().split('=', 1)
                decided[k] = yaml.safe_load(v)
    decided['trim_q'], decided['trim_maq'] = a.trim_q, a.trim_maq
    default = {'trim_5p_restrict_left': 30, 'adapter_ref_k': 23, 'trim_qtrim': True, 'trim_minlen': 10,
               'trim_maxns': 0 if se else -1}
    order = ['adapter_ref_5p', 'trim_5p_restrict_left', 'adapter_ref_3p', 'trim_mgi_adapter', 'adapter_ref_k',
             'trim_polya', 'trim_q', 'trim_maq', 'trim_qtrim', 'trim_minlen', 'trim_maxns']
    unused = {'trim_5p_restrict_left': decided.get('adapter_ref_5p') is False,
              'adapter_ref_k': decided.get('trim_mgi_adapter') is False}
    fmt = lambda v: str(v).lower() if isinstance(v, bool) else str(v)
    lines = ['    # trim parameters for config.yaml (' + a.label + ')']
    for k in order:
        if k in decided and not unused.get(k):
            v, src = decided[k], 'from data'
        elif k in cfg:
            v, src = cfg[k], 'from your config.yaml'
        else:
            v, src = default[k], 'pipeline default'
        note = src + (', unused here' if unused.get(k) else '') + (' -- ' + a.reason if k == 'trim_maq' and a.reason else '')
        lines.append('    {}: {}    # {}'.format(k, fmt(v), note))
    open(a.out, 'w').write('\n'.join(lines) + '\n')
    return {}


def cmd_setcfg(a):
    import yaml
    cfg = yaml.safe_load(open(a.config))
    for kv in a.kv:
        k, v = kv.split('=', 1)
        cfg.setdefault('params', {})[k] = yaml.safe_load(v)
    yaml.safe_dump(cfg, open(a.out, 'w'))
    return {}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest='cmd')
    sub.required = True     # add_subparsers(required=) needs Python 3.7
    for name, fn, mp in (('split', cmd_split, 3.0), ('mgi', cmd_mgi, 3.0), ('polya', cmd_polya, 1.0)):
        p = sub.add_parser(name)
        p.add_argument('--reads', required=True)
        p.add_argument('-n', type=int, default=20000)
        p.add_argument('--min_pct', type=float, default=mp)
        p.add_argument('--anchor_max', type=int, default=100)
        p.add_argument('--yaml_out', default='', help='write the decided keys as a yaml params snippet')
        p.set_defaults(fn=fn)
    p = sub.add_parser('params')
    p.add_argument('--config', required=True)
    p.add_argument('--kv', nargs='+', required=True)
    p.add_argument('--trim_q', type=int, required=True)
    p.add_argument('--trim_maq', type=int, required=True)
    p.add_argument('--reason', default='')
    p.add_argument('--label', default='')
    p.add_argument('--out', required=True)
    p.set_defaults(fn=cmd_params)
    p = sub.add_parser('setcfg')
    p.add_argument('--config', required=True)
    p.add_argument('--out', required=True)
    p.add_argument('kv', nargs='*')
    p.set_defaults(fn=cmd_setcfg)
    a = ap.parse_args()
    res = a.fn(a)
    if res:
        print('=> ' + ' '.join('{}: {}'.format(k, v) for k, v in res.items()))
        if getattr(a, 'yaml_out', ''):
            open(a.yaml_out, 'w').write(''.join('{}={}\n'.format(k, v) for k, v in res.items()))


if __name__ == '__main__':
    main()
