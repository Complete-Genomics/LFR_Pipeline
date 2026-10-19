#!/usr/bin/env python3
"""
Per-chromosome UMI salvage for the consensus step.

Reads in UMIs with < --min_reads reads (low UMIs) are moved into a valid UMI (>= --min_reads reads on
this chromosome) when BOTH hold:
  1. the UMI sequences differ by <= --max_mm substitutions (an N counts as a mismatch), and
  2. the read lies inside the valid UMI's fragment [min start, max end] (+/- --pad; 0 = never extend).
A read with more than one valid UMI at its best mismatch count is left alone (ambiguous).
Without the position check about 60% of such matches are random collisions (see memory/consensus_salvage.md).

--module map    name-sorted BAM + split_stat_read1.log -> table of reads to move + stats
--module apply  BAM + table -> BAM with the first and last '#' field (the UMI) rewritten; the caller
                must re-sort by name (consensus_fasta.py groups consecutive reads by the last '#' field)

Read name: <UMI>#<pos>#<readid>#<UMI>[/2]
"""
import argparse
import itertools
import re
import subprocess
import sys

import numpy as np
import pandas as pd

L = 15
SH = [2 * (L - 1 - i) for i in range(L)]
LUT = np.full(256, 4, np.uint8)
for _i, _b in enumerate(b"ACGT"):
    LUT[_b] = _i
CIG = re.compile(r"(\d+)([MIDNSHP=X])")


def log(msg):
    print(msg, file=sys.stderr, flush=True)


def encode(mat):
    code = np.zeros(len(mat), np.int64)
    for i in range(L):
        code |= mat[:, i].astype(np.int64) << SH[i]
    return code


def xor_masks(max_mm):
    masks = []
    for k in range(1, max_mm + 1):
        for pos in itertools.combinations(range(L), k):
            for ds in itertools.product((1, 2, 3), repeat=k):
                masks.append((sum(d << SH[p] for p, d in zip(pos, ds)), k))
    return masks


def load_bam(samtools, bam):
    """Primary mapped records: names, UMI (first '#' field), start, end."""
    names, umi, pos, end = [], [], [], []
    p = subprocess.Popen([samtools, "view", "-F", "2308", bam], stdout=subprocess.PIPE,
                         text=True, bufsize=1 << 20)
    for line in p.stdout:
        f = line.split("\t", 6)
        span = sum(int(k) for k, op in CIG.findall(f[5]) if op in "MDN=X")
        names.append(f[0])
        umi.append(f[0].split("#", 1)[0])
        pos.append(int(f[3]))
        end.append(int(f[3]) + span - 1)
    if p.wait() != 0:
        sys.exit(f"samtools view failed on {bam}")
    return names, np.array(umi, f"S{L}"), np.array(pos, np.int64), np.array(end, np.int64)


def log_counts(path, umis):
    """Genome-wide reads per UMI from split_stat_read1.log ('<idx>\\t<reads>\\t<UMI>', 4 header lines)."""
    want = set(umis)
    out = {}
    rd = pd.read_csv(path, sep="\t", header=None, skiprows=4, usecols=[1, 2], names=["c", "u"],
                     dtype={"c": np.int64, "u": str}, chunksize=5_000_000)
    for ch in rd:
        m = ch["u"].isin(want)
        for c, u in zip(ch.loc[m, "c"], ch.loc[m, "u"]):
            out[u] = out.get(u, 0) + c
    return np.array([out.get(u, 0) for u in umis], np.int64)


def spans(inv, pos, end):
    order = np.lexsort((pos, inv))
    first = np.r_[0, np.flatnonzero(np.diff(inv[order])) + 1]
    return (np.minimum.reduceat(pos[order], first), np.maximum.reduceat(end[order], first),
            inv[order][first])


def find_hits(codes, present, max_mm):
    r_all, c_all, k_all = [], [], []
    for m, k in xor_masks(max_mm):
        nb = codes ^ m
        idx = np.flatnonzero(present[nb])
        if idx.size:
            r_all.append(idx)
            c_all.append(nb[idx])
            k_all.append(np.full(idx.size, k, np.int8))
    if not r_all:
        z = np.zeros(0, np.int64)
        return z, z, z.astype(np.int8)
    return np.concatenate(r_all), np.concatenate(c_all), np.concatenate(k_all)


def assign(r, v, k):
    """Per read: lowest mismatch count, then the valid UMI and how many valid UMIs tie at that count."""
    if r.size == 0:
        z = np.zeros(0, np.int64)
        return z, z, z, z
    order = np.lexsort((v, k, r))
    r, v, k = r[order], v[order], k[order].astype(np.int64)
    grp = r * 4 + k
    gi = np.flatnonzero(np.r_[True, grp[1:] != grp[:-1]])
    n = np.diff(np.r_[gi, len(grp)])
    gr = r[gi]
    first = np.r_[True, gr[1:] != gr[:-1]]
    sel = gi[first]
    return r[sel], v[sel], k[sel], n[first]


def run_map(a):
    names, ru, pos, end = load_bam(a.samtools, a.input_bam)
    u, inv, bam_cnt = np.unique(ru, return_inverse=True, return_counts=True)
    log(f"{a.chrom}: {len(ru)} primary reads, {len(u)} UMIs")
    stats = dict(chrom=a.chrom, primary_reads=len(ru), valid_umis=0, low_reads=0, rescued_reads=0,
                 rescued_mm1=0, rescued_mm2=0, ambiguous_skipped=0, chance_rescued=0)
    moves = []
    if len(u):
        lg = log_counts(a.split_log, [x.decode() for x in u])
        mat = LUT[u.view(np.uint8).reshape(-1, L)]
        has_n = (mat == 4).any(axis=1)
        code = encode(np.where(mat == 4, 0, mat))
        valid = (lg >= a.min_reads) & (bam_cnt >= a.min_reads) & ~has_n
        low = (lg < a.min_reads) & ~has_n
        vmin, vmax, vu = spans(inv, pos, end)
        vs = valid[vu]
        v_umi = vu[vs]
        o = np.argsort(code[v_umi])
        v_umi, v_code, v_min, v_max = v_umi[o], code[v_umi][o], vmin[vs][o], vmax[vs][o]
        lowread = np.flatnonzero(low[inv])
        stats["valid_umis"] = len(v_code)
        stats["low_reads"] = len(lowread)
        if len(v_code) and len(lowread):
            present = np.zeros(4 ** L, bool)
            present[v_code] = True
            lc, lp, le = code[inv[lowread]], pos[lowread], end[lowread]
            r, nb, k = find_hits(lc, present, a.max_mm)
            vi = np.searchsorted(v_code, nb)

            def pick(p_, e_):
                ok = (e_[r] >= v_min[vi] - a.pad) & (p_[r] <= v_max[vi] + a.pad)
                return assign(r[ok], vi[ok], k[ok])

            rr, vv, kk, nn = pick(lp, le)
            uniq = nn == 1
            stats["ambiguous_skipped"] = int((~uniq).sum())
            stats["rescued_reads"] = int(uniq.sum())
            stats["rescued_mm1"] = int((uniq & (kk == 1)).sum())
            stats["rescued_mm2"] = int((uniq & (kk == 2)).sum())
            # chance level: shuffle positions among the low reads and repeat
            perm = np.random.default_rng(1).permutation(len(lowread))
            cr, _, _, cn = pick(lp[perm], le[perm])
            stats["chance_rescued"] = int((cn == 1).sum())
            for ri, vj in zip(rr[uniq], vv[uniq]):
                g = lowread[ri]
                moves.append((names[g], u[inv[g]].decode(), u[v_umi[vj]].decode()))
            kk_u = kk[uniq]
            moves = [m + (int(x),) for m, x in zip(moves, kk_u)]
    with open(a.output_reads, "w") as fh:
        fh.write("qname\told_umi\tnew_umi\tmismatch\n")
        for q, o_, n_, k_ in moves:
            fh.write(f"{q}\t{o_}\t{n_}\t{k_}\n")
    cols = list(stats)
    with open(a.output_stats, "w") as fh:
        fh.write("\t".join(cols) + "\n" + "\t".join(str(stats[c]) for c in cols) + "\n")
    log(f"{a.chrom}: rescued {stats['rescued_reads']} reads (chance {stats['chance_rescued']}), "
        f"{stats['ambiguous_skipped']} ambiguous skipped")


def rename(qname, old, new):
    parts = qname.split("#")
    parts[0] = new
    if len(parts) > 1 and parts[-1].startswith(old):
        parts[-1] = new + parts[-1][len(old):]
    return "#".join(parts)


def run_apply(a):
    import pysam
    tab = pd.read_csv(a.reads_table, sep="\t", dtype=str)
    new_of = {q: (o, n) for q, o, n in zip(tab["qname"], tab["old_umi"], tab["new_umi"])}
    n_renamed = 0
    with pysam.AlignmentFile(a.input_bam, "rb") as inp, \
            pysam.AlignmentFile(a.output_bam, "wb", template=inp) as out:
        for read in inp:
            hit = new_of.get(read.query_name)
            if hit:
                read.query_name = rename(read.query_name, *hit)
                n_renamed += 1
            out.write(read)
    log(f"renamed {n_renamed} records ({len(new_of)} reads in table)")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--module", required=True, choices=["map", "apply"])
    ap.add_argument("--input_bam", required=True)
    ap.add_argument("--chrom", default="")
    ap.add_argument("--split_log")
    ap.add_argument("--min_reads", type=int, default=25)
    ap.add_argument("--max_mm", type=int, default=2, choices=[1, 2])
    ap.add_argument("--pad", type=int, default=0)
    ap.add_argument("--samtools", default="samtools")
    ap.add_argument("--output_reads")
    ap.add_argument("--output_stats")
    ap.add_argument("--reads_table")
    ap.add_argument("--output_bam")
    a = ap.parse_args()
    if a.module == "map":
        if not (a.split_log and a.output_reads and a.output_stats):
            ap.error("map needs --split_log --output_reads --output_stats")
        run_map(a)
    else:
        if not (a.reads_table and a.output_bam):
            ap.error("apply needs --reads_table --output_bam")
        run_apply(a)


if __name__ == "__main__":
    main()
