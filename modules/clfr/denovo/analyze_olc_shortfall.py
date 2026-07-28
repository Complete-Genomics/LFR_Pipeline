#!/usr/bin/env python3
"""Classify UMI assemblies where MEGAHIT is >=1 kb but OLC is <1 kb.

The comparison is deliberately reference-free in the biological sense:
MEGAHIT is used only as a coordinate system to ask whether input reads bridge
an OLC endpoint.  A read is placed only when one orientation has at least
three co-linear exact k-mer anchors spanning 40 bp.  This is the same kind of
evidence used by the collective OLC rescue, but reported rather than acted on.
"""

import argparse
import csv
from collections import Counter, OrderedDict, defaultdict


BC_START = 5
BC_LEN = 15
RC = str.maketrans("ACGT", "TGCA")


def rc(seq):
    return seq.translate(RC)[::-1]


def barcode_from_header(header):
    return header[1:].split()[0].split(">", 1)[0]


def read_fasta_longest(path):
    result = {}
    name = None
    chunks = []
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if line.startswith(">"): 
                if name is not None:
                    seq = "".join(chunks)
                    if len(seq) > len(result.get(name, "")):
                        result[name] = seq
                name, chunks = barcode_from_header(line), []
            elif line:
                chunks.append(line)
    if name is not None:
        seq = "".join(chunks)
        if len(seq) > len(result.get(name, "")):
            result[name] = seq
    return result


def read_first_barcode_groups(path, n):
    groups = OrderedDict()
    with open(path) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2:
                continue
            barcode = fields[0][BC_START:BC_START + BC_LEN]
            if barcode not in groups:
                if len(groups) == n:
                    break
                groups[barcode] = []
            groups[barcode].append(fields[1])
    return groups


def kmer_index(seq, k):
    index = defaultdict(list)
    for pos in range(len(seq) - k + 1):
        index[seq[pos:pos + k]].append(pos)
    return index


def place(query, target, target_index, k=17, min_anchors=3, min_span=40):
    """Return (start, end, anchors, span), or None for no unambiguous placement."""
    offsets = defaultdict(list)
    for qpos in range(len(query) - k + 1):
        hits = target_index.get(query[qpos:qpos + k], ())
        if len(hits) > 12:
            continue
        for tpos in hits:
            offsets[tpos - qpos].append(qpos)
    candidates = []
    for offset, positions in offsets.items():
        unique = sorted(set(positions))
        if len(unique) >= min_anchors and unique[-1] - unique[0] >= min_span:
            candidates.append((len(unique), unique[-1] - unique[0], offset))
    if not candidates:
        return None
    candidates.sort(reverse=True)
    best = candidates[0]
    if len(candidates) > 1 and candidates[1][:2] == best[:2]:
        return None
    anchors, span, start = best
    return start, start + len(query), anchors, span


def place_either(query, target, target_index):
    placements = []
    for oriented in (query, rc(query)):
        result = place(oriented, target, target_index)
        if result is not None:
            placements.append(result)
    if len(placements) != 1:
        return None
    return placements[0]


def classify(barcode, reads, olc, megahit):
    index = kmer_index(megahit, 17)
    olc_placement = place_either(olc, megahit, index)
    row = {
        "barcode": barcode,
        "reads_raw": len(reads),
        "reads_unique": len(set(reads)),
        "megahit_len": len(megahit),
        "olc_len": len(olc),
        "gap_bp": len(megahit) - len(olc),
        "olc_start": "",
        "olc_end": "",
        "left_flank_bp": "",
        "right_flank_bp": "",
        "placed_reads": 0,
        "left_bridge_reads": 0,
        "right_bridge_reads": 0,
        "class": "olc_unplaced",
    }
    if olc_placement is None:
        return row
    start, end, _, _ = olc_placement
    row.update({
        "olc_start": start,
        "olc_end": end,
        "left_flank_bp": max(0, start),
        "right_flank_bp": max(0, len(megahit) - end),
    })
    left = set()
    right = set()
    placed = 0
    for read_id, read in enumerate(sorted(set(reads))):
        placement = place_either(read, megahit, index)
        if placement is None:
            continue
        placed += 1
        read_start, read_end, _, _ = placement
        # Require the same 20 bp minimum overlap on both sides of an endpoint.
        # A one-base overhang would otherwise look like a bridge despite being
        # unusable by the OLC overlap rule.
        if read_start <= start - 20 and read_end >= start + 20:
            left.add(read_id)
        if read_start <= end - 20 and read_end >= end + 20:
            right.add(read_id)
    row["placed_reads"] = placed
    row["left_bridge_reads"] = len(left)
    row["right_bridge_reads"] = len(right)
    left_supported = start > 0 and len(left) >= 2
    right_supported = end < len(megahit) and len(right) >= 2
    if start <= 0 and end >= len(megahit):
        row["class"] = "olc_covers_megahit"
    elif left_supported and right_supported:
        row["class"] = "bridged_both_ends"
    elif left_supported or right_supported:
        row["class"] = "bridged_one_end"
    else:
        row["class"] = "no_two_read_bridge"
    return row


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--r2", required=True)
    parser.add_argument("--megahit", required=True)
    parser.add_argument("--olc", required=True)
    parser.add_argument("--n", type=int, default=20000)
    parser.add_argument("--out", required=True)
    parser.add_argument("--replay-collective", action="store_true",
                        help="also reassemble each shortfall UMI with collective rescue off/on")
    args = parser.parse_args()

    groups = read_first_barcode_groups(args.r2, args.n)
    megahit = read_fasta_longest(args.megahit)
    olc = read_fasta_longest(args.olc)
    rows = []
    for barcode, reads in groups.items():
        if len(megahit.get(barcode, "")) < 1000 or not megahit.get(barcode):
            continue
        if len(olc.get(barcode, "")) >= 1000 or not olc.get(barcode):
            continue
        rows.append(classify(barcode, reads, olc[barcode], megahit[barcode]))

    if args.replay_collective:
        from denovo_seed_olc import assemble_umi
        for pos, row in enumerate(rows, 1):
            reads = groups[row["barcode"]]
            no_rescue = assemble_umi(reads, use_collective_rescue=False)
            with_rescue = assemble_umi(reads, use_collective_rescue=True)
            row["no_collective_len"] = max(map(len, no_rescue), default=0)
            row["collective_len"] = max(map(len, with_rescue), default=0)
            if pos % 100 == 0:
                print("replayed={}".format(pos), flush=True)

    fields = list(rows[0]) if rows else []
    with open(args.out, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    counts = Counter(row["class"] for row in rows)
    print("shortfall_umis={}".format(len(rows)))
    for name in sorted(counts):
        print("{}={}".format(name, counts[name]))
    if args.replay_collective:
        longer = sum(row["collective_len"] > row["no_collective_len"] for row in rows)
        reaches_1k = sum(row["no_collective_len"] < 1000 <= row["collective_len"]
                         for row in rows)
        print("collective_longer={}".format(longer))
        print("collective_reaches_1k={}".format(reaches_1k))


if __name__ == "__main__":
    main()
