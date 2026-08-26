# db — reference databases (paths only, not shipped in git)

Canonical, documented *location* for the reference databases `random_inspection.smk`
(and any future QA tooling) expects, so a standalone checkout has one obvious
place to point config at — but the actual files are NOT committed here. They
are large (see below), change rarely, and git (no LFS configured in this repo)
handles neither well. This directory exists so the path convention is fixed
and documented even though the bytes are not tracked.

## Expected files

| file | size | source | used by |
|---|---:|---|---|
| `gg.fna` | ~1.78 GB | Greengenes 13_5 (`ftp://greengenes.microbio.me/greengenes_release/gg_13_5/gg_13_5.fasta.gz`) | `random_inspection.smk`'s reference-assignment step (vsearch `--usearch_global` against it, assigning each barcode a "local ground truth" reference — see denovo.md sec 59) |

## Getting `gg.fna`

    wget ftp://greengenes.microbio.me/greengenes_release/gg_13_5/gg_13_5.fasta.gz
    gunzip -c gg_13_5.fasta.gz > modules/clfr/denovo/db/gg.fna

(The copy already in use as of 2026-08 lives at
`subprojects/olc/kraken2_db/k2_greengenes/library/gg.fna` on the dev machine
this was built on — see `kraken2_db/build_greengenes.log` there for the exact
fetch. Point config at that path, or copy/symlink it here, rather than
re-downloading if it's already present.)

## Why not git-lfs

Not set up for this repo yet. If reference data needs change more than
occasionally, that's the first thing to add rather than growing this directory
with more raw binaries.
