# globin_regions.bed

GRCh38 (chr-prefixed contigs, as in `GCA_000001405.15_GRCh38_no_alt_analysis_set`) regions of the human globin genes and
globin pseudogenes. `consensus_fasta.smk: filter_globin_bam` removes reads overlapping them before consensus.
BED, 0-based half-open, columns: contig, start, end, gene.

- 10 globin genes: MANE Select transcript span (HBA1 HBA2 HBB HBD HBE1 HBG1 HBG2 HBM HBQ1 HBZ).
  Not the gene span: GENCODE gives HBE1 / HBG2 ~240 kb read-through transcripts that would also cover the
  neighbouring olfactory receptor genes.
- 3 globin pseudogenes: gene span (HBBP1 HBZP1 HBAP1).
- Built from GENCODE v50 (`gencode.v50.annotation.gtf.gz`); identical from v49. Built once, not at run time.

Rebuild (e.g. for a new GENCODE release):

    gzip -dc gencode.vNN.annotation.gtf.gz | awk -F'\t' -v OFS='\t' \
      -v genes="HBA1 HBA2 HBB HBD HBE1 HBG1 HBG2 HBM HBQ1 HBZ" -v pseudo="HBBP1 HBZP1 HBAP1" '
      BEGIN{n=split(genes,a," "); for(i=1;i<=n;i++) g[a[i]]=1; n=split(pseudo,b," "); for(i=1;i<=n;i++) p[b[i]]=1}
      match($9,/gene_name "[^"]+"/) {
        nm=substr($9,RSTART+11,RLENGTH-12)
        if ($3=="transcript" && $9 ~ /tag "MANE_Select"/ && (nm in g)) print $1,$4-1,$5,nm
        else if ($3=="gene" && $9 ~ /pseudogene/ && (nm in p)) print $1,$4-1,$5,nm}' \
      | sort -k1,1 -k2,2n -u > globin_regions.bed

Lives next to `consensus_fasta.smk`, which reads it from `<src_dir>/modules/clfr/consensus_fasta/`. Expect 13 lines. For a reference with other contig names (e.g. `11` instead of `chr11`) or another species, set
`params.consensus_exclude_bed` to your own BED.
