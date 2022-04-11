# gff2tsv

Script to convert refseq GFF to TSV using CDS.

Only getting NC chromosomes i.e. `NC_` and not:

- NT_ --> Genomic Contig or scaffold, clone-based or WGSa
- NW_ --> Genomic Contig or scaffold, primarily WGSa

## How to run

Tested only with a Refseq gff file

```bash
python gff2tsv.py ${gff_file} [ -f ${flank} -o ${output_name} ]
```

## Outputs

TSV file with the CDS

Format:

```tsv
chrom   start   end HGNC:ID Refseq_transcript_id    exon_nb
```
