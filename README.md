# gff2tsv

Script to convert refseq GFF to TSV using CDS.

Info on the file:

- coordinates get converted from 1-based to 0-based format.
- exon number is infered using cds and exon coordinates (any overlap size)

CDS not gathered:

- NT_ --> Genomic Contig or scaffold, clone-based or WGSa
- NW_ --> Genomic Contig or scaffold, primarily WGSa
- Removing Y chromosome CDSes
- Removing duplicated exons

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
