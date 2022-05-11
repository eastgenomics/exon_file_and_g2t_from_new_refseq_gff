# gff2tsv

Script to convert refseq GFF to TSV using CDS.

Info on the file:

- coordinates get converted from 1-based to 0-based format.
- exon number is infered using cds and exon coordinates (any overlap size)

CDS not gathered:

- NT_ --> Genomic Contig or scaffold, clone-based or WGSa
- NW_ --> Genomic Contig or scaffold, primarily WGSa
- Removing Y chromosome CDSes because of transcripts spanning X and Y chromosomes causing exons to be duplicated:
  - NM_176786.2
  - NM_178129.5
- Removing duplicated exons:
  - NM_001134939.1
  - NM_001172437.2
  - NM_001184961.1
  - NM_001301020.1
  - NM_001301302.1
  - NM_001301371.1
  - NM_002537.3
  - NM_004152.3
  - NM_015068.3
  - NM_016178.2
- Mitochondrial genes:
  - NC_012920.1

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
