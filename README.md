# gff2tsv

Suite of scripts to convert refseq GFF to TSV using CDS.

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
# get all CDS with prelimenary filtering
python gff2tsv.py ${gff_file} [ -f ${flank} -o ${output_name} ]

# remove duplicated exons from output of gff2tsv.py and save it in new file
./remove_duplicated_exons.sh ${output_from_gff2tsv.py}
```

## Outputs

TSV file with the CDS

Format:

```tsv
chrom   start   end HGNC:ID Refseq_transcript_id    exon_nb
```
