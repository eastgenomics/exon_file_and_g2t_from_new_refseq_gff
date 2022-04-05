# gff2bed
Script to convert refseq gff to bed

## How to run

```bash
python gff2bed.py ${gff_file} [ -f ${flank} -o ${output_name} ]
```

## Outputs

BED file with the exons

Format:

```tsv
chrom   start   end HGNC:ID Refseq_transcript_id    exon_nb
```
