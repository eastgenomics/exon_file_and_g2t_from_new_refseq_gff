# Create exon file and g2t file from new Refseq GFF

## First script, gff2tsv

Script to convert refseq GFF to TSV using CDS.

Info on the file:

- coordinates get converted from 1-based to 0-based format.
- exon number is infered using cds and exon coordinates (any overlap size)

CDS not gathered:

- NT_ --> Genomic Contig or scaffold, clone-based or WGSa
- NW_ --> Genomic Contig or scaffold, primarily WGSa
- Removing duplicated exons entirely, except for exons present in X and Y chromosomes i.e. NM_000451.4 in that case keep the X copy:
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

### How to run

Tested only with a Refseq gff file

```bash
python gff2tsv.py ${gff_file} [ -f ${flank} -o ${output_name} ]
```

### Outputs

TSV file with the CDS

Format:

```tsv
chrom   start   end HGNC:ID Refseq_transcript_id    exon_nb
```

## Second script, refseq_g2t

Script to check that g2t from new refseq gff is correct and generate a g2t from exon file from gff2tsv

HGMD file was generated using the following command line on a MySQL database using a HGMD dump (project-Fz4Q15Q42Z9YjYk110b3vGYQ:file-Fz4Q46842Z9z2Q6ZBjy7jVPY):

```shell
sudo mysql -e "select markname.hgncID, concat(refcore, '.', refversion) from gene2refseq join markname on markname.gene_id=gene2refseq.hgmdID" hgmd_2020_3 | grep -v "concat(" > hgmd_transcripts.txt
```

The Haemonc file was provided to me by Aisha Dahir and contains a single column with transcripts without versions (https://github.com/eastgenomics/Haemonc_requests/tree/main/URA-97)

### How to run

```shell
python refseq_g2t.py -n ${nirvana_g2t} -r ${refseq_exon_file} [ check ${hgmd_file} ${haemonc_file} ]
```

Use the check flag to do checks:

- nirvana
- hgmd
- haemonc

### Outputs

G2T formatted file i.e.:

```tsv
HGNC:14825  NM_001005484.2  not_clinical_transcript not_canonical
HGNC:31275  NM_001005221.2  not_clinical_transcript not_canonical
HGNC:15079  NM_001005277.1  not_clinical_transcript not_canonical
HGNC:28706  NM_001385641.1  not_clinical_transcript not_canonical
HGNC:28706  NM_001385640.1  not_clinical_transcript not_canonical
HGNC:28706  NM_152486.4 not_clinical_transcript not_canonical
HGNC:24517  NM_015658.4 not_clinical_transcript not_canonical
HGNC:24023  NM_198317.3 not_clinical_transcript not_canonical
HGNC:25284  NM_032129.3 not_clinical_transcript not_canonical
HGNC:25284  NM_001367552.1  not_clinical_transcript not_canonical
HGNC:25284  NM_001160184.2  not_clinical_transcript not_canonical
HGNC:28208  NM_001291366.2  not_clinical_transcript not_canonical
```
