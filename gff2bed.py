import argparse
from pathlib import Path

import gffutils

# replace refseq chrom by "normal" chrom names
refseq_chrom = {
    "NC_000001.10": "1", "NC_000002.11": "2", "NC_000003.11": "3",
    "NC_000004.11":	"4", "NC_000005.9": "5", "NC_000006.11": "6",
    "NC_000007.13": "7", "NC_000008.10": "8", "NC_000009.11": "9",
    "NC_000010.10": "10", "NC_000011.9": "11", "NC_000012.11": "12",
    "NC_000013.10": "13", "NC_000014.8": "14", "NC_000015.9": "15",
    "NC_000016.9": "16", "NC_000017.10": "17", "NC_000018.9": "18",
    "NC_000019.9": "19", "NC_000020.10": "20", "NC_000021.8": "21",
    "NC_000022.10": "22", "NC_000023.10": "X", "NC_000024.9": "Y"
}


def parse_gff(gff, flank):
    """ Parse the gff data and get exons data

    Args:
        gff (str): Path to GFF
        flank (int): Flank to add in the exon start/end

    Returns:
        list: List of exons to write into BED format
    """

    gff_data = []

    db = gffutils.create_db(
        gff, ":memory:", verbose=True, merge_strategy="create_unique"
    )

    for exon in db.features_of_type("exon"):
        # check if the exon has a transcript_id in the attributes column
        if "transcript_id" in exon.attributes:
            # check if the exon chrom is a refseq chrom
            if exon.chrom in refseq_chrom:
                # get the hgnc id from the attributes column
                hgnc_list = [
                    i
                    for i in exon.attributes["Dbxref"]
                    if "HGNC" in i
                ]

                # if there is none or multiple we don't want it
                if len(hgnc_list) == 1:
                    # format of the hgnc id in attributes column: HGNC:HGNC:1
                    # get the only element in the hgnc list
                    # split on ":" i.e. ["HGNC", "HGNC:1"]
                    # get the last element i.e. get the actual HGNC id
                    hgnc_id = hgnc_list[0].split(":", 1)[-1]

                    transcript = ','.join(exon.attributes["transcript_id"])
                    exon_nb = exon.id.split("-")[-1]

                    gff_data.append(
                        f"{refseq_chrom[exon.chrom]}\t"
                        f"{exon.start - 1 - flank}\t{exon.end + flank}\t"
                        f"{hgnc_id}\t{transcript}\t{exon_nb}\n"
                    )
                else:
                    print(
                        f"{refseq_chrom[exon.chrom]}\t"
                        f"{exon.start}\t{exon.end}\t"
                        f"{exon.attributes['Dbxref']}"
                    )

    return gff_data


def write_bed(gff_data, gff, output_name=None):
    """ Write bed with exon data

    Args:
        gff_data (list): List with exon data
        gff (str): Path to gff file to extract default output name
        output_name (str, optional): Output name. Defaults to None.
    """

    if not output_name:
        path = Path(gff)
        name = str(path.name).replace(".gff", "").replace(".gz", "")
        output_name = f"{name}.bed"

    with open(output_name, "w") as f:
        for data in gff_data:
            f.write(data)


def main(gff, flank, output_name):
    gff_data = parse_gff(gff, flank)
    write_bed(gff_data, gff, output_name)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("gff", help="Refseq GFF file to parse")
    parser.add_argument(
        "-f", "--flank", default=0, help="Flank to add the exons"
    )
    parser.add_argument(
        "-o", "--output_name", help="Name of the output bed file"
    )
    args = parser.parse_args()
    main(args.gff, args.flank, args.output_name)
