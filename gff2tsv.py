import argparse
from pathlib import Path
import sqlite3

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
    "NC_000022.10": "22", "NC_000023.10": "X"
}


def parse_gff(gff):
    """ Parse the gff data

    Args:
        gff (str): Path to GFF

    Returns:
        FeatureDB: FeatureDB object for the gff
    """

    # try to create sqlite db
    try:
        db = gffutils.create_db(
            gff, "VEP_refseq.sqlite", verbose=True,
            merge_strategy="create_unique"
        )
    except sqlite3.OperationalError as e:
        # use existing db
        db = gffutils.FeatureDB("VEP_refseq.sqlite")

    return db


def get_parents2features(db, feature_type):
    """ Get a dict of parent of feature type with its children

    Args:
        db (gffutils.FeatureDB): FeatureDB object containing gff data
        feature_type (str): Feature type to get parents for

    Returns:
        dict: Dict of parents to feature of specified type
    """

    print(f"Getting {feature_type}...")

    parents2features = {}

    for feature in db.features_of_type(feature_type):
        # check if feature has multiple parents
        if len(feature.attributes["Parent"]) != 1:
            print(feature)
            continue

        # filter out some CDS because their parents were 
        # C_gene_segment/V_gene_segment, causing the transcript names to be
        # stuff like TRGV8 and TRDC
        if db[feature.attributes["Parent"][0]].featuretype != "mRNA":
            continue

        # filter features that don't have one HGNC id or have contig chrom
        if filter_out_features(feature):
            parent = feature.attributes["Parent"][0]
            parents2features.setdefault(parent, []).append(feature)

    return parents2features


def infer_exon_number(parents2cds, parents2exons):
    """ Infer exon number for cds

    Args:
        parents2cds (dict): Parent2cds dict
        parents2exons (dict): Parent2exon dict

    Returns:
        dict: cds2exon dict
    """

    print("Infering exon number for the CDS...")

    cds_w_exon_nb = {}

    for parent in parents2cds:
        for cds in parents2cds[parent]:
            for exon in parents2exons[parent]:
                # overlap exon: ---------
                #         cds : -------
                #         cds :    --------
                #         cds :        ------
                if cds.start >= exon.start and cds.start <= exon.end:
                    exon_nb = exon.id.split("-")[-1]

                # exon:    ---------
                # cds :      -------
                # cds :  -------
                # cds : -----
                elif cds.end <= exon.end and cds.end >= exon.start:
                    exon_nb = exon.id.split("-")[-1]

                else:
                    exon_nb = None

                if exon_nb is not None:
                    cds_w_exon_nb.setdefault(cds, []).append(exon)

            # check if exons for cds and check that there is only one
            if cds in cds_w_exon_nb and len(cds_w_exon_nb[cds]) != 1:
                print("multiple exons in cds")
                print(cds)

                for e in cds_w_exon_nb[cds]:
                    print(e)
                exit()

    return cds_w_exon_nb


def filter_out_features(feature):
    """ Filter out features that have:
        - contig chrom
        - have no or multiple HGNC ids

    Args:
        feature (gffutils.Feature): Feature object

    Returns:
        bool: Whether feature get filtered or not
    """

    # check if the feature chrom is a refseq chrom
    if feature.chrom not in refseq_chrom:
        return False

    # get the hgnc id from the attributes column
    hgnc_list = [
        i
        for i in feature.attributes["Dbxref"]
        if "HGNC" in i
    ]

    # if there is none or multiple we don't want it
    if len(hgnc_list) != 1:
        return False

    return True


def write_tsv(db, data, gff, flank, output_name=None):
    """ Write tsv

    Args:
        db (gffutils.FeatureDB): FeatureDB object
        data (dict): cds2exon dict
        gff (str): Name of gff file for default name of output tsv
        flank (int): Flank to add to regions. Default is 0
        output_name (str, optional): Output name. Defaults to None.
    """

    if not output_name:
        path = Path(gff)
        name = str(path.name).replace(".gff", "").replace(".gz", "")
        output_name = f"{name}.tsv"

    print(f"Writing in {output_name}...")

    with open(output_name, "w") as f:
        for feature in data:
            hgnc_list = [
                i
                for i in feature.attributes["Dbxref"]
                if "HGNC" in i
            ]
            # format of the hgnc id in attributes column: HGNC:HGNC:1
            # get the only element in the hgnc list
            # split on ":" i.e. ["HGNC", "HGNC:1"]
            # get the last element i.e. get the actual HGNC id
            hgnc_id = hgnc_list[0].split(":", 1)[-1]

            # get the parent id and extract transcript name from it
            parent = db[feature.attributes["Parent"][0]]
            transcript = parent.id.split("-")[1]

            feature_nb = data[feature][0].id.split("-")[-1]

            f.write(
                f"{refseq_chrom[feature.chrom]}\t"
                f"{feature.start - 1 - flank}\t{feature.end + flank}\t"
                f"{hgnc_id}\t{transcript}\t{feature_nb}\n"
            )


def main(gff, flank, output_name):
    gff_db = parse_gff(gff)
    parents2exons = get_parents2features(gff_db, "exon")
    parents2cds = get_parents2features(gff_db, "CDS")
    cds_exon_nb = infer_exon_number(parents2cds, parents2exons)
    write_tsv(gff_db, cds_exon_nb, gff, flank, output_name)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("gff", help="Refseq GFF file to parse")
    parser.add_argument(
        "-f", "--flank", default=0, type=int, help="Flank to add the features"
    )
    parser.add_argument(
        "-o", "--output_name", help="Name of the output tsv file"
    )
    args = parser.parse_args()
    main(args.gff, args.flank, args.output_name)
