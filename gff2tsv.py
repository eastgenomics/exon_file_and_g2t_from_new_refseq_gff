import argparse
from pathlib import Path
import sqlite3

import gffutils
from natsort import natsorted
import pandas as pd


def set_chromosome_numbers(genome_build):
    """
    Replace the RefSeq chromosome numbers with numerical values for chromosomes
    based on the genome build.
    Args:
        genome_build (str): Genome build of RefSeq gff

    Returns:
        refseq_chrom (dict): Dictionary mapping RefSeq chromosome numbers to
                                simple numeric values for chromosomes
    """

    if genome_build == "37":
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
    elif genome_build == "38":
        refseq_chrom = {
            "NC_000001.11": "1", "NC_000002.12": "2", "NC_000003.12": "3",
            "NC_000004.12":	"4", "NC_000005.10": "5", "NC_000006.12": "6",
            "NC_000007.14": "7", "NC_000008.11": "8", "NC_000009.12": "9",
            "NC_000010.11": "10", "NC_000011.10": "11", "NC_000012.12": "12",
            "NC_000013.11": "13", "NC_000014.9": "14", "NC_000015.10": "15",
            "NC_000016.10": "16", "NC_000017.11": "17", "NC_000018.10": "18",
            "NC_000019.10": "19", "NC_000020.11": "20", "NC_000021.9": "21",
            "NC_000022.11": "22", "NC_000023.11": "X", "NC_000024.10": "Y"
        }
    else:
        print('''
            Genome build not given as '37' or '38'. Unable to map RefSeq
            chromosome numbers (e.g. NC_000001.10) to simple chromosome numbers
            (e.g. 1)
              ''')
    return refseq_chrom


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


def get_parents2features(db, feature_type, refseq_chrom):
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
        if filter_out_features(feature, refseq_chrom):
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

    print("Inferring exon number for the CDS...")

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


def filter_out_features(feature, refseq_chrom):
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


def get_transcripts_to_remove(db, data):
    """ Get the transcripts to remove because of duplicated exons

    Args:
        db (FeatureDB object): Feature DB object of the GFF
        data (dict): Dict of cds to exons

    Returns:
        list: List of transcripts to remove
    """

    print("Gathering transcripts to remove...")

    list_of_transcripts_exons = []

    for feature in data:
        # get the parent id and extract transcript name from it
        parent = db[feature.attributes["Parent"][0]]
        transcript = parent.id.split("-")[1]
        feature_nb = data[feature][0].id.split("-")[-1]

        # get all transcripts+exons
        list_of_transcripts_exons.append({
            "transcript": transcript, "exon_nb": feature_nb
        })

    if list_of_transcripts_exons:
        df = pd.DataFrame(list_of_transcripts_exons)
        # get duplicated exons and get the corresponding transcripts into a list
        transcripts_to_remove = set(df[df.duplicated()]["transcript"].to_list())
    else:
        transcripts_to_remove = []

    return transcripts_to_remove


def write_sorted_data_to_file(filename, data):
    """
    Write sorted and filtered data from the GFF to file

    Args:
        filename (str): Name of file to write out
        data (list): List containing data to write
    """
    try:
        with open(filename, "w") as outfile:
            for line in data:
                str_line = [str(l) for l in line]
                outfile.write("\t".join(str_line) + "\n")
    except IOError as e:
        print(f"Error writing to {filename}: {e}")
        raise


def write_output_tsvs(
    db, data, transcripts_to_remove, gff, flank, refseq_chrom,
    hgnc_output_name=None, symbol_output_name=None
):
    """ Write exons TSV files, one with HGNC IDs and one with gene symbols

    Args:
        db (gffutils.FeatureDB): FeatureDB object
        data (dict): cds2exon dict
        gff (str): Name of gff file for default name of output tsv
        flank (int): Flank to add to regions. Default is 0
        hgnc_output_name (str, optional): Output name for exons TSV with HGNC
        IDs. Defaults to None.
        symbol_output_name (str, optional): Output name for exons TSV with
        gene symbols. Defaults to None.
    """

    if hgnc_output_name is None:
        path = Path(gff)
        name = str(path.name).replace(".gff", "").replace(".gz", "")
        hgnc_output_name = f"{name}.exon_{flank}bp.tsv"
    if symbol_output_name is None:
        name = str(path.name).replace(".gff", "").replace(".gz", "")
        symbol_output_name = f"{name}.symbols.exon_{flank}bp.tsv"

    data_to_write = []
    symbol_data_to_write = []

    print("Sorting data...")

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
        gene_symbols = feature.attributes.get("gene")
        gene_symbol = ','.join(gene_symbols) if gene_symbols else "Unknown"

        # get the parent id and extract transcript name from it
        parent = db[feature.attributes["Parent"][0]]
        transcript = parent.id.split("-")[1]

        if transcript not in transcripts_to_remove:
            feature_nb = data[feature][0].id.split("-")[-1]

            data_to_write.append([
                refseq_chrom[feature.chrom], feature.start - 1 - flank,
                feature.end + flank, hgnc_id, transcript, feature_nb
            ])
            symbol_data_to_write.append([
                refseq_chrom[feature.chrom], feature.start - 1 - flank,
                feature.end + flank, gene_symbol, transcript, feature_nb
            ])
        else:
            # some duplicated transcripts span X and Y, final decision is
            # to keep the X copy of the transcript
            if refseq_chrom[feature.chrom] == "X":
                feature_nb = data[feature][0].id.split("-")[-1]

                data_to_write.append([
                    refseq_chrom[feature.chrom], feature.start - 1 - flank,
                    feature.end + flank, hgnc_id, transcript, feature_nb
                ])
                symbol_data_to_write.append([
                    refseq_chrom[feature.chrom], feature.start - 1 - flank,
                    feature.end + flank, gene_symbol, transcript, feature_nb
                ])

    # sort data by chrom, start, end
    sorted_data = natsorted(data_to_write, key=lambda x: (x[0], x[1], x[2]))
    sorted_symbol = natsorted(
        symbol_data_to_write, key=lambda x: (x[0], x[1], x[2])
    )

    print(f"Writing in {hgnc_output_name}...")
    write_sorted_data_to_file(hgnc_output_name, sorted_data)

    print(f"Writing in {symbol_output_name}...")
    write_sorted_data_to_file(symbol_output_name, sorted_symbol)


def main(build, gff, flank, hgnc_output_name, symbols_output_name):
    refseq_chrom = set_chromosome_numbers(build)
    gff_db = parse_gff(gff)
    parents2exons = get_parents2features(gff_db, "exon", refseq_chrom)
    parents2cds = get_parents2features(gff_db, "CDS", refseq_chrom)
    cds_exon_nb = infer_exon_number(parents2cds, parents2exons)
    transcripts_to_remove = get_transcripts_to_remove(gff_db, cds_exon_nb)
    write_output_tsvs(
        gff_db, cds_exon_nb, transcripts_to_remove, gff, flank, refseq_chrom,
        hgnc_output_name, symbols_output_name
    )

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-g", "--gff", help="Refseq GFF file to parse")
    parser.add_argument(
        "-f", "--flank", default=0, type=int, help="Flank to add the features"
    )
    parser.add_argument(
        "-o", "--hgnc_output_name", help=(
            "Name of the output TSV file with HGNC IDs"
        )
    )
    parser.add_argument(
        "-s", "--symbols_output_name", help=(
            "Name of the output TSV file with gene symbols"
        )
    )
    parser.add_argument(
        "-b", "--build", help="Genome build of RefSeq GFF"
    )
    args = parser.parse_args()
    main(
        args.build, args.gff, args.flank, args.hgnc_output_name,
        args.symbols_output_name
    )
