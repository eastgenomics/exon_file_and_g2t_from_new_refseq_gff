import argparse
import datetime
from collections import defaultdict


def get_datetime():
    """ Get the datetime

    Returns:
        str: Datetime when function is called
    """

    now = datetime.datetime.now()
    run_datetime = now.strftime("%y%m%d")
    return run_datetime


def parse_nirvana_g2t(nirvana_g2t):
    """Parse nirvana genes2transcripts for comparing to refseq data

    Args:
        nirvana_g2t (str): File path to nirvana g2t

    Returns:
        dict: Dict of dict for genes and transcripts in nirvana g2t
    """

    nirvana_g2t_data = defaultdict(lambda: defaultdict(int))

    with open(nirvana_g2t) as f:
        for line in f:
            gene, transcript, clinical_transcript, canonical = line.strip().split()
            nirvana_g2t_data[gene][transcript] = (clinical_transcript, canonical)

    return nirvana_g2t_data


def parse_refseq_exons(refseq_exons):
    """ Parse refseq exon file

    Args:
        refseq_exons (str): Path to refseq exons file

    Returns:
        dict: Dict of data from refseq exon file
    """

    refseq_exon_data = {}

    with open(refseq_exons) as f:
        for line in f:
            chrom, start, end, gene, transcript, exon_nb = line.strip().split()
            refseq_exon_data.setdefault(gene, set()).add(transcript)

    return refseq_exon_data


def find_discrepancies(nirvana_data, refseq_data):
    """ Find discrepancies between nirvana and refseq data

    Args:
        nirvana_data (dict): Dict of data from nirvana g2t
        refseq_data (dict): Dict of data from refseq exon file
    """

    print("Checking nirvana...")

    nb_discrepancy_genes = []
    nb_discrepancy_clinical_transcript = []
    nb_discrepancy_transcript = []
    nb_discrepancy_transcript_version = []

    for gene in nirvana_data:
        if gene not in refseq_data:
            nb_discrepancy_genes.append(gene)
        else:
            for transcript in nirvana_data[gene]:
                if transcript not in refseq_data[gene]:
                    transcript_base, transcript_version = transcript.split(".")

                    refseq_base_transcripts = [
                        tx.split(".")[0] for tx in refseq_data[gene]
                    ]

                    if transcript_base not in refseq_base_transcripts:
                        nb_discrepancy_transcript.append(transcript)

                        if nirvana_data[gene][transcript][0] == "clinical_transcript":
                            nb_discrepancy_clinical_transcript.append(transcript)
                    else:
                        nb_discrepancy_transcript_version.append(transcript)

    print(f"Nb genes missing: {nb_discrepancy_genes}")
    print(f"Nb clinical transcript missing: {nb_discrepancy_clinical_transcript}")
    print(f"Nb transcripts missing: {nb_discrepancy_transcript}")
    print(f"Nb transcript version wrong: {nb_discrepancy_transcript_version}")


def check_HGMD_in_refseq(hgmd_file, refseq_data):
    """ Check if hgmd transcripts are in refseq data

    Args:
        refseq_data (dict): Dict of data from Refseq gff
    """

    print("Checking HGMD...")

    transcripts_OK = []
    transcripts_not_OK = []

    with open(hgmd_file) as f:
        for line in f:
            gene, tx = line.strip().split()
            hgnc_id = f"HGNC:{gene}"

            if hgnc_id in refseq_data:
                if tx in refseq_data[hgnc_id]:
                    transcripts_OK.append(tx)
                else:
                    tx_base, tx_version = tx.split(".")
                    refseq_base_transcripts = [
                        tx.split(".")[0] for tx in refseq_data[hgnc_id]
                    ]

                    if tx_base not in refseq_base_transcripts:
                        transcripts_not_OK.append(tx)
                    else:
                        transcripts_OK.append(tx)

    print(f"Transcripts OK: {len(transcripts_OK)}")
    print(transcripts_OK)
    print(f"Transcripts not OK: {len(transcripts_not_OK)}")
    print(transcripts_not_OK)


def check_haemonc_transcript_in_refseq(haemonc_file, refseq_data):
    """ Check haemonc transcripts are in refseq data

    Args:
        refseq_data (dict): Dict of refseq data
    """

    print("Checking Haemonc...")

    refseq_transcripts = [tx for gene, tx in refseq_data.items()]

    flatten_refseq_base_transcripts = [
        refseq_tx.split(".")[0]
        for refseq_sets in refseq_transcripts
        for refseq_tx in refseq_sets
    ]

    transcripts_OK = []
    transcripts_not_OK = []

    with open(haemonc_file) as f:
        for line in f:
            gene, transcripts = line.strip().split("\t")
            list_transcripts = [tx.strip() for tx in transcripts.split(",")]

            for ho_tx in list_transcripts:
                if ho_tx not in flatten_refseq_base_transcripts:
                    transcripts_not_OK.append(ho_tx)
                else:
                    transcripts_OK.append(ho_tx)

    print(f"Transcripts OK: {len(transcripts_OK)}")
    print(transcripts_OK)
    print(f"Transcripts not OK: {len(transcripts_not_OK)}")
    print(transcripts_not_OK)


def write_refseq_g2t(refseq_data, nirvana_data):
    """ Write refseq g2t

    Args:
        refseq_data (dict): Dict of data from refseq GFF
        nirvana_data (dict): Dict of data from nirvana g2t file
    """

    output_name = f"{get_datetime()}_g2t.tsv"

    with open(output_name, "w") as f:
        for gene in refseq_data:
            refseq_transcripts = refseq_data[gene]

            for transcript in refseq_transcripts:
                # check for nirvana data to gather clinical transcript and
                # canonical status
                if transcript in nirvana_data[gene]:
                    clinical_transcript, canonical = nirvana_data[gene][transcript]
                else:
                    # check if it's just a version difference
                    transcript_base, transcript_version = transcript.split(".")

                    nirvana_base_transcripts = [
                        tx.split(".")[0] for tx in nirvana_data[gene]
                    ]

                    if transcript_base in nirvana_base_transcripts:
                        nirvana_transcript = sorted([
                            tx for tx in nirvana_data[gene]
                            if tx.startswith(transcript_base)
                        ], reverse=True)
                        clinical_transcript, canonical = nirvana_data[gene][nirvana_transcript[0]]
                    else:
                        clinical_transcript = "not_clinical_transcript"
                        canonical = "not_canonical"

                data = (
                    f"{gene}\t{transcript}\t{clinical_transcript}\t"
                    f"{canonical}\n"
                )

                f.write(data)


def main(args):
    nirvana_data = parse_nirvana_g2t(args.nirvana_g2t)
    refseq_data = parse_refseq_exons(args.refseq_exons)

    if args.cmd == "check":
        find_discrepancies(nirvana_data, refseq_data)
        check_HGMD_in_refseq(args.hgmd, refseq_data)
        check_haemonc_transcript_in_refseq(args.haemonc, refseq_data)

    write_refseq_g2t(refseq_data, nirvana_data)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    subparser = parser.add_subparsers(dest="cmd")

    parser.add_argument(
        "-r", "--refseq_exons", required=True, help="refseq exons file"
    )
    parser.add_argument(
        "-n", "--nirvana_g2t", required=True, help="Nirvana g2t file"
    )

    check = subparser.add_parser("check")
    check.add_argument("hgmd", help="HGMD file containing HGMD transcripts")
    check.add_argument("haemonc", help="File containing Haemonc transcripts")

    args = parser.parse_args()
    main(args)
