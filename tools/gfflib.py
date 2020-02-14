import csv

import gffutils
import fastalib as fl

MARGIN_SIZE = 150
DONOR = 'GT'
ACCEPTOR = 'AG'


def transform_func(d):
    try:
        d['transcriptId'] = d['proteinId']
    except KeyError:
        pass
    return d


def parse_gff(gff_file: str):
    fn = gffutils.example_filename(gff_file)

    db = gffutils.create_db(
        fn, ":memory:",
        id_spec={'transcript': 'transcriptId', 'gene': 'name'},
        gtf_transcript_key='transcriptId',
        gtf_gene_key='name',
        transform=transform_func
    )

    return db


def extract_exon_positions(gff_file, assembly_fasta, output_csv: str, validate=False) -> None:
    """
    Extracts exon coordinates into a csv with columns [scaffold, start, end]
    Only from positive strand
    :param validate: Validate positions
    :param gff_file: GFF file that contains fungi gene annotations (exons and their positions)
    :param assembly_fasta: FASTA file with whole fungi genome. Serves only as validation of extraction
    :param output_csv: CSV file with exon positions (scaffold and location)
    :return: None
    """
    db = parse_gff(gff_file)

    if validate:
        with open(assembly_fasta, 'r') as f:
            scaffold_seq_dict = {desc: seq for desc, seq in fl.read_fasta(f)}

    exons = []
    for gene in db.all_features(featuretype='gene'):
        for f in db.children(gene['name'][0], featuretype='exon'):
            # Extract the exon described in the GFF file from the given FASTA.
            exon_seq = f.sequence(assembly_fasta)

            scaffold = f.chrom
            strand = f.strand

            if strand == "+":
                if validate:
                    scaffold_seq = scaffold_seq_dict[scaffold]
                    assert scaffold_seq[f.start - 1: f.end] == exon_seq

                exons.append([scaffold, f.start - 1, f.end])

    with open(output_csv, "w") as csv_file:
        writer = csv.writer(csv_file, delimiter=';')
        writer.writerow(['scaffold', 'start', 'end'])
        writer.writerows(exons)
