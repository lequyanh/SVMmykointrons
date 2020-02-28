import csv

import gffutils
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import fastalib as fl


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


def extract_exon_positions(
        gff_file,
        assembly_fasta,
        output_csv: str,
        output_fasta: str,
        validate=False
) -> None:
    """
    Extracts exon coordinates into a csv with columns [scaffold, start, end]
    Only from positive strand
    :param output_fasta: FASTA file with exon sequences and positions
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

    exon_positions = []
    exon_sequences = []
    for gene in db.all_features(featuretype='gene'):
        for f in db.children(gene['name'][0], featuretype='exon'):
            # Extract the exon described in the GFF file from the given FASTA.
            exon_seq = f.sequence(assembly_fasta)

            scaffold, strand = f.chrom, f.strand
            start, end = f.start - 1, f.end

            if strand == "+":
                if validate:
                    scaffold_seq = scaffold_seq_dict[scaffold]
                    assert scaffold_seq[start: end] == exon_seq

                exon_positions.append([scaffold, start, end])
                exon_sequences.append(
                    SeqRecord(id=' '.join([scaffold, strand, str(start), str(end)]), seq=Seq(exon_seq))
                )

            if strand == '-':
                if validate:
                    scaffold_seq = scaffold_seq_dict[scaffold]
                    assert Seq(scaffold_seq[start: end]).reverse_complement() == exon_seq

    with open(output_csv, "w") as csv_file:
        writer = csv.writer(csv_file, delimiter=';')
        writer.writerow(['scaffold', 'start', 'end'])
        writer.writerows(exon_positions)

    with open(output_fasta, "w") as out_exon_fasta:
        SeqIO.write(exon_sequences, out_exon_fasta, 'fasta')
