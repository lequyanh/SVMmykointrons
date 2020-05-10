import csv
import os
import sys

import gffutils
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import config
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


def extract_intron_positions(
        gff_file,
        fungi_name: str,
) -> None:
    """
    Extracts exon coordinates into a csv with columns [scaffold, start, end]
    Only from positive strand
    :param gff_file: GFF file that contains fungi gene annotations (exons and their positions)
    :param fungi_name: Fungi name
    :return: None
    """
    db = parse_gff(gff_file)
    gene_identifier = get_gff_identifier(gff_file)

    assembly_fasta = config.get_fungi_assembly(fungi_name)
    with open(assembly_fasta, 'r') as f:
        scaffold_seq_dict = {desc: seq for desc, seq in fl.read_fasta(f)}

    intron_positions = []
    intron_sequences = []

    for gene in db.all_features(featuretype='gene'):
        intron_start = None

        for i, f in enumerate(db.children(gene[gene_identifier][0], featuretype='exon')):
            scaffold, strand = f.chrom, f.strand
            start, end = f.start - 1, f.end

            # If first exon, just set the start/end of the following intron
            if i == 0 and strand == '+':
                intron_start = end  # If first exon and + strand, the exon's end is the intron start
                continue

            if i == 0 and strand == '-':
                # If strand - however, the exon order is reversed and we start with the last exon.
                # Hence the last exon's start is the intron's end
                intron_end = start
                continue

            if strand == '+':
                intron_end = start
            elif strand == '-':
                intron_start = end

            # Extract intron
            scaffold_seq = scaffold_seq_dict[scaffold]

            if strand == "+":
                intron_seq = Seq(scaffold_seq[intron_start: intron_end])
            elif strand == '-':
                intron_seq = Seq(scaffold_seq[intron_start: intron_end]).reverse_complement()

            intron_positions.append([scaffold, strand, start, end])
            intron_sequences.append(
                SeqRecord(id=' '.join([scaffold, strand, str(start), str(end)]), seq=intron_seq)
            )

            # Move cursor to the next intron
            if strand == "+":
                intron_start = end
            elif strand == '-':
                intron_end = start

    write_output(fungi_name, intron_positions, intron_sequences, type='intron')


def extract_exon_positions(
        gff_file,
        fungi_name: str,
        validate=False
) -> None:
    """
    Extracts exon coordinates into a csv with columns [scaffold, start, end]
    Only from positive strand
    :param validate: Validate positions
    :param gff_file: GFF file that contains fungi gene annotations (exons and their positions)
    :param fungi_name: fungi abbreviation.
    :return: None
    """
    db = parse_gff(gff_file)
    gene_identifier = get_gff_identifier(gff_file)

    assembly_fasta = config.get_fungi_assembly(fungi_name)
    if validate:
        with open(assembly_fasta, 'r') as f:
            scaffold_seq_dict = {desc: seq for desc, seq in fl.read_fasta(f)}

    exon_positions = []
    exon_sequences = []
    for gene in db.all_features(featuretype='gene'):
        for f in db.children(gene[gene_identifier][0], featuretype='exon'):
            # Extract the exon described in the GFF file from the given FASTA.
            exon_seq = f.sequence(assembly_fasta)

            scaffold, strand = f.chrom, f.strand
            start, end = f.start - 1, f.end

            if strand == "+":
                if validate:
                    scaffold_seq = scaffold_seq_dict[scaffold]
                    assert scaffold_seq[start: end] == exon_seq

                exon_positions.append([scaffold, strand, start, end])
                exon_sequences.append(
                    SeqRecord(id=' '.join([scaffold, strand, str(start), str(end)]), seq=Seq(exon_seq))
                )
            # TODO Extract also negative strand exons
            if strand == '-':
                if validate:
                    scaffold_seq = scaffold_seq_dict[scaffold]
                    assert Seq(scaffold_seq[start: end]).reverse_complement() == exon_seq

    write_output(fungi_name, exon_positions, exon_sequences, type='exon')


def write_output(fungi_name: str, positions: list, sequences: list, type: str):
    out_folder = f'{config.NEWSEQUENCES_LOC}/{fungi_name}'
    if not os.path.exists(out_folder):
        os.mkdir(out_folder)

    if type == 'exon':
        csv_suffix = config.EXON_CSV_SUFF
        fasta_suffix = config.EXON_FASTA_SUFF
    elif type == 'intron':
        csv_suffix = config.INTRON_CSV_SUFF
        fasta_suffix = config.INTRON_FASTA_SUFF
    else:
        raise ValueError("Type not recognized")

    output_csv = f'{out_folder}/{fungi_name}{csv_suffix}'
    with open(output_csv, "w") as csv_file:
        writer = csv.writer(csv_file, delimiter=';')
        writer.writerow(['scaffold', 'strand', 'start', 'end'])
        writer.writerows(positions)

    output_fasta = f'{out_folder}/{fungi_name}{fasta_suffix}'
    with open(output_fasta, "w") as out_exon_fasta:
        SeqIO.write(sequences, out_exon_fasta, 'fasta')


def get_gff_identifier(gff_file):
    filename, file_extension = os.path.splitext(gff_file)
    if file_extension == '.gff':
        gene_identifier = 'name'
    else:
        assert file_extension == '.gff3'
        gene_identifier = 'ID'
    return gene_identifier


if __name__ == '__main__':
    extract_intron_positions(sys.argv[1], sys.argv[2])
