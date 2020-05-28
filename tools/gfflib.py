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


def extract_gff_introns(
        gff_file,
        fungi_name: str,
) -> None:
    db = parse_gff(gff_file)
    gene_identifier = get_gff_identifier(gff_file)

    assembly_fasta = config.get_fungi_assembly(fungi_name)
    with open(assembly_fasta, 'r') as f:
        scaffold_seq_dict = {desc: seq for desc, seq in fl.read_fasta(f)}

    intron_positions = []
    intron_sequences = []

    for gene in db.all_features(featuretype='gene'):
        intron_start = None

        print_next = False
        for i, f in enumerate(db.children(gene[gene_identifier][0], featuretype='exon')):
            scaffold, strand, start, end = \
                f.chrom, f.strand, f.start - 1, f.end

            # If first exon, just set the start/end of the following intron
            if i == 0:
                intron_start = end  # If first exon and + strand, the exon's end is the intron start
                continue
            intron_end = start

            if intron_start > intron_end:
                print(f'Anomally at strand {strand}. Next introns of the gene are:')
                print_next = True
                intron_start = end
                continue

            scaffold_seq = scaffold_seq_dict[scaffold]
            intron_seq = Seq(scaffold_seq[intron_start: intron_end])

            if intron_end - intron_start != len(intron_seq):
                print('Valid coordinates but no intron')
                print([intron_end - intron_start, len(intron_seq)])
                intron_start = end
                continue

            if strand == "+":
                assert intron_end - intron_start == len(intron_seq)

                intron_positions.append([scaffold, strand, intron_start, intron_end])
                intron_sequences.append(
                    SeqRecord(id=' '.join([scaffold, strand, str(intron_start), str(intron_end)]),
                              seq=intron_seq)
                )
            elif strand == '-':
                intron_seq = intron_seq.reverse_complement()

                # Coordinates relative to the minus strand
                intron_start_conv = len(scaffold_seq) - intron_end
                intron_end_conv = len(scaffold_seq) - intron_start

                # Verify the sequence and intron locations are correct
                intron_seq2 = Seq(scaffold_seq).reverse_complement()[intron_start_conv:intron_end_conv]
                assert intron_seq == intron_seq2
                assert intron_end_conv - intron_start_conv == len(intron_seq)

                intron_positions.append([scaffold, strand, intron_start_conv, intron_end_conv])
                intron_sequences.append(
                    SeqRecord(id=' '.join([scaffold, strand, str(intron_start_conv), str(intron_end_conv)]),
                              seq=intron_seq)
                )

            if print_next:
                print(intron_seq)

            # Move cursor to the next intron
            intron_start = end

    write_output(fungi_name, intron_positions, intron_sequences, type='intron')


def extract_gff3_introns(
        gff_file,
        fungi_name: str,
) -> None:
    db = parse_gff(gff_file)
    gene_identifier = get_gff_identifier(gff_file)

    assembly_fasta = config.get_fungi_assembly(fungi_name)
    with open(assembly_fasta, 'r') as f:
        scaffold_seq_dict = {desc: seq for desc, seq in fl.read_fasta(f)}

    intron_positions = []
    intron_sequences = []

    for gene in db.all_features(featuretype='gene'):
        intron_start = None
        print_next = False

        for i, f in enumerate(db.children(gene[gene_identifier][0], featuretype='exon')):
            scaffold, strand = f.chrom, f.strand
            start, end = f.start - 1, f.end

            # If first exon, just set the start/end of the following intron
            if i == 0 and strand == '+':
                intron_start = end  # If first exon and + strand, the exon's end is the intron start
                continue

            if i == 0 and strand == '-':
                # If strand - however, the exon order is reversed and we start with the last exon. This is the way
                # gff3 files are written. Hence the last exon's start is the intron's end
                intron_end = start
                continue

            if strand == '+':
                intron_end = start
            elif strand == '-':
                intron_start = end

            if intron_start > intron_end:
                print(f'Anomally at strand {strand}. Next introns of the gene are:')
                print_next = True
                intron_start = end
                continue

            # Extract intron
            scaffold_seq = scaffold_seq_dict[scaffold]
            intron_seq = Seq(scaffold_seq[intron_start: intron_end])

            if intron_end - intron_start != len(intron_seq):
                print('Valid coordinates but no intron')
                print([intron_end - intron_start, len(intron_seq)])
                intron_start = end
                continue

            if strand == "+":
                assert intron_end - intron_start == len(intron_seq)

                intron_positions.append([scaffold, strand, intron_start, intron_end])
                intron_sequences.append(
                    SeqRecord(id=' '.join([scaffold, strand, str(intron_start), str(intron_end)]),
                              seq=intron_seq)
                )
            elif strand == '-':
                intron_seq = Seq(scaffold_seq[intron_start: intron_end]).reverse_complement()

                # Coordinates relative to the minus strand
                intron_start_conv = len(scaffold_seq) - intron_end
                intron_end_conv = len(scaffold_seq) - intron_start

                # Verify the sequence and intron locations are correct
                intron_seq2 = Seq(scaffold_seq).reverse_complement()[intron_start_conv:intron_end_conv]
                assert intron_seq == intron_seq2
                assert intron_end_conv - intron_start_conv == len(intron_seq)

                intron_positions.append([scaffold, strand, intron_start_conv, intron_end_conv])
                intron_sequences.append(
                    SeqRecord(id=' '.join([scaffold, strand, str(intron_end_conv), str(intron_end_conv)]),
                              seq=intron_seq)
                )

            if print_next:
                print(intron_seq)

            # Move cursor to the next intron
            if strand == "+":
                intron_start = end
            elif strand == '-':
                intron_end = start

    write_output(fungi_name, intron_positions, intron_sequences, type='intron')


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
    filename, file_extension = os.path.splitext(gff_file)
    if file_extension == '.gff':
        extract_gff_introns(gff_file, fungi_name)
    else:
        assert file_extension == '.gff3'
        extract_gff3_introns(gff_file, fungi_name)


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

            if strand == '+' and validate:
                scaffold_seq = scaffold_seq_dict[scaffold]
                assert scaffold_seq[start: end] == exon_seq

                exon_positions.append([scaffold, strand, start, end])
                exon_sequences.append(
                    SeqRecord(id=' '.join([scaffold, strand, str(start), str(end)]),
                              seq=Seq(exon_seq))
                )

            if strand == '-' and validate:
                scaffold_seq = scaffold_seq_dict[scaffold]

                exon_start_conv = len(scaffold_seq) - end
                exon_end_conv = len(scaffold_seq) - start

                # Verify the sequence and intron locations are correct
                exon_seq2 = Seq(scaffold_seq).reverse_complement()[exon_start_conv:exon_end_conv]
                assert exon_seq2 == exon_seq

                exon_positions.append([scaffold, strand, exon_start_conv, exon_end_conv])
                exon_sequences.append(
                    SeqRecord(id=' '.join([scaffold, strand, str(exon_start_conv), str(exon_end_conv)]),
                              seq=Seq(exon_seq))
                )

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
