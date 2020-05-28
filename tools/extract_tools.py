import glob
import logging
import random
import sys
from enum import Enum
from typing import List
from typing import Optional
from typing import Tuple

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import config
import fastalib as fl
import gfflib

DONOR = 'GT'
ACCEPTOR = 'AG'

logging.basicConfig(
    level=logging.INFO,
    filename='extract-tools.log',
    filemode='w'
)
random.seed(42)


class ExtractOptions(Enum):
    FALSE_INTRAGEN_SPLICESITES = 1
    TRUE_SPLICESITES = 2


def main():
    # option = ExtractOptions.TRUE_SPLICESITES
    # fungi_name = 'Tripe1'
    # margin_size = 200
    # examples_limit = 25000

    option = ExtractOptions(int(sys.argv[1]))
    fungi_name = sys.argv[2]
    strand = sys.argv[3]
    margin_size = int(sys.argv[4])
    examples_limit = int(sys.argv[5]) if len(sys.argv) == 6 else None

    if option == ExtractOptions.FALSE_INTRAGEN_SPLICESITES:
        donor_windows, acceptor_windows = \
            generate_false_exon_splicesites(fungi_name, margin_size, examples_limit, strand)

        fasta_donor_wins = config.get_fungi_false_wins_fasta(fungi_name, strand, 'donor')
        fasta_acc_wins = config.get_fungi_false_wins_fasta(fungi_name, strand, 'acceptor')
    elif option == ExtractOptions.TRUE_SPLICESITES:
        donor_windows, acceptor_windows = \
            generate_true_splice_sites(fungi_name, margin_size, strand)

        fasta_donor_wins = config.get_fungi_true_wins_fasta(fungi_name, strand, 'donor')
        fasta_acc_wins = config.get_fungi_true_wins_fasta(fungi_name, strand, 'acceptor')
    else:
        print(f'No valid option {option} for splice site windows extraction. Exiting...')
        return

    with open(fasta_donor_wins, 'w') as f:
        SeqIO.write(donor_windows, f, 'fasta')

    with open(fasta_acc_wins, 'w') as f:
        SeqIO.write(acceptor_windows, f, 'fasta')


def generate_true_splice_sites(
        fungi_name: str,
        margin_size: int,
        from_strand: str
) -> Tuple[List, List]:
    introns_fasta = config.get_fungi_intron_fasta(fungi_name)
    assembly_fasta = config.get_fungi_assembly(fungi_name)

    # Get donor and acceptor positions from the given strand
    donor_positions, acceptor_positions = true_donorac_positions(introns_fasta, from_strand)

    with open(assembly_fasta, 'r') as assembly_f:
        donor_wins = pd.Series()
        acceptor_wins = pd.Series()

        # loop through the whole assembly to extract windows on the given strand
        for scaffold_sr in SeqIO.parse(assembly_f, 'fasta'):
            scaffold = scaffold_sr.description
            # If negative strand, reverse complement the scaffold
            scaffold_seq = str(scaffold_sr.seq if from_strand == '+' else scaffold_sr.seq.reverse_complement())

            # Extract donor and acceptor windows from the current scaffold
            temp_don_wins = donor_positions.loc[donor_positions['scaffold'] == scaffold, 'position'] \
                .apply(lambda position: extract_window(scaffold_seq, position, margin_size, scaffold, from_strand))

            temp_acc_wins = acceptor_positions.loc[acceptor_positions['scaffold'] == scaffold, 'position'] \
                .apply(lambda position: extract_window(scaffold_seq, position, margin_size, scaffold, from_strand))

            donor_wins = pd.concat([donor_wins, temp_don_wins], ignore_index=True)
            acceptor_wins = pd.concat([acceptor_wins, temp_acc_wins], ignore_index=True)

        donor_wins.dropna(inplace=True)
        acceptor_wins.dropna(inplace=True)

        print(f'\t True splice site windows for fungi {fungi_name} written (windows/positions total): '
              f'{donor_wins.shape[0]}/{donor_positions.shape[0]} donor, '
              f'{acceptor_wins.shape[0]}/{acceptor_positions.shape[0]} acceptor windows ')

        return donor_wins.to_list(), acceptor_wins.to_list()


def generate_false_exon_splicesites(
        fungi_name: str,
        margin_size: int,
        examples_limit: int,
        from_strand: str
) -> Tuple[List, List]:
    """
    Extract false splice site windows from exon regions of the given fungi. The regions will be sampled
    The windows will be saved to files at the new_sequences folder.
    NOTE: This method requires exon positions file to exist. If not, run the @extract_exons() function first
          on the selected fungi

    :param from_strand: Chosen strand to generate false splice sites from
    :param examples_limit: Limit on the number of false windows to extract (for each donor/acceptor windows)
    :param margin_size: Width of the window
    :param fungi_name: Fungi name to extract false intra-genic splice site windows from
    :return: Two list of donor and acceptor windows
    """
    assembly_fasta = config.get_fungi_assembly(fungi_name)
    with open(assembly_fasta, 'r') as assembly_file:
        scaffolds_input = {desc: seq for desc, seq in fl.read_fasta(assembly_file)}

    exons_file = config.get_fungi_exons_fasta(fungi_name)
    with open(exons_file, 'r') as ef:
        exons = list(SeqIO.parse(ef, 'fasta'))  # type: List[SeqRecord]

    random.shuffle(exons)

    intragenic_false_donors = list()
    intragenic_false_acceptors = list()

    def yield_false_windows():
        print(f'False intragenic splice site windows for fungi {fungi_name} written: '
              f'{len(intragenic_false_donors)} donor, {len(intragenic_false_acceptors)} acceptor '
              f'(limit {examples_limit})')

        return intragenic_false_donors, intragenic_false_acceptors

    for exon_seqrec in exons:  # type: SeqRecord
        desc = exon_seqrec.description.split(' ')

        scaffold, strand, start, end = \
            desc[0], desc[1], int(desc[2]), int(desc[3])

        if strand != from_strand:
            continue

        if strand == '+':
            scaffold_seq = scaffolds_input[scaffold]
        else:
            scaffold_seq = str(Seq(scaffolds_input[scaffold]).reverse_complement())
        exon_seq = scaffold_seq[start:end]

        # Compare the exon we get from the original fasta with the sequence in the exon file.
        assert exon_seq == str(exon_seqrec.seq)

        # Extract windows from GT/AG dimers inside exons
        for position, dimer in fl.dimers(exon_seq):
            if len(intragenic_false_donors) >= examples_limit:
                return yield_false_windows()

            if dimer == DONOR or dimer == ACCEPTOR:
                dimer_pos = start + position - 1  # The dimers() function adds 1 to the index. Correct for it
                false_splicesite = extract_window(scaffold_seq, dimer_pos, margin_size, scaffold, strand)

                if dimer == DONOR and false_splicesite:
                    intragenic_false_donors.append(false_splicesite)
                elif dimer == ACCEPTOR and false_splicesite:
                    intragenic_false_acceptors.append(false_splicesite)

    return yield_false_windows()


def extract_intron_exons(
        fungi_name: str,
        exon_flag: bool = True
) -> None:
    """
    Parses GFF file of a fungi (file location specified in config)
    Extracts exon/intron positions to a CSV file and their sequences to a FASTA file. Run by:

        fungi_names = f'../ascomycota.txt'
        exons = False
        with open(fungi_names, 'r') as f:
            names = f.readlines()
            names = map(lambda name: name.rstrip(), names)
            list(map(lambda name: extract_intron_exons(name, exons), names))

    :param exon_flag: If the flag is True, exons are extracted. If false, introns are extracted
    :param fungi_name: Fungi to extract exons/introns for
    :return: None
    """
    gff = glob.glob(f'{config.GFFS_LOC}/{fungi_name}_GeneCatalog_genes_*.gff')
    if len(gff) == 0:
        print(f'No GFF file for fungi {fungi_name}')
        return
    if len(gff) > 1:
        print(f'More than one GFF file for fungi {fungi_name}. Taking first.')

    gff_file = gff[0]
    print(f'Extracting {"exons" if exon_flag else "introns"} for fungi {fungi_name}')

    if exon_flag:
        gfflib.extract_exon_positions(gff_file, fungi_name, validate=True)
    else:
        gfflib.extract_intron_positions(gff_file, fungi_name)


def true_donorac_positions(
        introns_fasta: str,
        strand: str
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Reads fungi's intron.fasta file and extracts positions of true donors and acceptors (only canonic GT and AG ones).
    Records number of non-canonic splice sites.

    :param strand: DNA strand (+/-)
    :param introns_fasta: FASTA file with intron sequences and their positions
    :return: Two dictionaries. One with donor positions and one with acceptor positions.
                Keys are scaffold IDs, values list with coordinates
    """
    with open(introns_fasta, 'r') as intron_fasta_f:
        introns_seqrecords = list(SeqIO.parse(intron_fasta_f, 'fasta'))

    donor_data = list()
    acceptor_data = list()

    non_canonic_acc, non_canonic_don = 0, 0
    for isr in introns_seqrecords:
        scaff_id, intron_strand, start, end = isr.description.split(' ')[0:4]

        # Skip if the intron is not from strand we are interested
        if intron_strand != strand:
            continue

        if str(isr.seq[:2]) != DONOR:
            non_canonic_don += 1
        else:
            donor_start = int(start)
            donor_data.append([scaff_id, donor_start])

        if str(isr.seq[-2:]) != ACCEPTOR:
            non_canonic_acc += 1
        else:
            acceptor_start = int(end) - 2
            acceptor_data.append([scaff_id, acceptor_start])

    columns = ['scaffold', 'position']
    donor_positions = pd.DataFrame(donor_data, columns=columns)
    acceptor_positions = pd.DataFrame(acceptor_data, columns=columns)

    print(f'Getting ({strand}) strand splice-site positions. '
          f'Non-canonic splice-site dimers: {non_canonic_don} donors and {non_canonic_acc} acceptors.')

    return donor_positions, acceptor_positions


def extract_window(
        sequence: str,
        position: int,
        margin_size: int,
        scaffold: str,
        strand: str
) -> Optional[SeqRecord]:
    # compute start and end for the window
    win_start = position - margin_size
    win_end = position + margin_size + 2

    if win_start < 0:
        return None

    # extract the window
    window = sequence[win_start:win_end]

    dimer = window[margin_size:margin_size + 2]
    if dimer not in {DONOR, ACCEPTOR}:
        return None

    if len(window) < 2 * margin_size + 2:
        logging.warning("Windows too short, skipping...")
        return None

    return SeqRecord(
        id=' '.join([scaffold, strand, str(win_start), str(win_end)]),
        seq=Seq(window)
    )


if __name__ == '__main__':
    main()
