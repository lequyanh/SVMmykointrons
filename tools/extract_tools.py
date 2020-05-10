import glob
import logging
import random
import sys
from enum import Enum
from pathlib import Path
from typing import Optional, Tuple, List

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
    FALSE_SPLICESITES = 1
    FALSE_INTRAGEN_SPLICESITES = 2
    TRUE_SPLICESITES = 3


def main():
    # option = ExtractOptions.TRUE_SPLICESITES
    # fungi_name = 'Tripe1'
    # margin_size = 200
    # examples_limit = 25000

    option = ExtractOptions(int(sys.argv[1]))
    fungi_name = sys.argv[2]
    margin_size = int(sys.argv[3])
    examples_limit = int(sys.argv[4])

    if option == ExtractOptions.FALSE_SPLICESITES:
        donor_windows, acceptor_windows = generate_false_splicesites(fungi_name, margin_size, examples_limit)
        suffix = 'false'
    elif option == ExtractOptions.FALSE_INTRAGEN_SPLICESITES:
        donor_windows, acceptor_windows = generate_false_exon_splicesites(fungi_name, margin_size, examples_limit)
        suffix = 'false-intragenic'
    elif option == ExtractOptions.TRUE_SPLICESITES:
        donor_windows, acceptor_windows = generate_true_splice_sites(fungi_name, margin_size)
        suffix = 'true'
    else:
        print(f'No valid option {option} for splice site windows extraction. Exiting...')
        return

    with open(f'{config.NEWSEQUENCES_LOC}/{fungi_name}/{fungi_name}-donor-{suffix}.fasta', 'w') as f:
        SeqIO.write(donor_windows, f, 'fasta')

    with open(f'{config.NEWSEQUENCES_LOC}/{fungi_name}/{fungi_name}-acceptor-{suffix}.fasta', 'w') as f:
        SeqIO.write(acceptor_windows, f, 'fasta')


def generate_true_splice_sites(
        fungi_name: str,
        margin_size: int,
) -> Tuple[List, List]:
    introns_fasta = config.get_fungi_intron_fasta(fungi_name)
    donor_positions, acceptor_positions = true_donorac_positions(introns_fasta)

    assembly_fasta = config.get_fungi_assembly(fungi_name)
    with open(assembly_fasta, 'r') as assembly_f:
        true_donors = list()
        true_acceptors = list()

        # loop through the whole assembly
        for seq_rec in SeqIO.parse(assembly_f, 'fasta'):
            scaffold = seq_rec.description
            sequence = str(seq_rec.seq)

            for donor_pos in donor_positions.get(scaffold, []):
                window = extract_window(sequence, donor_pos, margin_size, scaffold)

                if window and str(window.seq[margin_size:margin_size + 2]) == DONOR:
                    true_donors.append(window)

            for acc_pos in acceptor_positions.get(scaffold, []):
                window = extract_window(sequence, acc_pos, margin_size, scaffold)

                if window and str(window.seq[margin_size:margin_size + 2]) == ACCEPTOR:
                    true_acceptors.append(window)

        print(f'\t True splice site windows for fungi {fungi_name} written: '
              f'{len(true_donors)} donor, {len(true_acceptors)} acceptor windows')

        return true_donors, true_acceptors


def generate_false_splicesites(
        fungi_name: str,
        margin_size: int,
        examples_limit: int
) -> Tuple[List, List]:
    """
    Extract false splice site windows anywhere in the genome of the given fungi. Scaffold, from which the windows will
    be extracted are shuffled and sampled.
    The windows will be saved to files at the new_sequences folder.

    :param examples_limit: Limit on the number of false windows to extract (for each donor/acceptor windows)
    :param margin_size: Width of the window
    :param fungi_name: Fungi name to extract false intra-genic splice site windows from
    :return: Two list of donor and acceptor windows
    """
    assembly_fasta = config.get_fungi_assembly(fungi_name)
    with open(assembly_fasta, 'r') as assembly_f:
        scaffold_records = list(SeqIO.parse(assembly_f, 'fasta'))  # type: List[SeqRecord]
        random.shuffle(scaffold_records)  # shuffle scaffolds for uniform sampling

    introns_fasta = config.get_fungi_intron_fasta(fungi_name)
    donor_positions, acceptor_positions = true_donorac_positions(introns_fasta)

    false_donors = list()
    false_acceptors = list()

    def yield_false_windows():
        print(f'False splice site windows for fungi {fungi_name} written:'
              f'{len(false_donors)} donor, {len(false_acceptors)} acceptor windows / {examples_limit} limit')

        return false_donors, false_acceptors

    for seq_record in scaffold_records:
        sequence = str(seq_record.seq)
        scaffold = seq_record.description

        for position, dimer in fl.dimers(sequence):
            if len(false_donors) >= examples_limit:
                return yield_false_windows()

            if dimer == DONOR and position not in donor_positions.get(scaffold, []):
                window = extract_window(sequence, position, margin_size, scaffold)
                false_donors = false_donors + [window] if window else false_donors

            elif dimer == ACCEPTOR and position not in acceptor_positions.get(scaffold, []):
                window = extract_window(sequence, position, margin_size, scaffold)
                false_acceptors = false_acceptors + [window] if window else false_acceptors
            else:
                logging.info(f'Dimer is a splice site, skip it')
                continue

    return yield_false_windows()


def generate_false_exon_splicesites(
        fungi_name: str,
        margin_size: int,
        examples_limit: int
) -> Tuple[List, List]:
    """
    Extract false splice site windows from exon regions of the given fungi. The regions will be sampled
    The windows will be saved to files at the new_sequences folder.
    NOTE: This method requires exon positions file to exist. If not, run the @extract_exons() function first
          on the selected fungi

    :param examples_limit: Limit on the number of false windows to extract (for each donor/acceptor windows)
    :param margin_size: Width of the window
    :param fungi_name: Fungi name to extract false intra-genic splice site windows from
    :return: Two list of donor and acceptor windows
    """
    assembly_fasta = config.get_fungi_assembly(fungi_name)
    with open(assembly_fasta, 'r') as assembly_file:
        scaffolds_input = {desc: seq for desc, seq in fl.read_fasta(assembly_file)}

    exons_file = config.get_fungi_exons_fasta(fungi_name)
    if not Path(exons_file).is_file():
        extract_exons(fungi_name)

    with open(exons_file, 'r') as ef:
        exons = list(SeqIO.parse(ef, 'fasta'))  # type: List[SeqRecord]
        random.shuffle(exons)

    intragenic_false_donors = list()
    intragenic_false_acceptors = list()

    def yield_false_windows():
        print(f'False intragenic splice site windows for fungi {fungi_name} written:'
              f'{len(intragenic_false_donors)} donor, {len(intragenic_false_acceptors)} acceptor windows '
              f'/ {examples_limit} limit')

        return intragenic_false_donors, intragenic_false_acceptors

    def window_coords(dimer_relative_position, exon_start):
        # First, compute relative start and end for the window within the exon
        win_start = dimer_relative_position - margin_size - 1
        win_end = dimer_relative_position + margin_size + 1

        # Adjust the window coordinates by the position of the exon
        return exon_start + win_start, exon_start + win_end

    def extract_exon_window(dimer_exon_position, exon_start):
        win_start, win_end = window_coords(dimer_exon_position, exon_start)

        window = scaffolds_input[scaffold][win_start:win_end]
        if len(window) != 2 * margin_size + 2:
            logging.warning("Windows too short, skipping...")
            return None

        win_mid = int(len(window) / 2)
        assert window[win_mid - 1: win_mid + 1] in {DONOR, ACCEPTOR}

        return SeqRecord(
            id=' '.join([scaffold, strand, str(win_start), str(win_end)]),
            seq=Seq(window)
        )

    for exon_seqrec in exons:  # type: SeqRecord
        desc = exon_seqrec.description.split(' ')

        scaffold, strand, start, end = desc[0], desc[1], int(desc[2]), int(desc[3])

        # Compare the exon we get from the original fasta with the sequence in the exon file.
        # We can't use directly sequences from the exon file, as we might need the surrounding nucleotides
        exon_seq = scaffolds_input[scaffold][start:end]
        assert exon_seq == str(exon_seqrec.seq)

        for position, dimer in fl.dimers(exon_seq):
            if len(intragenic_false_donors) >= examples_limit:
                return yield_false_windows()

            if dimer == DONOR or dimer == ACCEPTOR:
                false_splicesite = extract_exon_window(position, start)

                if dimer == DONOR and false_splicesite:
                    intragenic_false_donors.append(false_splicesite)
                elif dimer == ACCEPTOR and false_splicesite:
                    intragenic_false_acceptors.append(false_splicesite)

    return yield_false_windows()


def extract_exons(
        fungi_name: str
) -> None:
    """
    Extracts exon positions to a CSV file and exon sequences to FASTA file. Run by:

        fungi_names = f'../basidiomycota.txt'
        with open(fungi_names, 'r') as f:
            names = f.readlines()
            names = map(lambda name: name.rstrip(), names)
            list(map(extract_exons, names))

    :param fungi_name: Fungi to extract exons for
    :return: None
    """
    gff = glob.glob(f'{config.GFFS_LOC}/{fungi_name}_GeneCatalog_genes_*.gff')
    if len(gff) == 0:
        print(f'No GFF file for fungi {fungi_name}')
        return
    if len(gff) > 1:
        print(f'More than one GFF file for fungi {fungi_name}. Taking first.')

    gff_file = gff[0]
    print(f'Extracting exons for fungi {fungi_name}')

    gfflib.extract_exon_positions(gff_file, fungi_name, validate=True)


def true_donorac_positions(
        introns_fasta: str
) -> Tuple[dict, dict]:
    """
    Extracts positions of true donors and acceptors (only canonic GT and AG ones).
    Record number of non-canonic splice sites.

    :param introns_fasta: FASTA file with intron sequences and their positions
    :return: Two dictionaries. One with donor positions and one with acceptor positions.
                Keys are scaffold IDs, values list with coordinates
    """
    with open(introns_fasta, 'r') as intron_fasta_f:
        donor_positions = dict()
        acceptor_positions = dict()

        non_canonic_acc = 0
        non_canonic_don = 0

        introns_seqrecords = list(SeqIO.parse(intron_fasta_f, 'fasta'))
        for isr in introns_seqrecords:
            scaff_id, strand, start, end = isr.description.split(' ')

            if strand == '+':
                if str(isr.seq[:2]) != DONOR:
                    non_canonic_don += 1
                else:
                    donor_start = int(start)

                    don_pos = donor_positions.get(scaff_id, [])
                    don_pos.append(donor_start)
                    donor_positions[scaff_id] = don_pos

                if str(isr.seq[-2:]) != ACCEPTOR:
                    non_canonic_acc += 1
                else:
                    acceptor_start = int(end) - 1
                    acc_pos = acceptor_positions.get(scaff_id, [])
                    acc_pos.append(acceptor_start)
                    acceptor_positions[scaff_id] = acc_pos

        print(f'Parsed intron file {introns_fasta} '
              f'Non-canonic dimers: {non_canonic_don} donors and {non_canonic_acc} acceptors.')

    return donor_positions, acceptor_positions


def extract_window(
        sequence: str,
        position: int,
        margin_size: int,
        scaffold: str
) -> Optional[SeqRecord]:
    # compute start and end for the window
    win_start = position - margin_size - 1
    win_end = position + margin_size + 1

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
        id=' '.join([scaffold, '+', str(win_start), str(win_end)]),
        seq=Seq(window)
    )


if __name__ == '__main__':
    main()
