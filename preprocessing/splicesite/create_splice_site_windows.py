import csv
import logging
import os
import sys
from pathlib import Path
from typing import List, Iterable

import numpy as np

from fastalib import complementary, read_fasta

# ==== Use standalone to extract train data for a single shroom, or in a batch mode using a script ====
# >bash create_train_test_csvs.sh $SPLICE_TRAIN_NAMES ~/Desktop/mykointrons-data 150 150 ../data train

base_loc = sys.argv[1]
shroom_name = sys.argv[2]

donor_lwindow, donor_rwindow = int(sys.argv[3]), int(sys.argv[4])
acceptor_lwindow, acceptor_rwindow = int(sys.argv[5]), int(sys.argv[6])

csv_target_folder = sys.argv[7]
test_train = sys.argv[8]

# base_loc = '/home/anhvu/Desktop/mykointrons-data'
# shroom_name = 'Ramac1'
# csv_target_folder = '../data/'
#
# donor_lwindow, donor_rwindow = 150, 150
# acceptor_lwindow, acceptor_rwindow = 150, 150
#
# test_train = 'train'

assembly = f'{base_loc}/data/Assembly/{shroom_name}_AssemblyScaffolds.fasta'
introns_locs = f'{base_loc}/new-sequences/{shroom_name}/{shroom_name}-introns.fasta'

false_donor_file = f'{base_loc}/new-sequences/{shroom_name}/{shroom_name}-donor-false.fasta'
false_acceptor_file = f'{base_loc}/new-sequences/{shroom_name}/{shroom_name}-acceptor-false.fasta'

donor_dir = f'{csv_target_folder}/{test_train}/donor'
acceptor_dir = f'{csv_target_folder}/{test_train}/acceptor'

if not os.path.isdir(donor_dir):
    os.makedirs(donor_dir)
if not os.path.isdir(acceptor_dir):
    os.makedirs(acceptor_dir)

donor_csv = f'{donor_dir}/{shroom_name}-donor-windows.csv'
acceptor_csv = f'{acceptor_dir}/{shroom_name}-acceptor-windows.csv'
logging.info(f'Splice site training/testing CSV will be saved to {donor_dir} and {acceptor_dir}')

logging.getLogger().setLevel(logging.INFO)

POSITIVE_LABEL = '+1'
NEGATIVE_LABEL = '-1'

DONOR_DIMER = 'GT'
ACCEPTOR_DIMER = 'AG'


def get_donor_acceptor_windows_from_intron_locations(position_infos, scaffold_dict):
    donor_windows = []
    acceptor_windows = []

    def reverse_complement(s: str):
        return complementary(s)[::-1]

    for pos_info in position_infos:
        scaffold, strand, start, end, intron = pos_info[0], pos_info[1], int(pos_info[2]), int(pos_info[3]), pos_info[4]

        seq = scaffold_dict[scaffold[1:]]  # remove the > at the beginning

        if strand == '+':
            donor_wind_interval = slice(start - donor_lwindow - 1, start + donor_rwindow + 1)
            acc_wind_interval = slice(end - acceptor_lwindow - 2, end + acceptor_rwindow)

            acceptor_window = seq[acc_wind_interval]
            donor_window = seq[donor_wind_interval]
        else:
            donor_wind_interval = slice(end - donor_rwindow - 2, end + donor_lwindow)
            acc_wind_interval = slice(start - acceptor_rwindow - 1, start + acceptor_lwindow + 1)

            donor_window = reverse_complement(seq[donor_wind_interval])
            acceptor_window = reverse_complement(seq[acc_wind_interval])

        if len(donor_window) < donor_lwindow + donor_rwindow + 2 or \
                len(acceptor_window) < acceptor_lwindow + acceptor_rwindow + 2:
            # Skipping window due to premature scaffold end (windows not long enough)
            continue

        donor_dimer = donor_window[donor_lwindow: donor_lwindow + 2]
        acceptor_dimer = acceptor_window[acceptor_lwindow: acceptor_lwindow + 2]

        try:
            assert donor_dimer == intron[0:2]
            assert acceptor_dimer == intron[-2:]
        except AssertionError:
            if not (donor_dimer.isupper() and acceptor_dimer.isupper()):
                continue
            if len(intron) < 10:
                continue
            print(f'Error on parsing assembly of {shroom_name}')
            exit(1)

        # We know the dimers are OK thanks to asserts above
        if donor_dimer != DONOR_DIMER or acceptor_dimer != ACCEPTOR_DIMER:
            logging.debug(f'Unusual donor-acceptor dimer of intron: \n {intron}')
            continue

        donor_windows.append(donor_window)
        acceptor_windows.append(acceptor_window)

    return donor_windows, acceptor_windows


def get_introns_position_info(introns_locs_fasta: str) -> Iterable[List]:
    # Fasta file with intron locations and sequences
    with open(introns_locs_fasta, 'r') as introns_locs_f:
        lines = introns_locs_f.readlines()

        # Scaffold, strand, start, end
        position_infos = lines[::2]
        # Introns sequence
        intron_seqs = lines[1::2]

        def parse_rec(rec):
            info_fields = rec[0].rstrip().split(' ')
            intron_seq = rec[1].rstrip()

            return info_fields + [intron_seq]

        # Join parsed location info with the sequence itself
        position_infos = map(
            parse_rec,
            zip(position_infos, intron_seqs)
        )

        return position_infos


def write_true_splice_site_windows(
        train_windows: List[str],
        splice_windows_train_csv: str
):
    train_data_positive = zip(
        train_windows,
        [POSITIVE_LABEL] * len(train_windows)
    )

    with open(splice_windows_train_csv, 'w') as f:
        w = csv.writer(f, delimiter=';')
        w.writerow(("sequence", "label"))

        for row in train_data_positive:
            w.writerow(row)

    logging.info(f'Wrote {len(train_windows)} positive splice site windows to {splice_windows_train_csv}')


def append_true_splice_site_windows(
        true_donor_acceptor_file: str,
        target_csv: str,
        lwindow: int,
        rwindow: int,
):
    append_examples(true_donor_acceptor_file, target_csv, lwindow, rwindow, POSITIVE_LABEL)


def append_false_splice_site_windows(
        false_donor_acceptor_file: str,
        target_csv: str,
        lwindow: int,
        rwindow: int,
):
    append_examples(false_donor_acceptor_file, target_csv, lwindow, rwindow, NEGATIVE_LABEL)


def append_examples(
        donor_acceptor_file: str,
        target_csv: str,
        lwindow: int,
        rwindow: int,
        label: str,
):
    with open(donor_acceptor_file, 'r') as f:
        # Read sequences from FASTA file
        false_splice_windows = f.readlines()[1::2]

        # TODO Ramac1 is empty and failing here
        window_mid = int(0.5 * (len(false_splice_windows[1]) - 2))  # -2 for GT/AG dimers

        # Adjust the sequences so they have the same length as positive examples
        adjusted_range = slice(window_mid - lwindow, window_mid + rwindow + 2)
        false_splice_windows = map(
            lambda false_splice_window: false_splice_window.rstrip()[adjusted_range],
            false_splice_windows
        )
        false_splice_windows = list(false_splice_windows)

        splice_dimers_check = map(
            lambda false_splice_window: false_splice_window[lwindow:lwindow + 2] in {DONOR_DIMER, ACCEPTOR_DIMER},
            false_splice_windows
        )
        assert np.all(splice_dimers_check)

        # Create CSV rows
        data_negative = zip(
            false_splice_windows,
            [label] * len(false_splice_windows)
        )

        with open(target_csv, 'a') as f:
            w = csv.writer(f, delimiter=';')

            for row in data_negative:
                w.writerow(row)

        logging.info(f'Appended {len(false_splice_windows)} false splice site windows to {target_csv}')


if not Path(assembly).is_file():
    logging.warning(f'Assembly data for shroom {shroom_name} not found. Directory {assembly}')
    exit(1)

if os.path.isfile(donor_csv) and os.path.isfile(acceptor_csv):
    logging.warning(f'Files {donor_csv} and {acceptor_csv} already exist')
    exit(0)

scaffold_dict = {}
with open(assembly, 'r') as assembly_f:
    for desc, sequence in read_fasta(assembly_f):
        scaffold_dict[desc] = sequence

intron_positions = get_introns_position_info(introns_locs)

d_windows, a_windows = get_donor_acceptor_windows_from_intron_locations(intron_positions, scaffold_dict)

write_true_splice_site_windows(d_windows, donor_csv)
append_false_splice_site_windows(false_donor_file, donor_csv, donor_lwindow, donor_rwindow)

write_true_splice_site_windows(a_windows, acceptor_csv)
append_false_splice_site_windows(false_acceptor_file, acceptor_csv, acceptor_lwindow, acceptor_rwindow)

# append_positive_train_examples(true_acceptor_file, csv_acceptor_train, acceptor_lwindow, acceptor_rwindow)
# append_positive_train_examples(true_donor_file, csv_donors_train, donor_lwindow, donor_rwindow)


# tar -zcvf Agahy1-donor-train.csv.tar.gz ../data/splice-sites-train/donor/Agahy1-donor-dataset-splice-site-train.csv
# scp Agahy1-donor-train.csv.tar.gz lequyanh@skirit.metacentrum.cz:/storage/praha1/home/lequyanh
# qsub -l walltime=24:0:0 -l select=1:ncpus=10:mem=4gb:scratch_local=2gb -v DEGREE=15,LWINDOW=60,RWINDOW=70,C=1,DATAFILE=Armga1-donor-train.csv.tar.gz,CPU=10 train-splice-site.sh
