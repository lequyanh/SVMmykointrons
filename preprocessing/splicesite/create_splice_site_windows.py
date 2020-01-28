import csv
import logging
import os
import sys
from pathlib import Path

import numpy as np
from Bio import SeqIO

from extract_tools import ACCEPTOR, DONOR
from fastalib import read_fasta

POSITIVE_LABEL = '+1'
NEGATIVE_LABEL = '-1'


def main():
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

    logging.getLogger().setLevel(logging.INFO)

    assembly = f'{base_loc}/data/Assembly/{shroom_name}_AssemblyScaffolds.fasta'
    if not Path(assembly).is_file():
        logging.warning(f'Assembly data for shroom {shroom_name} not found. Directory {assembly}')
        exit(1)

    false_donor_file = f'{base_loc}/new-sequences/{shroom_name}/{shroom_name}-donor-false.fasta'
    false_acceptor_file = f'{base_loc}/new-sequences/{shroom_name}/{shroom_name}-acceptor-false.fasta'

    out_donor_csv, out_acceptor_csv = prepare_output(shroom_name, csv_target_folder, test_train)

    append_false_splice_site_windows(false_donor_file, out_donor_csv, donor_lwindow, donor_rwindow)
    append_false_splice_site_windows(false_acceptor_file, out_acceptor_csv, acceptor_lwindow, acceptor_rwindow)

    append_true_splice_site_windows(true_acceptor_file, out_acceptor_csv, acceptor_lwindow, acceptor_rwindow)
    append_true_splice_site_windows(true_donor_file, out_donor_csv, donor_lwindow, donor_rwindow)


def append_true_splice_site_windows(
        true_donor_acceptor_file: str,
        target_csv: str,
        lwindow: int,
        rwindow: int,
):
    append_splice_site_windows_(true_donor_acceptor_file, target_csv, lwindow, rwindow, POSITIVE_LABEL)


def append_false_splice_site_windows(
        false_donor_acceptor_file: str,
        target_csv: str,
        lwindow: int,
        rwindow: int,
):
    append_splice_site_windows_(false_donor_acceptor_file, target_csv, lwindow, rwindow, NEGATIVE_LABEL)


def append_splice_site_windows_(
        donor_acceptor_file: str,
        target_csv: str,
        lwindow: int,
        rwindow: int,
        label: str,
):
    with open(donor_acceptor_file, 'r') as f:
        # Read sequences from FASTA file
        splice_windows = [str(sr.seq) for sr in SeqIO.parse(f, 'fasta')]

        window_mid = int(0.5 * (len(splice_windows[1]) - 2))  # -2 for GT/AG dimers

        # Adjust the range of sequences
        adjusted_range = slice(window_mid - lwindow, window_mid + rwindow + 2)
        splice_windows = list(
            map(lambda window: window.rstrip()[adjusted_range],
                splice_windows)
        )

        # Validate the windows are correct
        dimers_check = map(
            lambda false_splice_window: false_splice_window[lwindow:lwindow + 2] in {DONOR, ACCEPTOR},
            splice_windows
        )
        assert np.all(dimers_check)

        # Create CSV rows
        labeled_windows = zip(
            splice_windows,
            [label] * len(splice_windows)
        )

        with open(target_csv, 'a') as f:
            w = csv.writer(f, delimiter=';')

            for row in labeled_windows:
                w.writerow(row)

        logging.info(f'Appended {len(splice_windows)} splice site windows to {target_csv}')


def prepare_output(shroom_name: str, target_folder: str, test_train: str):
    out_donor_dir = f'{target_folder}/{test_train}/donor'
    out_acceptor_dir = f'{target_folder}/{test_train}/acceptor'

    if not os.path.isdir(out_donor_dir):
        os.makedirs(out_donor_dir)
    if not os.path.isdir(out_acceptor_dir):
        os.makedirs(out_acceptor_dir)

    out_donor_csv = f'{out_donor_dir}/{shroom_name}-donor-windows.csv'
    out_acceptor_csv = f'{out_acceptor_dir}/{shroom_name}-acceptor-windows.csv'

    if os.path.isfile(out_donor_csv) and os.path.isfile(out_acceptor_csv):
        logging.warning(f'Files {out_donor_csv} and {out_acceptor_csv} already exist')
        exit(0)

    logging.info(f'Splice site training/testing CSV will be saved to {out_donor_dir} and {out_acceptor_dir}')
    return out_donor_csv, out_acceptor_csv


if __name__ == "__main__":
    main()
