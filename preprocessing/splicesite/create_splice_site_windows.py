import csv
import logging
import os
import random
import sys
from pathlib import Path
from typing import Tuple, List

import numpy as np
from Bio import SeqIO

import config
from extract_tools import ACCEPTOR, DONOR

POSITIVE_LABEL = '+1'
NEGATIVE_LABEL = '-1'

logging.getLogger().setLevel(logging.INFO)
random.seed(42)


def main():
    """
    For a given fungi species, generate a CSV file with true&false splice site windows.
        - Each window will be @xx_lindow and @yy_rwindow long
    The script takes 2 files - e.g. true-donors and false-donors and mix them into the resulting CSV
        - Expect the true/false files to exist. Run "extract_false/true_splice_sites.py" first
        - Can specify the ratio between positive and negative examples by arguments @max_pos_samples, @max_neg_samples
    The file can be either for testing/training, which is specified by @test_train parameter
    :return: 2 files - donor-train/test.csv acceptor-train/test.csv
    """
    fungi_name = sys.argv[1]

    donor_lwindow, donor_rwindow = int(sys.argv[2]), int(sys.argv[3])
    acceptor_lwindow, acceptor_rwindow = int(sys.argv[4]), int(sys.argv[5])

    csv_target_folder = sys.argv[6]
    test_train = sys.argv[7]

    max_pos_samples = int(sys.argv[8])
    max_neg_samples = int(sys.argv[9])

    # fungi_name = 'Mycreb1'
    # csv_target_folder = '../data/'
    #
    # donor_lwindow, donor_rwindow = 200, 200
    # acceptor_lwindow, acceptor_rwindow = 200, 200
    #
    # test_train = 'train'
    # max_pos_samples = 100000
    # max_neg_samples = 140000

    assembly = config.get_fungi_assembly(fungi_name)
    if not Path(assembly).is_file():
        logging.warning(f'Assembly data for shroom {fungi_name} not found. Directory {assembly}')
        exit(1)

    suffix = 'false-intragenic'
    false_donor_file = f'{config.NEWSEQUENCES_LOC}/{fungi_name}/{fungi_name}-donor-{suffix}.fasta'
    false_acceptor_file = f'{config.NEWSEQUENCES_LOC}/{fungi_name}/{fungi_name}-acceptor-{suffix}.fasta'

    true_donor_file = f'{config.NEWSEQUENCES_LOC}/{fungi_name}/{fungi_name}-donor-true.fasta'
    true_acceptor_file = f'{config.NEWSEQUENCES_LOC}/{fungi_name}/{fungi_name}-acceptor-true.fasta'

    out_donor_csv, out_acceptor_csv = prepare_output(fungi_name, csv_target_folder, test_train)

    dwins_pos = get_positive_windows(true_donor_file, donor_lwindow, donor_rwindow, max_pos_samples)
    dwins_neg = get_negative_windows(false_donor_file, donor_lwindow, donor_rwindow, max_neg_samples)

    create_csv(out_donor_csv, dwins_pos, dwins_neg)

    awins_neg = get_negative_windows(false_acceptor_file, acceptor_lwindow, acceptor_rwindow, max_neg_samples)
    awins_pos = get_positive_windows(true_acceptor_file, acceptor_lwindow, acceptor_rwindow, max_pos_samples)

    create_csv(out_acceptor_csv, awins_pos, awins_neg)


def create_csv(target_csv: str, positive_windows: list, negative_windows: list):
    with open(target_csv, 'w') as f:
        w = csv.writer(f, delimiter=';')
        w.writerow(("sequence", "label"))

        for row in positive_windows:
            w.writerow(row)

        for row in negative_windows:
            w.writerow(row)

    logging.info(f'Written {len(positive_windows)}/{len(negative_windows)} positive/negative'
                 f'splice site windows to {target_csv}')


def get_positive_windows(
        true_donor_acceptor_file: str,
        lwindow: int,
        rwindow: int,
        max_samples: int
) -> List[Tuple]:
    splice_windows = get_splice_site_windows(true_donor_acceptor_file, lwindow, rwindow, max_samples)

    labeled_windows = zip(
        splice_windows,
        [POSITIVE_LABEL] * len(splice_windows)
    )

    return list(labeled_windows)


def get_negative_windows(
        false_donor_acceptor_file: str,
        lwindow: int,
        rwindow: int,
        max_samples: int
) -> List[Tuple]:
    splice_windows = get_splice_site_windows(false_donor_acceptor_file, lwindow, rwindow, max_samples)

    labeled_windows = zip(
        splice_windows,
        [NEGATIVE_LABEL] * len(splice_windows)
    )

    return list(labeled_windows)


def get_splice_site_windows(
        donor_acceptor_file: str,
        lwindow: int,
        rwindow: int,
        max_samples: int,
):
    with open(donor_acceptor_file, 'r') as f:
        # Read sequences from FASTA file
        all_windows = [str(sr.seq) for sr in SeqIO.parse(f, 'fasta')]
        random.shuffle(all_windows)

        # Take either max_samples or all available if there is not enough windows
        splice_windows = all_windows[0:min(max_samples, len(all_windows))]

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

        return splice_windows


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
