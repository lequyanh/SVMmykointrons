import pandas as pd
import numpy as np
import sys

# acceptor_file = sys.argv[1]
# donor_results_file = sys.argv[2]
# intron_min_len, intron_max_len = sys.argv[3:4]

acceptor_file = "splice-site-acceptor-dataset.csv"
donor_results_file = "splice-site-donor-result.csv"
intron_min_len, intron_max_len = 40, 100

acc_cand_df = pd.read_csv(acceptor_file, sep=';')
donor_res_df = pd.read_csv(donor_results_file, sep=';')


def is_acc_pairable(don_pos, acc_pos):
    return don_pos + 40 < acc_pos < don_pos + 100


def next_acceptor(acc_positions, acc_idx):
    acc_idx += 1
    acc_pos = acc_positions[acc_idx] if acc_idx < len(acc_positions) else None

    return acc_idx, acc_pos


def get_applicable_acceptors(acc_positions, donor_positions):
    acc_idx, don_idx = 0, 0
    acc_pos = acc_positions[0]

    inrange_accs = np.empty(len(acc_positions))

    while don_idx < len(donor_positions) and acc_idx < len(acc_positions):

        don_pos = donor_positions['position'][don_idx]
        if is_acc_pairable(don_pos, acc_pos):
            while is_acc_pairable(don_pos, acc_pos) and acc_idx < len(acc_positions):
                # If acceptor candidate in valid range, flag as OK and shift to next one
                inrange_accs[acc_idx] = True
                acc_idx, acc_pos = next_acceptor(acc_positions, acc_idx)

            # Now the acceptor is too far. We can move to the next donor
            # .....d.....|..........|..a..
            don_idx += 1

        # Not pairable cases
        elif don_pos < acc_pos and not is_acc_pairable(don_pos, acc_pos):
            # case 1: donor preceding acceptor (which is OK), but distance too short
            # .....d..a..|..........|.....
            # Label other too close acceptor candidates as false
            while don_pos < acc_pos and not is_acc_pairable(don_pos, acc_pos) and acc_idx < len(acc_positions):
                inrange_accs[acc_idx] = False
                acc_idx, acc_pos = next_acceptor(acc_positions, acc_idx)
        elif don_pos > acc_pos:
            # case 2: donor is ahead of acceptor (should be vice versa)
            # a........d....|..........|.....
            # Label all acceptor candidates between the current donor and previous acceptor as false
            while don_pos > acc_pos and acc_idx < len(acc_positions):
                inrange_accs[acc_idx] = False
                acc_idx, acc_pos = next_acceptor(acc_positions, acc_idx)
        else:
            raise ValueError


for scaffold, donor_positions_df in donor_res_df.groupby(by='scaffold'):
    acc_positions = acc_cand_df[acc_cand_df['scaffold'] == scaffold]['position'].values
    donor_positions = donor_positions_df['position'].values

    applicable_accs = get_applicable_acceptors(acc_positions, donor_positions)

    print(f'Number of acceptor candidates / survived candidates: {len(acc_positions) / sum(inrange_accs)}')
