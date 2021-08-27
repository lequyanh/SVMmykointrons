import pandas as pd
import numpy as np
import sys

acceptor_file = sys.argv[1]
donor_results_file = sys.argv[2]
INTRON_MIN_LEN, INTRON_MAX_LEN = int(sys.argv[3]), int(sys.argv[4])

acc_cand_df = pd.read_csv(acceptor_file, sep=';')
donor_res_df = pd.read_csv(donor_results_file, sep=';')


def main():
    print('Removing acceptor candidates with no donor match')

    # Dataframe that will hold the reduced set of acceptor candidates
    new_acc_cand_df = pd.DataFrame(columns=['scaffold', 'position', 'sequence'])

    # For each scaffold, remove all infeasible acceptors (those that cannot have a matching donor)
    for scaffold, donor_positions_df in donor_res_df.groupby(by='scaffold'):
        acceptor_positions_df = acc_cand_df[acc_cand_df['scaffold'] == scaffold]

        donor_positions = donor_positions_df['position'].values
        acc_positions = acceptor_positions_df['position'].values

        if acc_positions.shape[0] == 0:
            continue

        inrange_accs_mask = get_applicable_acceptors(acc_positions, donor_positions)  # Boolean mask of chosen accs

        new_acc_cand_df = pd.concat(
            [new_acc_cand_df, acceptor_positions_df.loc[inrange_accs_mask]])

        # print(f'Portion of surviving acceptor candidates: {sum(inrange_accs_mask) / len(acc_positions)}. '
        #       f'Scaffold {scaffold}')

    new_acc_cand_df.to_csv(acceptor_file, sep=';', index=False)
    print(f'Number of acceptors before: {acc_cand_df.shape[0]}\nNumber of acceptor candidatess after filtering: '
          f'{new_acc_cand_df.shape[0]}')


def is_pairable(don_pos, acc_pos):
    return don_pos + INTRON_MIN_LEN <= acc_pos <= don_pos + INTRON_MAX_LEN


def next_acceptor(acc_positions, acc_idx):
    acc_idx += 1
    acc_pos = acc_positions[acc_idx] if acc_idx < len(acc_positions) else None

    return acc_idx, acc_pos


def get_applicable_acceptors(acc_positions, donor_positions):
    acc_idx, don_idx = 0, 0
    acc_pos = acc_positions[acc_idx]

    inrange_accs = np.zeros(len(acc_positions), dtype=int)

    # Loop until we run out of donors or acceptors
    while don_idx < len(donor_positions) and acc_idx < len(acc_positions):
        don_pos = donor_positions[don_idx]
        if is_pairable(don_pos, acc_pos):
            while acc_idx < len(acc_positions) and is_pairable(don_pos, acc_pos):
                # If acceptor candidate in valid range, flag as OK and shift to next one
                inrange_accs[acc_idx] = True
                acc_idx, acc_pos = next_acceptor(acc_positions, acc_idx)

            # Outside the loop. Now the acceptor is too far. We can move to the next donor
            # .....d.....|..........|..a..
            don_idx += 1

        # -------------------------------  Not pairable cases ---------------------------------
        # Case 1: donor preceding acceptor (which is OK), but distance too long
        # .....d.....|..........|..a..
        elif don_pos + INTRON_MAX_LEN < acc_pos:
            # Just shift to next donor
            don_idx += 1

        # Case 2: donor preceding acceptor (which is OK), but distance too short
        # .....d..a..|..........|.....
        elif don_pos + INTRON_MIN_LEN > acc_pos:
            # Until a pairable acceptor is found, all other ones are labeled as false
            while acc_idx < len(acc_positions) and don_pos + INTRON_MIN_LEN > acc_pos:
                acc_idx, acc_pos = next_acceptor(acc_positions, acc_idx)

        # Case 3: donor is ahead of acceptor (should be vice versa)
        # a........d....|..........|.....
        elif don_pos > acc_pos:
            # Label all acceptor candidates between the current donor and previous acceptor as false
            while don_pos > acc_pos and acc_idx < len(acc_positions):
                acc_idx, acc_pos = next_acceptor(acc_positions, acc_idx)
        else:
            raise ValueError

    return inrange_accs.astype(bool)


if __name__ == '__main__':
    main()
