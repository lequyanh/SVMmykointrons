import sys
from bisect import bisect_left, bisect_right

import pandas as pd


def get_near_positions(positions, min_pos, max_pos):
    """
    Returns all positions P such that `min_pos` <= P <= `max_pos`.

    Parameters
    ----------
    positions : [int]
        an ascending list of available positions
    min_pos
        an integer specifying the lowest position to consider
    max_pos
        an integer specifying the greatest position to consider
    """
    start_index = bisect_left(positions, min_pos)
    end_index = bisect_right(positions, max_pos)
    return [positions[i] for i in range(start_index, end_index)]


def read_data(filename) -> pd.DataFrame:
    """
    Reads data and returns them sorted by `scaffold` and `position`.

    Parameters
    ----------
    filename: str
        Name of a CSV file with `scaffold` and `position` columns. 
        The file should contain positions of positively classified splice sites only.
    """
    return pd.read_csv(filename, sep=';') \
        .sort_values(by=['scaffold', 'position']) \
        .reset_index(drop=True)


if __name__ == '__main__':
    # CSV file containing positions of positively classified donors
    donor_predictions = sys.argv[1]
    # CSV file containing positions of positively classified acceptors
    acceptor_predictions = sys.argv[2]
    # minimal intron length to consider
    intron_min_len = int(sys.argv[3])
    # maximal intron length to consider
    intron_max_len = int(sys.argv[4])

    assert intron_min_len <= intron_max_len, 'the minimal length must not be greater than the maximal length'

    donors = read_data(donor_predictions)
    acceptors = read_data(acceptor_predictions).set_index('scaffold')

    for _, row in donors.iterrows():
        scaffold, donor_position = row[['scaffold', 'position']]
        acceptor_positions = acceptors.loc[[scaffold], 'position']

        # for each donor find acceptors in allowed distance

        acceptor_positions = get_near_positions(acceptor_positions.array,
                                                donor_position + intron_min_len,
                                                donor_position + intron_max_len)

        for acceptor_position in acceptor_positions:
            print(scaffold, donor_position, acceptor_position + 1)
