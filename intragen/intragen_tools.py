import itertools
from bisect import bisect
from typing import Tuple, List

import numpy as np
import pandas as pd
from pandas import DataFrame


def intraexonic_cuts_in_scaffold(scaff_exon_positions: DataFrame, query_positions: DataFrame) -> Tuple[List[int], int]:
    """
    Determines which of @query_positions lie inside an exon (within a particular scaffold)

    :param scaff_exon_positions: Exon positions in a single scaffold (SINGLE)
    :param query_positions: Positions we want to evaluate
    :return: List of True/False indicators, the length of exon
    """
    start_end_tuplelist = zip(scaff_exon_positions.start, scaff_exon_positions.end)
    start_end_tuplelist = sorted(start_end_tuplelist, key=lambda tup: tup[0])

    start_end_flattened = list(itertools.chain(*start_end_tuplelist))
    exon_length = sum([se[1] - se[0] for se in start_end_tuplelist])

    def inside_exon_cut(cut_start) -> int:
        i = bisect(start_end_flattened, cut_start)
        return i % 2 == 1

    cut_starts = query_positions['start']
    exon_cuts = [inside_exon_cut(cut_start) for cut_start in cut_starts]

    return exon_cuts, exon_length


def join_on_position(cut_coords_file: str, intron_annotations_file: str) -> DataFrame:
    """
    Joins a table of intron cuts (their coordinates) and a table of labeled intron candidates.
    The goal is to be able to see for each candidate his label, prediction and cut decision.

    :param cut_coords_file: CSV file with coordinates of performed cuts
    :param intron_annotations_file: CSV file of labeled intron candidates as well their prediction flags
    :return: DataFrame with columns [scaffold, start, end, label, prediction, cut]
    """
    cut_coords_df = pd.read_csv(cut_coords_file, delimiter=';')
    intron_coords_df = pd.read_csv(intron_annotations_file, delimiter=';')

    cut_coords_df['cut'] = np.ones(len(cut_coords_df))

    joined = intron_coords_df.merge(cut_coords_df, how="left")
    joined[['cut']] = joined[['cut']].fillna(value=-1)

    return joined


def splice_site_positions_dict_to_df(splicesite_positions: dict) -> DataFrame:
    """
    Convert dictionary in form of {scaffold: [123, 12, 34], ...} to a data frame.

    :param splicesite_positions: Dictionary to convert
    :return: DataFrame with 'scaffold' column
    """
    per_scaffold_dfs = []
    for scaffold, positions in splicesite_positions.items():
        data = {'scaffold': [scaffold] * len(positions), 'position': positions}
        per_scaffold_dfs.append(pd.DataFrame(data))

    return pd.concat(per_scaffold_dfs)
