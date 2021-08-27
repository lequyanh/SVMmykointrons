#!/usr/bin/env python3


"""Classify introns with SVM.

Usage:
  classify-introns-stub.py <data-file-name>
"""
import logging
import sys

import pandas as pd
from docopt import docopt

from tools import performance_metrics

if __name__ == "__main__":
    arguments = docopt(__doc__, version='1.0')

    data_file = arguments['<data-file-name>']
    data = pd.read_csv(data_file, sep=';')

    predictions = [1] * data.shape[0]
    data.assign(prediction=pd.Series(predictions)).to_csv(sys.stdout, sep=';', index=False)

    if 'label' in data:
        logging.basicConfig(
            level=logging.INFO,
            filename=f'intron-accuracy.log',
            filemode='w'
        )
        logging.info(f'Validation of NN models on testset {data_file}')

        metrics_data = performance_metrics(
            data.label,
            predictions,
            imbalance_rat=None
        )

        metrics_str = '\n'.join(metrics_data)
        logging.info(f'Performance metrics:\n'
                     f'{metrics_str}')
