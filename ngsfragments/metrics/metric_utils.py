from typing import Dict
from intervalframe import IntervalFrame
import numpy as np
import pandas as pd


def score_window(scores: pd.Series,
                 windows: IntervalFrame,
                 strand: np.ndarray | None = None,
                 window_size: int = 2000):
    """
    """

    # Iterate over windows
    i = 0
    values = np.zeros((windows.shape[0], window_size + 1))
    for interval in windows.index:
        if (interval.end - interval.start) == (window_size + 1):
            value = scores.loc[interval.start:interval.end-1].values
            if len(value) == (window_size + 1):
                values[i,:] = value
            else:
                #values[i,:] = np.nan
                values[i,:] = 0
        else:
            #values[i,:] = np.nan
            values[i,:] = 0
        i += 1

    # Account for strand
    if strand is not None:
        values[strand == "-", :] = values[strand == "-", ::-1]

    return values