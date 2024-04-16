from typing import Dict
from intervalframe import IntervalFrame
import numpy as np
import pandas as pd


def nfr(values, center_x=1000, scale=False):
    """
    Calculate Nucleosome Free Region score

    Params
    ------
        values
            numpy.ndarray
        center_x
            int
        smooth
            bool
        scale
            bool

    Returns
    -------
        nfr_score
            float
        pos_score
            float
    """
    
    # Scale values between 0 and 1
    if scale:
        values = (values - np.min(values)) / (np.max(values) - np.min(values))

    # Calculate NFR
    nf = np.sum(values[center_x - 50:center_x + 50])
    n1 = np.sum(values[center_x - 200:center_x - 50])
    n2 = np.sum(values[center_x + 50:center_x + 200])
    nfr_score = np.log2(nf) - np.log2((n1 + n2) / 2)
    
    # Find position score
    control_mean = (np.sum(values[:100]) + np.sum(values[-100:])) / 200
    pos_score = np.log2(values[center_x] / control_mean)		
    
    return nfr_score, pos_score


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