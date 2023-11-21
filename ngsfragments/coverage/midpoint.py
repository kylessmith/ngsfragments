import numpy as np
import pandas as pd
from ailist import LabeledIntervalArray

# Local imports
from ..fragments import Fragments


def midpoint(intervals: Fragments | LabeledIntervalArray,
            min_length: int | None = None,
            max_length: int | None = None):
    """
    Calculate midpoint

    Parameters
    ----------
        intervals : Fragments
            Fragments object to calculate midpoint for
        min_length : int
            Minimum length of intervals to include [default = None]
        max_length : int
            Maximum length of intervals to include [default = None]
        
    Returns
    -------
        mp : pandas.Series
            Midpoint coverage
    """

    # Get intervals
    if isinstance(intervals, Fragments):
        intervals = intervals.frags

    # Calculate
    mp = intervals.midpoint_coverage(min_length, max_length)
    

    return mp