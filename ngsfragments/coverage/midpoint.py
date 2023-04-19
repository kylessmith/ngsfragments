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
    """

    # Get intervals
    if isinstance(intervals, Fragments):
        intervals = intervals.frags

    # Calculate
    mp = intervals.midpoint_coverage(min_length, max_length)
    

    return mp