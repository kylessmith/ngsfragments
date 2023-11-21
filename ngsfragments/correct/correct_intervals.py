import pandas as pd
import numpy as np
from intervalframe import IntervalFrame

# Local imports
from .correction import correct_counts


def calculate_interval_bias(intervals: IntervalFrame,
                            column: str,
                            cnv_bins: IntervalFrame | None = None,
                            genome_version: str = "hg19",
                            include_blacklist: bool = True,
                            include_repeat: bool = True,
                            include_gc: bool = True) -> None:
    """
    Calculate bias per interval
    
    Parameters
    ----------
        intervals : IntervalFrame
            Labeled intervals
        column : str
            Column to use for bias correction
        genome_version : str
            Genome version name
        include_blacklist : bool
            Flag to include blacklist
        include_repeat : bool
            Flag to include repeat
        include_gc : bool
            Flag to include gc

    Returns
    ----------
        None
    """

    # Assign genome
    if genome_version == "hg19":
        from hg19genome import calculate_bias
    elif genome_version == "hg38":
        from hg38genome import calculate_bias

    # Initialize bias records
    bias_record = calculate_bias(intervals.index,
                                 include_blacklist,
                                 include_repeat,
                                 include_gc)
    
    # Remove blacklist
    chosen = bias_record.df.loc[:,"blacklist"].values < 0.1
    bias_record = bias_record.iloc[chosen,:]
    intervals = intervals.iloc[chosen,:]
    bias_record.drop_columns(["blacklist"])
    
    # Calculate cnv
    if cnv_bins is not None:
        bias_record.annotate(cnv_bins,
                             column = "ratios", 
                             method = "mean",
                             column_name = "cnv_mean")
        # Filter nans
        chosen = ~pd.isnull(bias_record.df.loc[:,"cnv_mean"].values)
        bias_record = bias_record.iloc[chosen,:]
        intervals = intervals.iloc[chosen,:]

    # Correct
    intervals.df.loc[:,"corrected_values"] = correct_counts(intervals.df.loc[:,column].values,
                                                            bias_record.df)


    return intervals