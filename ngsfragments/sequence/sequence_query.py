import pysam
import numpy as np
import pandas as pd
from ailist import LabeledIntervalArray


def query_locations(filename: str,
                    intervals: LabeledIntervalArray):
    """
    """

    samfile = pysam.AlignmentFile(filename, "rb")
    it = samfile.pileup(reference=region_chrom, start=region_start, end=region_end,
							truncate=False, max_depth=80000, ignore_orphans=True, min_base_quality=min_qual,
							min_mapping_quality=min_mapq)