import numpy as np
import pandas as pd
from ailist import LabeledIntervalArray

# Local imports
from ..fragments import Fragments


def kmer_query(intervals: Fragments | LabeledIntervalArray,
               genome_version: str = "hg19",
               k: int = 2,
               last_n: int = 0) -> pd.DataFrame:
    """
    Query k-mers in intervals

    Parameters
    ----------
        intervals : :class:`~fragments.fragments`
            :class:`~fragments.fragments` object
        genome_version : str
            Genome version
        k : int
            K-mer length
        last_n : int
            Last n bases to query

    Returns
    -------
        kmer_df : :class:`~pandas.DataFrame`
            K-mer query results
    """

    # Import genome
    if isinstance(intervals, Fragments):
        genome = intervals.genome
        intervals = intervals.frags
    else:
        import genome_info
        genome = genome_info.GenomeInfo(genome_version)

    # Calculate k-mers
    kmer_dict = genome.kmers(intervals, k=k, last_n=last_n)

    # Convert to DataFrame
    kmer_df = pd.DataFrame(kmer_dict.values(), index=kmer_dict.keys(), columns=["counts"])
    kmer_df.loc[:,"proportion"] = kmer_df.loc[:,"counts"].values / kmer_df.loc[:,"counts"].values.sum()

    return kmer_df


