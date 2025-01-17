from ailist import LabeledIntervalArray
import gzip
from collections import defaultdict

# Local imports
from ..fragments import Fragments


def read_fragments(filename: str,
                   n_frags: int = 10000,
                   genome_version: str = "hg38",
                   verbose: bool = False):
    """
    """

    # Initialize fragment record
    record = defaultdict(LabeledIntervalArray)

    # Iterate over fragments file
    if verbose: print("Reading "+filename, flush=True)
    with gzip.open(filename, 'rt') as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.split("\t")
            record[fields[3]].add(int(fields[1]), int(fields[2]), fields[0])

    # Convert to fragments
    if verbose: print("Converting to fragments", flush=True)
    keys = list(record.keys())
    for key in keys:
        if len(record[key]) < n_frags:
            record[key].close()
            del record[key]
        else:
            record[key] = Fragments(record[key], genome_version = genome_version)

    return record