import os
import numpy as np
import pandas as pd
import pysam
from .read_sam import read_fragments
from intervalframe import IntervalFrame
from collections import OrderedDict


def get_genome(sam_file_fn):
	"""
	Read genome file
	"""

	# Open SAM file
	sam_file = pysam.Samfile(sam_file_fn)

	# Read chromosomes in SAM file
	genome = {}
	genome.update(list(zip(sam_file.references, sam_file.lengths)))

	# Close SAM file
	sam_file.close()

	return genome


def calculate_chromosome_shift(genome):
    """
    """
    
    chrom_shift = OrderedDict()
    
    shift = 0
    for chrom in genome:
        chrom_shift[chrom] = shift
        shift += genome[chrom]
        
    return chrom_shift
	

class fragments:
	"""
	Python class to store DNA fragments in Augmented Interval List

	:class:`~fragments.fragments` stores a list of NGS fragments
	
	Params
	------
		None
	"""
	def __init__(self, sam_fn=None, chrom=None, sam_file=None, min_size=1, max_size=1000,
				 paired=True, qcfail=False, mapq_cutoff=25, verbose=False,
				 n_jobs=1, proportion=1.0):
		"""
		Initialize fragments class

		Params
		------
			sam_fn
				str
			chrom
			sam_file
			min_size
			max_size
			paired
			qcfail
			mapq_cutoff
			verbose
			n_jobs
			proportion
		"""

		# Initialize attributes
		self.sam_file = ""
		self.genome = {}
		self.proportion = proportion
		self.verbose = verbose

		# Set data directory
		self.data_dir = os.path.split(os.path.split(os.path.realpath(__file__))[0])[0]
		self.data_dir = os.path.join(self.data_dir, "data")

		# If data is a string open as file
		if sam_fn is not None:
			# Assign bam file
			self.sam_file = sam_fn
			self.genome = get_genome(self.sam_file)
			self.chrom_shift = calculate_chromosome_shift(self.genome)

			# Add fragment intervals
			if verbose: print("Reading")
			frags = read_fragments(self.sam_file, min_size, max_size, paired, qcfail, mapq_cutoff, proportion)
			frags.construct()
			frags.freeze()

			# Assign
			self.frags = frags

	
	@property
	def size(self):

		return len(self.frags)

	@property
	def n_chroms(self):
		return len(self.frags.unique_labels)

	@property
	def chroms(self):
		return self.frags.unique_labels

	
	def __len__(self):
		"""
		Get length from self.frags
		"""
		
		return self.size


	def downsample(self, n=None, proportion=None):
		"""
		"""

		if n is not None:
			proportion = int(self.frags.size * proportion)

		self.frags = self.frags.downsample(proportion)


	def length_dist(self):
		"""
		Calculate fragment length distribution

		Parameters
		----------
			None

		Returns
		-------
			distribution : numpy.ndarray
				Fragment length distribution

		"""
		
		# Calculate fragment length distribution
		distribution = self.frags.length_dist()

		return distribution


	def bin_counts(self, bin_size=100000, min_length=None, max_length=None):
		"""
		Calculate coverage of fragments

		Parameters
		----------
			bin_size : int
				Size of each bin
			min_length : int
				Minimum length of each fragment to consider
			max_length : int
				Maximum length of each fragment to consider

		Returns
		-------
			bins : IntervalFrame
				Bins containing number of hits
		"""

		# Determine nhits per bin
		bins = self.frags.bin_nhits(bin_size, min_length, max_length)

		# Create interval frame
		bins_iframe = IntervalFrame(bins[0])
		bins_iframe.df.loc[:,"counts"] = bins[1]

		return bins_iframe

	
	def wps(self, chrom=None, protection=60, min_length=None, max_length=None):
		"""
		Calculate Window Protection Score for each position in AIList range

		Parameters
		----------
			protection : int
				Protection window to use
			min_length : int
				Minimum length of intervals to include [default = None]
			max_length : int
				Maximum length of intervals to include [default = None]
			label : str
				Label for hierarchical indexing

		Returns
		-------
			scores : dict of pandas.Series or pandas.Series
				Position on index and WPS as values

		"""

		if chrom is None:
			# Calculate WPS without label
			scores = self.frags.wps(protection, min_length, max_length)

			# Chop up into dictionaries
			scores_dict = {}
			ranges = self.frags.label_ranges
			count = 0
			for chrom in self.frags.unique_labels:
				length = ranges[chrom][1] - ranges[chrom][0]
				scores_dict[str(chrom)] = pd.Series(scores[count:count+length],
													index=range(ranges[chrom][0],
																ranges[chrom][1]))
				count += length
		
		else:
			# Calculate WPS without label
			scores = self.frags.get(chrom).wps(protection, min_length, max_length)

			# Chop up into dictionaries
			scores_dict = {}
			ranges = self.frags.label_ranges
			scores_dict[chrom] = pd.Series(scores,
										   index=range(ranges[chrom][0],
													   ranges[chrom][1]))
		
		return scores_dict