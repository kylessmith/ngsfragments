import os
from typing import Dict
import numpy as np
import pandas as pd
from ailist import LabeledIntervalArray
from collections import OrderedDict


def calculate_chromosome_shift(genome: Dict):
    """
    """
    
	# Initialize
    chrom_shift = OrderedDict()
    
	# Record shift across genome
    shift = 0
    for chrom in genome:
        chrom_shift[chrom] = shift
        shift += genome[chrom]
        
    return chrom_shift
	

class Fragments(object):
	"""
	Python class to store DNA fragments in Augmented Interval List

	:class:`~fragments.fragments` stores a list of NGS fragments
	
	Params
	------
		None
	"""

	def __init__(self,
				 frags: LabeledIntervalArray | str = "",
				 sam_file: str = "",
				 genome: str = "hg19"):
		"""
		Initialize fragments class

		Parameters
		----------
			frags : :class:`~ailist.LabeledIntervalArray`
				:class:`~ailist.LabeledIntervalArray`
			sam_file : str
				SAM file
			genome : str
				Genome file

		Returns
		-------
			None
		"""

		# Set file name
		self.sam_file = sam_file

		# Set fragments
		if isinstance(frags, LabeledIntervalArray):
			# Enforce construction
			self.frags = frags
			self.frags.construct()
			self.frags.freeze()
		
		#self.sam_file = sam_file
		if genome == "hg19":
			try:
				import hg19genome
				self.genome = hg19genome.Hg19Genome()
				self.genome_version = "hg19"
			except ImportError:
				print("hg19genome not installed. Please install hg19genome: 'pip install hg19genome'")
		else:
			raise NotImplementedError("Only hg19 genome is currently supported")

	def __repr__(self):
		"""
		Represent fragments
		"""

		repr_string = "Fragments: {} fragments\ngenome: {}".format(self.size, self.genome_version)

		return repr_string

	
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


	def intersect(self,
				  chrom: str,
				  start: int,
				  end: int):
		"""
		Intersect fragments with interval
		"""

		# Overlap
		frags = self.frags.intersect(start, end, chrom)

		# Build
		fragments = Fragments(frags,
							  self.sam_file,
							  self.genome)
		
		return fragments

	
	def filter(self,
				min_length: int | None = None,
				max_length: int | None = None):
		"""
		Filter fragments

		Parameters
		----------
			min_length : int
				Minimum fragment length
			max_length : int
				Maximum fragment length
		
		Returns
		-------
			filtered_frags : :class:`~fragments.Fragments`
				:class:`~fragments.Fragments`
		"""

		# Filter
		filtered_intervals = self.frags.filter(min_length, max_length)

		# Build
		filtered_frags = Fragments(filtered_intervals,
							  		self.sam_file,
							  		self.genome)
		
		return filtered_frags


	def downsample(self,
				   n: int = None,
				   proportion: float = None):
		"""
		Downsample fragments

		Parameters
		----------
			n : int
				Number of fragments to use
			proportion : float
				Proportion of fragments to use

		Returns
		-------
			None
		"""

		# Calculate proportion
		if n is not None:
			proportion = int(self.frags.size * proportion)

		# Downsample
		if proportion != 1:
			self.frags = self.frags.downsample(proportion)


	def length_dist(self):
		"""
		Calculate fragment length distribution

		Parameters
		----------
			None

		Returns
		-------
			distribution : pd.Series
				Fragment length distribution

		"""
		
		# Calculate fragment length distribution
		distribution = self.frags.length_dist()
		distribution = pd.Series(distribution)
		distribution.name = self.sam_file

		return distribution


	def simulate(self):
		"""
		Simulate fragments
		"""

		sim = Fragments()

		# Initialize attributes
		sim.sam_file = "Simulation"
		sim.genome = self.genome
		sim.proportion = self.proportion
		sim.verbose = self.verbose
		sim.data_dir = self.data_dir
		sim.chrom_shift = self.chrom_shift

		# Simulate
		sim.frags = self.frags.simulate()

		return sim