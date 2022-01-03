#cython: language_level=3

import numpy as np
cimport numpy as np
cimport cython
from ailist.LabeledIntervalArray_core cimport LabeledIntervalArray, labeled_aiarray_t, labeled_aiarray_add, labeled_aiarray_init
from .RunningMedian cimport running_median
from .RunningMean cimport running_mean
from scipy.signal import savgol_filter

from libc.stdint cimport uint32_t, uint8_t, uint64_t, int64_t, uint16_t
from libc.stdio cimport printf


@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing.
cdef labeled_aiarray_t *above_zero_regions(double[::1] signal, char *label_name, long start_pos, int shift, int min_length):
	"""
	Find regions of above zero
	"""
	# Initialize variables
	cdef labeled_aiarray_t *cregions = labeled_aiarray_init()
	cdef int start = -1
	cdef int end
	cdef int length
	cdef int i
	cdef bint above

	# Iterate over signal to find regions above zero
	for i in range(signal.size):
		# Check if i is above 0
		if signal[i] > 0:
			above = True
			if start == -1:
				start = i + shift
		else:
			above = False
			# Check if region ends
			if start != -1:
				end = i + shift
				length = end - start
				if length >= min_length:
					labeled_aiarray_add(cregions, start_pos + start, start_pos + end, label_name)
			# Reset start
			start = -1

	return cregions
		

@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing.
cpdef np.ndarray normalize_signal(double[::1] signal, int window=1000, bint smooth=True,
								 int smooth_window=21, int polyorder=2, int n_threads=1, bint use_mean=True):
	"""
	"""

	# Determine if signal is long enough
	if signal.size < window + 1:
		return np.asarray(signal)

	# Intitialize variables
	cdef np.ndarray smooth_signal

	# Determine rolling median/mean
	cdef np.ndarray rolling_values
	if use_mean:
		rolling_values = running_mean(signal, window)
	else:
		rolling_values = running_median(signal, window, n_threads=n_threads)
	
	# Subtract median
	cdef np.ndarray norm_signal = signal - rolling_values
	
	# Smooth with savgol_filter
	if smooth:
		norm_signal = savgol_filter(norm_signal, window_length=smooth_window, polyorder=polyorder)
	
	return norm_signal

	
@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing.
cpdef LabeledIntervalArray call_wps_peaks(np.ndarray wps, long start_pos, str label, int shift=0, int merge_distance=5,
										  int min_length=50, int max_length=150):
	"""
	"""
	# Initialize variables
	cdef LabeledIntervalArray peaks = LabeledIntervalArray()
	cdef labeled_aiarray_t *cpeaks
	
	# Determine regions above zero
	cdef bytes label_name = label.encode()
	cpeaks = above_zero_regions(wps, label_name, start_pos, shift, min_length)
	peaks.set_list(cpeaks)
	
	# Merge nearby regions
	if peaks.size > 1:
		peaks = peaks.merge(merge_distance)
	
	# Filter peaks by length
	peaks = peaks.filter(min_length, max_length)

	return peaks
	

@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing.
cpdef LabeledIntervalArray call_auc_peaks(np.ndarray coverage, long start_pos, str label, int shift=0, int merge_distance=5,
										  int min_length=50, int max_length=150):
	"""
	Call peaks similar to SEACR
	"""
	# Initialize variables
	cdef LabeledIntervalArray peaks = LabeledIntervalArray()
	cdef labeled_aiarray_t *cpeaks
	
	# Determine regions above zero
	cdef bytes label_name = label.decode()
	cpeaks = above_zero_regions(coverage, label_name, start_pos, shift, 0)
	peaks.set_list(cpeaks)
	
	# Merge nearby regions
	peaks = peaks.merge(merge_distance)
	
	# Filter peaks by length
	peaks = peaks.filter(min_length, max_length)

	return peaks