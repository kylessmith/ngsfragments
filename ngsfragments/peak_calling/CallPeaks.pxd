# cython: language_level=3

import numpy as np
cimport numpy as np
from ailist.LabeledIntervalArray_core cimport LabeledIntervalArray, labeled_aiarray_t, labeled_aiarray_add, labeled_aiarray_init


cdef labeled_aiarray_t *above_zero_regions(double[::1] signal, char *label_name, long start_pos, int shift, int min_length)

cpdef np.ndarray normalize_signal(double[::1] signal, int window=?, bint smooth=?,
								 int smooth_window=?, int polyorder=?, int n_threads=?,
								 bint use_mean=?)

cpdef LabeledIntervalArray call_wps_peaks(np.ndarray wps, long start_pos, str label, int shift=?, int merge_distance=?,
								 int min_length=?, int max_length=?)

cpdef LabeledIntervalArray call_auc_peaks(np.ndarray coverage, long start_pos, str label, int shift=?, int merge_distance=?,
								 int min_length=?, int max_length=?)