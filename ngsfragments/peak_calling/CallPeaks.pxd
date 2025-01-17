# cython: language_level=3

import numpy as np
cimport numpy as np
np.import_array()
from ailist.LabeledIntervalArray_core cimport LabeledIntervalArray, labeled_aiarray_t, labeled_aiarray_add, labeled_aiarray_init


cdef labeled_aiarray_t *above_zero_regions(double[::1] signal, const char *label_name, long start_pos)

cpdef np.ndarray normalize_signal(double[::1] signal, int window=?, bint smooth=?,
								 int smooth_window=?, int polyorder=?, int n_threads=?,
								 bint use_mean=?)
