# cython: language_level=3

from cython.parallel import prange
import numpy as np
cimport numpy as np
cimport cython


cdef extern from "sklist_pd.c":
	# C is include here so that it doesn't need to be compiled externally
	pass

cdef extern from "sklist_pd.h":
	void rolling_median(const double values[], double medians[], int n, int window) nogil
	void rolling_median_shifted(const double values[], double medians[], int n, int window, int start) nogil
	

cpdef np.ndarray running_median(const double[:] values, int window, int n_threads=?, int n_chunks=?)