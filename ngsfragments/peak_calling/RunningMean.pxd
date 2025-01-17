# cython: language_level=3

from cython.parallel import prange
import numpy as np
cimport numpy as np
np.import_array()
cimport cython


cdef extern from "running_mean.c":
	# C is include here so that it doesn't need to be compiled externally
	pass

cdef extern from "running_mean.h":
	void rolling_mean_array(double values[], double means[], int length, int window_size) nogil
	

cdef void c_running_mean(const double[::1] values, double[::1] means, int n, int window)
cpdef np.ndarray running_mean(const double[::1] values, int window)