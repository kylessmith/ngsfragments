#cython: language_level=3

from cython.parallel import prange
import numpy as np
cimport numpy as np
cimport cython
	

cdef void c_running_mean(const double[::1] values, double[::1] means, int n, int window):
	rolling_mean_array(&values[0], &means[0], n, window)


cpdef np.ndarray running_mean(const double[::1] values, int window):
	"""
	Calculate a running mean across values
	
	Arguments
		values: double array/numpy array
		window: int
		
	Returns
		means: numpy array
	"""
	# Initialize variables
	cdef int n = len(values)
	cdef np.ndarray means = np.zeros(n, dtype=np.double)
	cdef double[::1] means_mem = means
	
	c_running_mean(values, means_mem, n, window)

	return means