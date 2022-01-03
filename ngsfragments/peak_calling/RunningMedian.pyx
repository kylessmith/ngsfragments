#cython: language_level=3

from cython.parallel import prange
import numpy as np
cimport numpy as np
cimport cython
	

@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing.
@cython.initializedcheck(False) # Deactivating memoryview init checking
cdef void c_running_median(const double[:] values, double[:] medians, int n, int window) nogil:
	rolling_median(&values[0], &medians[0], n, window)
	

@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing.
@cython.initializedcheck(False) # Deactivating memoryview init checking
cdef void c_running_median_shifted(const double[:] values, double[:] medians,
								   int window, int n_threads, int n_chunks) nogil:
	"""
	"""
	# Initialize variables
	cdef int chunk_size
	cdef int chunk_start
	cdef int chunk_end
	cdef int i
	cdef int n = len(values)
	
	# Calculate number of chunks
	chunk_size = int(n / n_chunks)
	
	# Iterate over chunks and calculate medians
	for i in prange(0, n_chunks, nogil=True, num_threads=n_threads, schedule='static'):
		chunk_start = (i * chunk_size)
		chunk_end = min(n, chunk_start + chunk_size + window)
		rolling_median_shifted(&values[0], &medians[0], chunk_end, window, chunk_start)


@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing.
cpdef np.ndarray running_median(const double[:] values, int window, int n_threads=1, int n_chunks=10):
	"""
	Calculate a running median across values
	
	Arguments
		values: double array/numpy array
		window: int
		n_threads: int
		n_chunks: int
		
	Returns
		medians: numpy array
	"""
	# Initialize variables
	cdef int n = len(values)
	cdef double[:] medians = np.zeros(n, dtype=np.double)
	
	if n_threads == 1:
		c_running_median(values, medians, n, window)
	else:
		c_running_median_shifted(values, medians, window, n_threads, n_chunks)

	return np.asarray(medians)