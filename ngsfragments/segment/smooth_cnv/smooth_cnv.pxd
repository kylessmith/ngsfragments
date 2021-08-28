import numpy as np
cimport numpy as np
cimport cython


cdef extern from "smoothCNV.c":
	# C is include here so that it doesn't need to be compiled externally
	pass

cdef extern from "smoothCNV.h":
	void smoothRL(int n, double *gdat, int nchr, long *cfrq, double *sgdat, int k, double oSD, double sSD)


