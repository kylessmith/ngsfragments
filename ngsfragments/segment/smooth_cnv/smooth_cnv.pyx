#cython: embedsignature=True
#cython: profile=False
#cython: language_level=3

import os
import numpy as np
cimport numpy as np
np.import_array()
cimport cython
from scipy.stats import norm


def trimmed_variance(np.ndarray genomedat, double trim=0.025):
    """
    """

    n = len(genomedat)
    n_keep = int(np.round((1 - 2 * trim) * (n - 1)))

    return inflfact(trim) * np.sum((np.sort(np.abs(np.diff(genomedat)))[:n_keep])**2 / (2 * n_keep))


def inflfact(double trim):
    """
    """

    a = norm.ppf(1 - trim)
    x = np.arange(-a, a, (a*2) / 10001)
    x1 = (x[:10000] + x[1:]) / 2
    
    return 1 / (np.sum(x1**2 * norm.pdf(x1) / (1 - 2 * trim)) * (2 * a / 10000))


def smoothCNV(np.ndarray chroms, np.ndarray genomedat, int smooth_region=10, int outlier_SD_scale=4, int smooth_SD_scale=2, double trim=0.025):
    """
    """

    trimmed_SD = np.sqrt(trimmed_variance(genomedat, trim))
    cdef double outlier_SD = outlier_SD_scale * trimmed_SD
    cdef double smooth_SD = smooth_SD_scale * trimmed_SD
    cdef int k = smooth_region
    cdef int n = len(genomedat)

    cdef long[:] cfrq = np.diff(np.append(np.sort(np.unique(chroms, return_index=True)[1]), np.array([n]))) 
    cdef int nchr = len(cfrq) # to allow for some chrom with all missing

    cdef double[:] genomedat_mem = genomedat
    cdef double[:] sgdat = np.zeros(n, dtype=np.double)

    smoothRL(n, &genomedat_mem[0], nchr, &cfrq[0], &sgdat[0], k, outlier_SD, smooth_SD)

    return np.asarray(sgdat)