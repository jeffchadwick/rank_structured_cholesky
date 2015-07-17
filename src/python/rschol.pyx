import numpy as np
cimport numpy as np
from scipy.sparse import csr_matrix

cdef extern from "pcg_pysolver.h":
    void pcg_pysolve(int nrow, const int* jc, const int* ir, const double* pr,
                     const double* rhs, double* solution)

def solve(A, np.ndarray[np.float_t, ndim=1] b):
    cdef np.ndarray[np.int32_t, ndim=1] indices, indptr
    cdef np.ndarray[np.float64_t, ndim=1] data
    if not isinstance(A, csr_matrix):
        A = csr_matrix(A)
    indices = A.indices.astype(np.int32)
    indptr = A.indptr.astype(np.int32)
    data = A.data.astype(np.float64)
    cdef np.ndarray x = b.copy()
    pcg_pysolve(A.shape[0],
                <int*> indptr.data,
                <int*> indices.data,
                <double*> data.data,
                <double*> b.data,
                <double*> x.data)
    return x
