# Copyright (c) 2015, Jeff Chadwick and David Bindel
# All rights reserved.
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
# 
# 1. Redistributions of source code must retain the above copyright notice,
# this list of conditions and the following disclaimer.
# 
# 2. Redistributions in binary form must reproduce the above copyright notice,
# this list of conditions and the following disclaimer in the documentation
# and/or other materials provided with the distribution.
# 
# 3. Neither the name of the copyright holder nor the names of its
# contributors may be used to endorse or promote products derived from this
# software without specific prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
# THE POSSIBILITY OF SUCH DAMAGE.

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
