#!/usr/bin/env python

# NB: I would have liked to use gzip.open, but kept getting an
# error "invalid compressed data--length error" on decompression.

import numpy as np
import scipy.io as sio

m = sio.loadmat('1138_bus.mat')
A = m['Problem']['A'][0,0]

# nnz, m, n are uint64; offsets/indices are int32; data is double
with open('1138_bus_system.bcsm', 'wb') as f:
    f.write(np.array([A.nnz, A.shape[0], A.shape[1]], dtype=np.uint64))
    f.write(A.indptr)
    f.write(A.indices)
    f.write(A.data)

with open('1138_bus_rhs.vector', 'wb') as f:
    f.write(np.array([A.shape[0]], dtype=np.int32))
    f.write(np.ones((A.shape[0],)))
