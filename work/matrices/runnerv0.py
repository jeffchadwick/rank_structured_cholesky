#!/usr/bin/env python

import os
import sys

import numpy as np
import scipy.io as sio
import gzip


# BUILD_PATH="../../../../rsc-build"
BUILD_PATH="/Users/dbindel/local/rsc/bin"


def mat_to_bcsm(basename, write_rhs=True):
    """Write UF sparse collection .mat file to internal format.

    NB: We probably ought to just hook up a Python front-end to the solver.
    """

    m = sio.loadmat('{0}.mat'.format(basename))
    problem = m['Problem']
    A = problem['A'][0,0]

    # nnz, m, n are uint64; offsets/indices are int32; data is double
    with gzip.open('{0}_system.bcsm.gz'.format(basename), 'wb') as f:
        hdr = np.array([A.nnz, A.shape[0], A.shape[1]], dtype=np.uint64)
        f.write(hdr.tostring())
        f.write(A.indptr.tostring())
        f.write(A.indices.tostring())
        f.write(A.data.tostring())

    if write_rhs:
        if 'b' in problem:
            b = problem['b'][0,0][:,0]
        else:
            b = np.ones((A.shape[0],))
        with open('{0}_rhs.vector'.format(basename), 'wb') as f:
            hdr = np.array([A.shape[0]], dtype=np.int32)
            f.write(hdr)
            f.write(b)


def download(collection, basename):
    baseURL = 'http://www.cise.ufl.edu/research/sparse/mat/'
    matname = '{0}.mat'.format(basename)
    if not os.path.exists(matname):
        assert collection is not None, "Invalid collection for download"
        os.system('wget {0}/{1}/{2}.mat'.format(baseURL, collection, basename))


def setup(collection, basename):
    download(collection, basename)
    if not os.path.exists('{0}_system.bcsm.gz'.format(basename)):
        mat_to_bcsm(basename)


def run(collection, basename):
    setup(collection, basename)
    rundir = '{0}_run'.format(basename)
    if not os.path.exists(rundir):
        os.mkdir(rundir)
    os.chdir(rundir)
    os.system('{0}/pcg_solver ../{1}'.format(BUILD_PATH, basename))


if __name__ == "__main__":
    cmds = {'run': run, 'download': download, 'setup': setup}
    if len(sys.argv) == 3:
        cmd, name = sys.argv[1:]
        name_parts = name.split('/')
        assert len(name_parts) <= 2, 'Invalid name'
        if len(name_parts) == 2:
            collection = name_parts[0]
            basename = name_parts[1]
        else:
            collection = None
            basename = name_parts[0]
        cmds[cmd](collection, basename)
    else:
        print("Format: translate.py command matname")
