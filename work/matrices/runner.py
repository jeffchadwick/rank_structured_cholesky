#!/usr/bin/env python

import os
import sys

import numpy as np
import numpy.linalg as la
import scipy.io as sio
import gzip

import rschol


try:
    from wget import wget
except:
    def wget(url):
        fname = os.system('wget {0}'.format(url))
        assert os.path.exists(fname), "Download failed"


def download(collection, basename):
    baseURL = 'http://www.cise.ufl.edu/research/sparse/mat/'
    matname = '{0}.mat'.format(basename)
    if not os.path.exists(matname):
        assert collection is not None, "Invalid collection for download"
        wget('{0}/{1}/{2}.mat'.format(baseURL, collection, basename))


def setup(collection, basename):
    download(collection, basename)
    m = sio.loadmat('{0}.mat'.format(basename))
    problem = m['Problem']
    A = problem['A'][0,0]
    if 'b' in problem:
        b = problem['b'][0,0][:,0]
    else:
        b = np.ones((A.shape[0],))
    return A, b


def run(collection, basename):
    A, b = setup(collection, basename)
    x = rschol.solve(A, b)
    r = A*x-b
    print("Relative residual norm: {0}".format(la.norm(r)/la.norm(b)))


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
        print("Format: runner.py command matname")
