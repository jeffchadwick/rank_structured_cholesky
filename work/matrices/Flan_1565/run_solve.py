#!/usr/bin/python

import sys
import os

systemURL = 'http://www.cise.ufl.edu/research/sparse/mat/Janna/Flan_1565.mat';
systemFile = 'Flan_1565.mat';

matrixFile = 'Flan_1565_system.bcsm';
rhsFile = 'Flan_1565_rhs.vector';

buildPath = '../../../gcc-build';

# Optionally overwrite the build path
if ( len( sys.argv ) > 1 ):
  buildPath = sys.argv[1];

if not os.path.exists( systemFile ):
  # Download the file
  os.system( 'wget %s' % systemURL );

if not os.path.exists( matrixFile ) or not os.path.exists( rhsFile ):
  # Run the matlab conversion script
  os.system( 'matlab -r convert' );

if not os.path.exists( '%s.gz' % matrixFile ):
  # Compress
  os.system( 'gzip %s' % matrixFile );

# Run the solver
os.system( '%s/pcg_solver Flan_1565' % buildPath );
