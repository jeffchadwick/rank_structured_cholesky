//////////////////////////////////////////////////////////////////////
// Copyright (c) 2015, Jeff Chadwick and David Bindel
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
// 
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
// 
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
// 
// 3. Neither the name of the copyright holder nor the names of its
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
// FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
// COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
// INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
// BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
// ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////
// MathUtil.cpp: Implementation of the MathUtil class
//
//////////////////////////////////////////////////////////////////////

#include "MathUtil.h"

#include "MERSENNETWISTER.h"

//MERSENNETWISTER MathUtil::GENERATOR = MERSENNETWISTER( 12345 );
MERSENNETWISTER MathUtil::GENERATOR = MERSENNETWISTER( 152 );

//////////////////////////////////////////////////////////////////////
// Generates a random vector with entries in the given range
//////////////////////////////////////////////////////////////////////
VECTOR MathUtil::randomVector( int sz, Real rangeMin, Real rangeMax )
{
  VECTOR               result( sz );

  MERSENNETWISTER      generator;

  for ( int i = 0; i < sz; i++ )
  {
    result( i ) = rangeMin + generator.rand( rangeMax - rangeMin );
  }

  return result;
}

//////////////////////////////////////////////////////////////////////
// Generates a random matrix with entries in the given range
//////////////////////////////////////////////////////////////////////
MATRIX MathUtil::randomMatrix( int nRows, int nCols,
                               Real rangeMin, Real rangeMax )
{
  MATRIX               result( nRows, nCols );

  MERSENNETWISTER      generator;

  for ( int i = 0; i < nRows; i++ )
  for ( int j = 0; j < nCols; j++ )
  {
    result( i, j ) = rangeMin + generator.rand( rangeMax - rangeMin );
  }

  return result;
}

//////////////////////////////////////////////////////////////////////
// Generates a random Gaussian matrix with given mean and
// standard deviation
//////////////////////////////////////////////////////////////////////
MATRIX MathUtil::randomGaussianMatrix( int nRows, int nCols,
                                       Real mean, Real stddev )
{
  MATRIX               result( nRows, nCols );

  randomGaussianMatrix( nRows, nCols, result.data(), mean, stddev );

  return result;
}

//////////////////////////////////////////////////////////////////////
// Generates a random Gaussian matrix with given mean and
// standard deviation
//////////////////////////////////////////////////////////////////////
void MathUtil::randomGaussianMatrix( int nRows, int nCols,
                                     Real *matrix,
                                     Real mean, Real stddev )
{
  Real                 variance = stddev * stddev;

  for ( int i = 0; i < nRows; i++ )
  for ( int j = 0; j < nCols; j++ )
  {
    MATRIX::access( matrix, nRows, nCols, i, j )
                = GENERATOR.randNorm( mean, variance );
  }
}

//////////////////////////////////////////////////////////////////////
// Constructor
//////////////////////////////////////////////////////////////////////
MathUtil::MathUtil()
{
}

//////////////////////////////////////////////////////////////////////
// Destructor
//////////////////////////////////////////////////////////////////////
MathUtil::~MathUtil()
{
}

