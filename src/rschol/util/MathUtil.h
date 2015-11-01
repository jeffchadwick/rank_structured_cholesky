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
// MathUtil.h: Interface for the MathUtil class
//
//////////////////////////////////////////////////////////////////////

#ifndef MATH_UTIL_H
#define MATH_UTIL_H

#include <rschol/linearalgebra/MATRIX.h>
#include <rschol/linearalgebra/MATRIX3.h>
#include <rschol/linearalgebra/VECTOR.h>

#include <rschol/SETTINGS.h>
#include <rschol/TYPES.h>

#include <rschol/util/MERSENNETWISTER.h>

//////////////////////////////////////////////////////////////////////
// MathUtil class
//
// General utility functions for matrices, vectors, etc.
//////////////////////////////////////////////////////////////////////
class MathUtil {
	public:
    // Generates a random vector with entries in the given range
    static VECTOR randomVector( int sz, Real rangeMin, Real rangeMax );

    // Generates a random matrix with entries in the given range
    static MATRIX randomMatrix( int nRows, int nCols,
                                Real rangeMin, Real rangeMax );

    // Generates a random Gaussian matrix with given mean and
    // standard deviation
    static MATRIX randomGaussianMatrix( int nRows, int nCols,
                                        Real mean = 0.0, Real stddev = 1.0 );

    // Generates a random Gaussian matrix with given mean and
    // standard deviation
    static void randomGaussianMatrix( int nRows, int nCols,
                                      Real *matrix,
                                      Real mean = 0.0, Real stddev = 1.0 );

    static MERSENNETWISTER             GENERATOR;

	private:
    // Constructor
		MathUtil();

		// Destructor
		virtual ~MathUtil();

};

#endif
