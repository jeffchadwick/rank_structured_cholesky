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
