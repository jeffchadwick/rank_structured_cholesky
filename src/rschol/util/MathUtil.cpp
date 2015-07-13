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

