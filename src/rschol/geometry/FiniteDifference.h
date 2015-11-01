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
// FiniteDifference.h: Interface for the FiniteDifference class
//
//////////////////////////////////////////////////////////////////////

#ifndef FINITE_DIFFERENCE_H
#define FINITE_DIFFERENCE_H

#include <rschol/linearalgebra/MATRIX.h>
#include <rschol/linearalgebra/MATRIX3.h>
#include <rschol/linearalgebra/SPARSE_MATRIX.h>
#include <rschol/linearalgebra/VECTOR.h>

#include <rschol/util/Evaluator.h>

#include <rschol/SETTINGS.h>
#include <rschol/TYPES.h>

#include <rschol/geometry/BoundingBox.h>

#include <rschol/util/IO.h>

#define INDEX_3D_GRID(gridDivisions, idx) \
  idx[2] + gridDivisions[2] * ( idx[1] + gridDivisions[1] * idx[0] ) \

#define FOR_ALL_3D_POINTS(gridDivisions, idx) \
  for ( Tuple3i idx(0,0,0); idx[0] < gridDivisions[0]; idx[0]++ ) \
  for ( idx[1] = 0; idx[1] < gridDivisions[1]; idx[1]++ ) \
  for ( idx[2] = 0; idx[2] < gridDivisions[2]; idx[2]++ ) \

#define FOR_ALL_3D_POINTS_IN_RANGE(lowBound, highBound, idx) \
  for ( Tuple3i idx(lowBound); idx[0] <= highBound[0]; idx[0]++ ) \
  for ( idx[1] = lowBound[1]; idx[1] <= highBound[1]; idx[1]++ ) \
  for ( idx[2] = lowBound[2]; idx[2] <= highBound[2]; idx[2]++ ) \

#define GRID_SIZE(gridDivisions) \
  gridDivisions[0] * gridDivisions[1] * gridDivisions[2] \

//////////////////////////////////////////////////////////////////////
// FiniteDifference class
//
// Functions for generating finite difference grids, etc.
//////////////////////////////////////////////////////////////////////
class FiniteDifference {
	public:
    // Generates a finite difference Laplacian matrix with Dirichlet
    // boundary conditions, assuming a constant grid cell spacing of h
    static void generateLaplacian( SPARSE_MATRIX &L,
                                   Tuple3i gridDivisions,
                                   Real h );

    // Evaluates a function at all points of a finite difference grid
    static void evaluateGridEntries( const Tuple3i &gridDivisions, Real h,
                                     const ScalarEvaluator &f,
                                     VECTOR &x );

    // Models a finite difference grid point.  Useful for mapping
    // indices to grid indices, etc.
    class GridPoint {
      public:
        GridPoint( int index, const Tuple3i &gridDivisions,
                   Real h = 1.0, VEC3F origin = VEC3F( 0.0, 0.0, 0.0 ) )
          : _index( index )
        {
          // Figure out the x, y and z indices for the point
          int              x_idx, y_idx, z_idx;
          int              remainder;

          x_idx = index / ( gridDivisions[1] * gridDivisions[2] );

          remainder = index % ( gridDivisions[1] * gridDivisions[2] );

          y_idx = remainder / ( gridDivisions[2] );
          z_idx = remainder % gridDivisions[2];

          TRACE_ASSERT( x_idx >= 0 && x_idx < gridDivisions[0]
                     && y_idx >= 0 && y_idx < gridDivisions[1]
                     && z_idx >= 0 && z_idx < gridDivisions[2],
                     "Index out of range" );

          Tuple3i idx( x_idx, y_idx, z_idx );

          VEC3F x = assignPoint( idx, h );

          _bbox.setMinBound( x );
          _bbox.setMaxBound( x );
        }

        int index() const { return _index; }

        const BoundingBox &bbox() const { return _bbox; }

      private:
        int                _index;

        // For KD tree compatibility
        BoundingBox        _bbox;
    };

  private:
		FiniteDifference() {}

		// Destructor
		virtual ~FiniteDifference() {}

    // Checks to make sure an index is valid
    static bool checkFDIndex( Tuple3i gridDivisions, Tuple3i idx )
    {
      return ( idx[0] >= 0 && idx[1] >= 0 && idx[2] >= 0 &&
               idx[0] < gridDivisions[0] &&
               idx[1] < gridDivisions[1] &&
               idx[2] < gridDivisions[2] );
    }

	protected:

  private:
    static inline VEC3F assignPoint( const Tuple3i &idx, Real h )
    {
      VEC3F p;

      for ( int i = 0; i < 3; i++ )
      {
        p[ i ] = h * (Real)( idx[ i ] + 1 );
      }

      return p;
    }

	private:
    // Stencils for finite difference laplacian operator
    static const Real LAPLACIAN_STENCIL[3];
    static const int LAPLACIAN_STENCIL_INDEX[3];
    static const int LAPLACIAN_STENCIL_SZ;
};

#endif
