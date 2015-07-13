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
