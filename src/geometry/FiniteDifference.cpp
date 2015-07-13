//////////////////////////////////////////////////////////////////////
// FiniteDifference.cpp: Implementation
//
//////////////////////////////////////////////////////////////////////

#include "FiniteDifference.h"

//////////////////////////////////////////////////////////////////////
// Stencils for finite difference laplacian operator
//////////////////////////////////////////////////////////////////////
const Real FiniteDifference::LAPLACIAN_STENCIL[] = { -1.0, 2.0, -1.0 };
const int FiniteDifference::LAPLACIAN_STENCIL_INDEX[] = { -1, 0, 1 };
const int FiniteDifference::LAPLACIAN_STENCIL_SZ = 3;

//////////////////////////////////////////////////////////////////////
// Generates a finite difference Laplacian matrix with Dirichlet
// boundary conditions, assuming a constant grid cell spacing of h
//////////////////////////////////////////////////////////////////////
void FiniteDifference::generateLaplacian( SPARSE_MATRIX &L,
                                          Tuple3i gridDivisions,
                                          Real h )
{
  int sz = GRID_SIZE( gridDivisions );
  Real inv_h2 = 1.0 / ( h * h );

  L.resize( sz, sz );
  L.clear();

  FOR_ALL_3D_POINTS( gridDivisions, idx )
  {
    // Apply the Laplacian along each dimension
    for ( int d = 0; d < 3; d++ )
    {
      for ( int lap_idx = 0; lap_idx < LAPLACIAN_STENCIL_SZ; lap_idx++ )
      {
        Tuple3i entry_idx = idx;

        int row_idx = INDEX_3D_GRID( gridDivisions, idx );

        int offset = LAPLACIAN_STENCIL_INDEX[ lap_idx ];

        entry_idx[d] += offset;

        if ( !checkFDIndex( gridDivisions, entry_idx ) )
          continue;

        int col_idx = INDEX_3D_GRID( gridDivisions, entry_idx );

        L( row_idx, col_idx ) += LAPLACIAN_STENCIL[ lap_idx ] * inv_h2;
      }
    }
  }
}

// Evaluates a function at all points of a finite difference grid
void FiniteDifference::evaluateGridEntries( const Tuple3i &gridDivisions, Real h,
                                            const ScalarEvaluator &f,
                                            VECTOR &x )
{
  int sz = GRID_SIZE( gridDivisions );
  VEC3F p;

  int x_idx = 0;

  x.resizeAndWipe( sz );

  FOR_ALL_3D_POINTS( gridDivisions, idx )
  {
    p = assignPoint( idx, h );

    x( x_idx ) = f( p );

    x_idx++;
  }
}
