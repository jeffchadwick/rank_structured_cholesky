//////////////////////////////////////////////////////////////////////
// CHOLMOD_Supernodal.cpp: Implementation
//
//////////////////////////////////////////////////////////////////////

#include "CHOLMOD_Supernodal.h"

#include "CHOLMOD_Environment.h"

#ifdef USE_MKL
#include <mkl_lapack.h>
#include <mkl_types.h>
#include <mkl_cblas.h>
#else
// Just use standard blas libraries
extern "C" {
#include <cblas.h>
#include <lapacke.h>
}
#endif

#include <set>

#include <util/IO.h>

#include <cholmod_blas.h>

using namespace std;

//////////////////////////////////////////////////////////////////////
// Static instantiator, which is handy so that we can
// make sure the factor satisfies certain conditions.
//////////////////////////////////////////////////////////////////////
CHOLMOD_Supernodal *CHOLMOD_Supernodal::Build( cholmod_common &common,
                                               cholmod_factor *factor )
{
  if ( !factor->is_super )
  {
    return NULL;
  }
  else
  {
    return new CHOLMOD_Supernodal( common, factor );
  }
}

//////////////////////////////////////////////////////////////////////
// Constructor:
// We need the factorization itself, and whatever cholmod
// workspace it is associated with.
//////////////////////////////////////////////////////////////////////
CHOLMOD_Supernodal::CHOLMOD_Supernodal( cholmod_common &common,
                                         cholmod_factor *factor )
  : _common( common ),
    _factor( factor ),
    _indexToSupernode( factor->n ),
    _p( NULL ),
    _i( NULL )
{
  buildIndexToSupernode();
  buildCompressedSymbolicFactor();
}

//////////////////////////////////////////////////////////////////////
// Destructor
//////////////////////////////////////////////////////////////////////
CHOLMOD_Supernodal::~CHOLMOD_Supernodal()
{
  delete[] _p;
  delete[] _i;
}

//////////////////////////////////////////////////////////////////////
// Writes some debugging output for the compressed factorization
//////////////////////////////////////////////////////////////////////
void CHOLMOD_Supernodal::debugPrintCompressed()
{
  int *factorNodeColumns = (int *)_factor->super;
  int *factorNodeIntPointers = (int *)_factor->pi;

  cout << "CHOLMOD_Supernodal::debugPrintCompressed" << endl;
  cout << SDUMP( _factor->nsuper ) << endl;
  for ( int i = 0; i < (int)_factor->nsuper; i++ )
  {
    int nrows = _p[ i+1 ] - _p[ i ];
    int start_idx = _p[ i ];

    int fullNumCols = factorNodeColumns[ i+1 ] - factorNodeColumns[ i ];
    int fullNumRows = factorNodeIntPointers[ i+1 ] - factorNodeIntPointers[ i ];

    cout << "Node " << i << " size: " << SDUMP( fullNumRows ) << SDUMP( fullNumCols ) << endl;

    cout << "Node " << i << " entries: " << endl;
    for ( int j = 0; j < nrows; j++ )
    {
      cout << "\t" << _i[ start_idx + j ] << endl;
    }
    cout << endl;
  }
}

//////////////////////////////////////////////////////////////////////
// Constructs the list of supernodes associated with
// the non-zero entries of a sparse vector.  The provided
// list should be at least as large as the most dense
// column of the input matrix.
// The output list may have duplicate supernodal indices,
// though this should not be an issue if we do DFS with
// markers.
//////////////////////////////////////////////////////////////////////
int CHOLMOD_Supernodal::getAssociatedSupernodes( const cholmod_sparse *B,
                                                  int k, /* column to consider */
                                                  int *nodes ) const
{
  const int   *Bp = (const int *)B->p;
  const int   *Bi = (const int *)B->i;

  int numNodes = 0;

  for ( int i = Bp[ k ]; i < Bp[ k+1 ]; i++ )
  {
    nodes[ numNodes ] = _indexToSupernode[ Bi[ i ] ];
    numNodes++;
  }

  return numNodes;
}

//////////////////////////////////////////////////////////////////////
// Solves the given sparse right hand side system using this
// factorization.  This only handles unpermuted lower triangular
// solves.
//    E is workspace - should have size at least _factor->maxesize
//    marker_list is just a copy of _p for marking purposes
//    k is the column of B to use as the RHS
//    xi is the set of supernode indices that will be
//      filled in to the solution
//    x is the dense solution vector.  Does not need
//      to be initialized.
//////////////////////////////////////////////////////////////////////
void CHOLMOD_Supernodal::spsolve( const cholmod_sparse *B,
                                  cholmod_dense *E,
                                  int *marker_list,
                                  int k, int *nodes, int *xi, double *x )
{
  // The logic here mostly follows cholmod_super_lsolve,
  // with the exception that we only iterate over a
  // subset of supernodes.
  double      *Lx, *Ex, *Bx;
  int         *Bp, *Bi;
  double       minus_one[ 2 ], one[ 2 ];
  int         *Lpi, *Lpx, *Ls, *Super;
  int          k1, k2, psi, psend, psx, nsrow, nscol, ii, s;
  int          nsrow2, n, ps2, j, i, nsuper;
  int          numNodes, top, nodeIndex;
  int          maxNumRows = 0;

  // Get inputs
  Ex = (double *)E->x;
  n = _factor->n;

  nsuper = _factor->nsuper;
  Lpi = (int *)_factor->pi;
  Lpx = (int *)_factor->px;
  Ls = (int *)_factor->s;
  Super = (int *)_factor->super;
  Lx = (double *)_factor->x;

  Bx = (double *)B->x;
  Bp = (int *)B->p;
  Bi = (int *)B->i;

  minus_one[ 0 ] = -1.0;
  minus_one[ 1 ] = 0.0;
  one[ 0 ] = 1.0;
  one[ 1 ] = 0.0;

  // Now, we need to figure out which supernodes will
  // be involved in the solution vector.
  CHOLMOD_Environment::CholmodWrapper G( this );

  _reachTimer.tick();
  // Get the list of nodes associated with this column
  numNodes = getAssociatedSupernodes( B, k, nodes );

  // Figure out which nodes will appear in the solution
  top = CHOLMOD_Environment::node_reach( &G, nodes, numNodes,
                                         marker_list, xi, NULL );
  _reachTimer.tock();

  _solveTimer.tick();

  // Start by cleaning the relevant portions of x, to make
  // sure we aren't polluted by a previous solve.
  for ( nodeIndex = top; nodeIndex < nsuper; nodeIndex++ )
  {
    s = xi[ nodeIndex ];

    k1 = Super[ s ];
    k2 = Super[ s+1 ];
    psi = Lpi[ s ];
    nscol = k2 - k1;

    for ( ii = 0; ii < nscol; ii++ )
    {
      x[ Ls[ psi + ii ] ] = 0.0;
    }
  }

  // Scatter the right hand side in to x
  for ( ii = Bp[ k ]; ii < Bp[ k+1 ]; ii++ )
  {
    x[ Bi[ ii ] ] = Bx[ ii ];
  }

  // Now we are ready to proceed with the solution process
  for ( nodeIndex = top; nodeIndex < nsuper; nodeIndex++ )
  {
    s = xi[ nodeIndex ];

    k1 = Super[ s ];
    k2 = Super[ s+1 ];
    psi = Lpi[ s ];
    psend = Lpi[ s+1 ];
    psx = Lpx[ s ];
    nsrow = psend - psi;
    nscol = k2 - k1;
    nsrow2 = nsrow - nscol;
    ps2 = psi + nscol;

    maxNumRows = max( nsrow2, maxNumRows );

    // L1 is nscol-by-nscol, lower triangular with non-unit diagonal.
    // L2 is nsrow2-by-nscol.  L1 and L2 have leading dimension of
    // nsrow.  x1 is nscol-by-nsrow, with leading dimension n.
    // E is nsrow2-by-1, with leading dimension nsrow2.

    // Gather input in to worksapce
    for ( ii = 0; ii < nsrow2; ii++ )
    {
      Ex[ ii ] = x[ Ls[ ps2 + ii ] ];
    }

    // Solve L1 * x1 (that is, x1 = L1 \ x1)
#if 0
    BLAS_dtrsv( "L", "N", "N",
                nscol,            // N: L1 is nscol x nscol
                Lx + psx, nsrow,  // A, LDA:  L1
                x + k1, 1 );      // x, incx = 1
#endif
    cblas_dtrsv( CblasColMajor, CblasLower, CblasNoTrans, CblasNonUnit,
                 nscol,            // N: L1 is nscol x nscol
                 Lx + psx, nsrow,  // A, LDA:  L1
                 x + k1, 1 );      // x, incx = 1

    // E = E - L2 * x1
#if 0
    BLAS_dgemv ("N",
                nsrow2, nscol,		    /* M, N:    L2 is nsrow2-by-nscol */
                minus_one,		    /* ALPHA:   -1 */
                Lx + ENTRY_SIZE*(psx + nscol),   /* A, LDA:  L2 */
                nsrow,
                Xx + ENTRY_SIZE*k1, 1,	    /* X, INCX: x1 */
                one,			    /* BETA:    1 */
                Ex, 1) ;		    /* Y, INCY: E */
#endif
    cblas_dgemv( CblasColMajor, CblasNoTrans,
                 nsrow2, nscol,       // M, N:    L2 is nsrow2 x nscol
                 minus_one[0],        // ALPHA:   -1
                 Lx + psx + nscol,    // A:       L2
                 nsrow,               // LDA:     nsrow
                 x + k1, 1,           // X, INCX: x1
                 one[0],              // BETA: 1
                 Ex, 1 );             // Y, INCY: E

    // Scatter E back to X
    for ( ii = 0; ii < nsrow2; ii++ )
    {
      x[ Ls[ ps2 + ii ] ] = Ex[ ii ];
    }
  }

  // Clean whatever part of the workspace we used
  memset( (void *)Ex, 0, maxNumRows * sizeof( double ) );

  _solveTimer.tock();
}

//////////////////////////////////////////////////////////////////////
// Builds a map of matrix indices to supernodes
//////////////////////////////////////////////////////////////////////
void CHOLMOD_Supernodal::buildIndexToSupernode()
{
  int *factorNodeColumns = (int *)_factor->super;
  int *factorNodeIntPointers = (int *)_factor->pi;
  int *factorMapToIndices = (int *)_factor->s;
  
  int  startIndex;
  int  ncols;

  for ( int i = 0; i < (int)_factor->nsuper; i++ )
  {
    ncols = factorNodeColumns[ i+1 ] - factorNodeColumns[ i ];

    startIndex = factorNodeIntPointers[ i ];

    for ( int j = 0; j < ncols; j++ )
    {
      // Assign the appropriate supernode to this index
      _indexToSupernode[ factorMapToIndices[ startIndex + j ] ] = i;
    }
  }
}

//////////////////////////////////////////////////////////////////////
// Compresses a indices in to supernodes and produces a
// symbolic "factorization" which encodes the graph connecting
// supernodes.
//////////////////////////////////////////////////////////////////////
void CHOLMOD_Supernodal::buildCompressedSymbolicFactor()
{
  // We know how many columns are in the symbolic factor - namely,
  // the number of supernodes
  int        nCol = _factor->nsuper;
  int        nzmax = 0;
  int        nColNode, nRowNode;
  int        startIndex;
  int        superNodeID;
  int        nColEntries;
  int        currentSize;
  int        oldSize;

  int       *factorNodeColumns = (int *)_factor->super;
  int       *factorNodeIntPointers = (int *)_factor->pi;
  int       *factorMapToIndices = (int *)_factor->s;

  // Need to maintain a set of supernodes encountered so far.
  // This isn't great, and implies M log M complexity, where
  // M is the number of supernodes.
  set<int>   nodeSet;

  _p = new int[ nCol + 1 ];

  // First pass: figure out how much space to allocate
  for ( int i = 0; i < nCol; i++ )
  {
    nColNode = factorNodeColumns[ i+1 ] - factorNodeColumns[ i ];
    nRowNode = factorNodeIntPointers[ i+1 ] - factorNodeIntPointers[ i ];
    startIndex = factorNodeIntPointers[ i ];

    // We know that this node is associated with itself,
    // so for now just get a count of additional nodes.
    nodeSet.insert( i );

    for ( int j = nColNode; j < nRowNode; j++ )
    {
      superNodeID = _indexToSupernode[ factorMapToIndices[ startIndex + j ] ];

      nodeSet.insert( superNodeID );
    }

    _p[ i ] = nodeSet.size();
    nzmax += _p[ i ];
    nodeSet.clear();
  }

  _p[ nCol ] = 0;

  // Initialize the index data
  _i = new int[ nzmax ];

  nzmax = 0;

  // Second pass: fill in the index data
  for ( int i = 0; i < nCol; i++ )
  {
    nColEntries = _p[ i ];

    nColNode = factorNodeColumns[ i+1 ] - factorNodeColumns[ i ];
    nRowNode = factorNodeIntPointers[ i+1 ] - factorNodeIntPointers[ i ];
    startIndex = factorNodeIntPointers[ i ];

    // We know that this node is associated with itself,
    // so for now just get a count of additional nodes.
    nodeSet.insert( i );

    for ( int j = nColNode; j < nRowNode; j++ )
    {
      superNodeID = _indexToSupernode[ factorMapToIndices[ startIndex + j ] ];

      nodeSet.insert( superNodeID );
    }

    for ( set<int>::iterator j = nodeSet.begin(); j != nodeSet.end(); j++ )
    {
      _i[ nzmax ] = *j;
      nzmax++;
    }

    nodeSet.clear();
  }

  // Third pass - fix the values stored in p so that
  // they are proper column pointers
  currentSize = 0;
  for ( int i = 0; i < nCol + 1; i++ )
  {
    oldSize = currentSize;
    currentSize += _p[ i ];
    _p[ i ] = oldSize;
  }
}
