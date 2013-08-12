//////////////////////////////////////////////////////////////////////
// CHOLMOD_Environment.cpp: Implementation
//
//////////////////////////////////////////////////////////////////////

#include "CHOLMOD_Environment.h"

#include <stdlib.h>

#include <datastructure/ThreadSpecificData.h>

#include <util/IO.h>
#include <util/STLUtil.h>

#ifdef USE_OMP
#include <omp.h>
#endif

//////////////////////////////////////////////////////////////////////
// Constructor
//////////////////////////////////////////////////////////////////////
CHOLMOD_Environment::CHOLMOD_Environment( int ordering )
  : _matrix( NULL ),
    _matrixTriplet( NULL ),
    _factor( NULL ),
    _factorSparse( NULL ),
    _factorTriplet( NULL ),
    _factorTimer( "Factor timer" ),
    _solveTimer( "Solve timer" ),
    _testSolveTimer( "Test Solve Timer" )
{
  cholmod_start( &_common );

#if 0
  // Use nested dissection for now
  _common.nmethods = 1;
  _common.method[0].ordering = ordering;

#if 0
  _common.supernodal = CHOLMOD_SIMPLICIAL;
#endif
  _common.supernodal = CHOLMOD_SUPERNODAL;
#endif
#if 0
  _common.nmethods = 1;
  _common.method[0].ordering = CHOLMOD_NESDIS;
#endif

  _common.final_ll = 1;
  //_common.supernodal = CHOLMOD_SIMPLICIAL;

  _common.error_handler = CHOLMOD_Environment::handleError;
}

//////////////////////////////////////////////////////////////////////
// Destructor
//////////////////////////////////////////////////////////////////////
CHOLMOD_Environment::~CHOLMOD_Environment()
{
  freeSparse( _matrix );
  freeSparse( _factorSparse );

  freeTriplet( _matrixTriplet );
  freeTriplet( _factorTriplet );

  freeFactor( _factor );

  cholmod_finish( &_common );
}

//////////////////////////////////////////////////////////////////////
// Sets a new matrix to factorize
//////////////////////////////////////////////////////////////////////
void CHOLMOD_Environment::setMatrix( const SPARSE_MATRIX &M )
{
  clear();

  _factor = NULL;

  cout << SDUMP( M(0,0) ) << endl;

  allocateSparse( M, _matrix, _matrixTriplet, true );
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void CHOLMOD_Environment::setMatrix( const SPARSE_MATRIX::SparseColumnMatrix &M )
{
  clear();

  _factor = NULL;

  allocateSparseColumnCopy( M, _matrix, _matrixTriplet, true );
}

//////////////////////////////////////////////////////////////////////
// Sets a matrix from a stored sparse matrix file
//
// format = 0 (SPARSE_MATRIX object written to file)
// format = 1 (SPARSE_MATRIX::SparseColumnCopy object written to file)
//////////////////////////////////////////////////////////////////////
void CHOLMOD_Environment::setMatrix( const char *filename, int format )
{
  clear();

  _factor = NULL;

  allocateSparse( filename, format, _matrix, _matrixTriplet, true );
}

//////////////////////////////////////////////////////////////////////
// Factors the current matrix
//////////////////////////////////////////////////////////////////////
void CHOLMOD_Environment::computeCholeskyFactor( bool numeric, bool symbolic,
                                                 int *userPerm )
{
  if ( !_matrix )
    return;

#if 0
  // FIXME
  userPerm = new int[ _matrix->nrow ];
  for ( int i = 0; i < _matrix->nrow; i++ )
  {
    userPerm[ i ] = i;
  }
#endif

  if ( symbolic )
  {
    if ( userPerm )
    {
      _factor = cholmod_analyze_p( _matrix, userPerm, NULL, 0, &_common );
    }
    else
    {
      cout << "Analyzing" << endl;
      _factor = cholmod_analyze( _matrix, &_common );
    }
  }

  if ( numeric )
  {
    _factorTimer.tick();
    cholmod_factorize( _matrix, _factor, &_common );
    _factorTimer.tock();
  }

  cholmod_print_sparse( _matrix, "Cholesky Factor Input", &_common );
  cholmod_print_factor( _factor, "Cholesky Factor Output", &_common );

  if ( _common.method[ _common.selected ].ordering == CHOLMOD_AMD )
  {
    cout << "AMD used!!!" << endl;
  }
  else if ( _common.method[ _common.selected ].ordering == CHOLMOD_METIS )
  {
    cout << "METIS used!!!" << endl;
  }
  else
  {
    cout << "Some other ordering was used :(" << endl;
  }

  cout << SDUMP( _factor->nsuper ) << endl;
  cout << SDUMP( _factor->nzmax ) << endl;

  // FIXME
  long int totalNz = 0;
  int *super = (int *)_factor->super;
  int *pi = (int *)_factor->pi;

  for ( int super_idx = 0; super_idx < _factor->nsuper; super_idx++ )
  {
    int nCols = super[ super_idx + 1 ] - super[ super_idx ];
    int nRows = pi[ super_idx + 1 ] - pi[ super_idx ];

    totalNz += nCols * nRows;
    //totalNz -= nCols * ( nCols - 1 ) / 2;
  }
  cout << SDUMP( totalNz ) << endl;
  cout << SDUMP( _common.nrelax[ 0 ] ) << endl;
  cout << SDUMP( _common.nrelax[ 1 ] ) << endl;
  cout << SDUMP( _common.nrelax[ 2 ] ) << endl;
  cout << SDUMP( _common.zrelax[ 0 ] ) << endl;
  cout << SDUMP( _common.zrelax[ 1 ] ) << endl;
  cout << SDUMP( _common.zrelax[ 2 ] ) << endl;
#if 0
#endif

#if 0
  int *Lperm = (int *)_factor->Perm;

  for ( int i = 0; i < _matrix->nrow; i++ )
  {
    TRACE_ASSERT( i == Lperm[ i ], "Non identity permutation" );
    //printf( "  Permutation[ %d ] = %d\n", i, Lperm[ i ] );
  }
#endif
}

//////////////////////////////////////////////////////////////////////
// Puts the factor in the given sparse matrix
//////////////////////////////////////////////////////////////////////
void CHOLMOD_Environment::getFactor( SPARSE_MATRIX &G )
{
  if ( !_factor )
    return;

  freeSparse( _factorSparse );
  freeTriplet( _factorTriplet );

#if 0
  G.resize( _factor->n, _factor->n );
  G.clear();
#endif

  cholmod_factor *factorCopy = cholmod_copy_factor( _factor, &_common );

  _factorSparse = cholmod_factor_to_sparse( factorCopy, &_common );
  _factorTriplet = cholmod_sparse_to_triplet( _factorSparse, &_common );

  cholmod_free_factor( &factorCopy, &_common );

#if 0
  int *row_indices = (int *)_factorTriplet->i;
  int *col_indices = (int *)_factorTriplet->j;
  double *matrix_entries = (double *)_factorTriplet->x;

  for ( int i = 0; i < _factorTriplet->nnz; i++ )
  {
    int row_idx = row_indices[ i ];
    int col_idx = col_indices[ i ];

    double value = matrix_entries[ i ];

    G( row_idx, col_idx ) = value;
  }
#endif

  tripletToSparse( _factorTriplet, G );

  // If we have an LDL' factorization, convert it to a
  // standard Cholesky factor.
  if ( !_factor->is_ll )
  {
    VECTOR diagonals( _factor->n );

    for ( int i = 0; i < _factor->n; i++ )
    {
      diagonals(i) = sqrt( G(i, i) );
      G(i, i) = diagonals(i);
    }

    cout << "Converting LDL' factorization to Cholesky" << endl;

    for ( SPARSE_MATRIX::iterator i = G.begin(); i != G.end(); i++ )
    {
      int row_idx = i->first.first;
      int col_idx = i->first.second;

      if ( row_idx == col_idx )
        continue;

      i->second *= diagonals( col_idx );
    }
  }
  else
  {
    cout << "Factorization already in Cholesky format: no conversion needed" << endl;
  }
}

//////////////////////////////////////////////////////////////////////
    // Creates an empty factor
//////////////////////////////////////////////////////////////////////
void CHOLMOD_Environment::generateBlankFactor()
{
  if ( _factor )
  {
    cholmod_free_factor( &_factor, &_common );
  }

  _factor = cholmod_allocate_factor( _matrix->nrow, &_common );
}

//////////////////////////////////////////////////////////////////////
// Builds an ordered list of separators resulting from CHOLMOD'S
// nexted dissection routine
//////////////////////////////////////////////////////////////////////
void CHOLMOD_Environment::nestedDissection( SeparatorList &separators,
                                            IntArray &permutation,
                                            IntArray &separatorParents,
                                            IntArray &separatorLevels,
                                            IntArray &separatorMap )
{
  int                       *perm;
  int                       *cParent;
  int                       *cMember;

  IntArray                   inversePermutation;

  UF_long                    numComponents;

  int                        currentComponent = 0;

  TRACE_ASSERT( _matrix != NULL, "No matrix assigned" );

  perm = new int[ _matrix->nrow ];
  cParent = new int[ _matrix->nrow ];
  cMember = new int[ _matrix->nrow ];

  numComponents = cholmod_nested_dissection( _matrix, NULL, 0,
                                             perm, cParent, cMember,
                                             &_common );

  permutation.resize( _matrix->nrow );
  separatorMap.resize( _matrix->nrow );
  for ( int i = 0; i < _matrix->nrow; i++ )
  {
    permutation[ i ] = perm[ i ];

    separatorMap[ i ] = cMember[ i ];
  }
  invertIntArray( permutation, inversePermutation );

  separators.resize( numComponents );
  separatorParents.resize( numComponents );
  for ( int i = 0; i < numComponents; i++ )
  {
    separators[ i ]._size = 0;
    separators[ i ]._numChildren = 0;

    separatorParents[ i ] = cParent[ i ];
  }

  for ( int i = 0; i < numComponents; i++ )
  {
    int parent = cParent[ i ];
    if ( parent >= 0 )
    {
      separators.at( parent )._numChildren += 1;
    }
  }

  // Determine the level of each separator in the nested dissection
  // tree
  getSeparatorLevels( separatorParents, separatorLevels );

  // Check to see if this component list works the way I think it does
  for ( int i = 0; i < inversePermutation.size(); i++ )
  {
    int newIndex = permutation[ i ];
    int componentNumber = cMember[ newIndex ];

    //printf( "Node[ %d ] component number = %d\n", newIndex, componentNumber );

    if ( componentNumber == currentComponent + 1 )
    {
      currentComponent += 1;
    }
    else if ( componentNumber > currentComponent + 1
           || componentNumber < currentComponent )
    {
      cerr << "Bad component number" << endl;
      abort();
    }

    separators[ componentNumber ]._size += 1;
  }

  separators[ 0 ]._columnRange.first = 0;
  separators[ 0 ]._columnRange.second = separators[ 0 ]._size - 1;
  for ( int i = 1; i < numComponents; i++ )
  {
    separators[ i ]._columnRange.first
        = separators[ i - 1 ]._columnRange.second + 1;
    separators[ i ]._columnRange.second
        = separators[ i - 1 ]._columnRange.second + separators[ i ]._size;
  }

#if 0
  // Print out the separators
  for ( int i = 0; i < numComponents; i++ )
  {
    printf( "separator[ %08d ]: size = %08d, range = [ %08d, %08d ],"
            " children = [ %d ]\n",
            i, separators[ i ]._size, separators[ i ]._columnRange.first,
            separators[ i ]._columnRange.second,
            separators[ i ]._numChildren );
  }
#endif

  delete[] perm;
  delete[] cParent;
  delete[] cMember;
}

//////////////////////////////////////////////////////////////////////
// Returns a wrapper for a CHOLMOD supernodal factorization
//////////////////////////////////////////////////////////////////////
bool CHOLMOD_Environment::getSupernodalFactor(
                                  CHOLMOD_Environment::FactorWrapper &factor )
{
  if ( !_factor || !_factor->is_ll || !_factor->is_super )
    return false;

  if ( _factor->dtype != CHOLMOD_DOUBLE )
  {
    cerr << "Single precision factors not supported!" << endl;
    abort();
  }

  factor._n = _factor->n;
  factor._nsuper = _factor->nsuper;
  factor._ssize = _factor->ssize;
  factor._xsize = _factor->xsize;
  factor._maxcsize = _factor->maxcsize;
  factor._maxesize = _factor->maxesize;

  // Construct copies of the fields in _factor
  factor._perm = (int *)malloc( factor._n * sizeof( int ) );
  factor._x = (double *)malloc( factor._xsize * sizeof( double ) );
  factor._super = (int *)malloc( (factor._nsuper + 1) * sizeof( int ) );
  factor._pi = (int *)malloc( (factor._nsuper + 1) * sizeof( int ) );
  factor._px = (int *)malloc( (factor._nsuper + 1) * sizeof( int ) );
  factor._s = (int *)malloc( factor._ssize * sizeof( int ) );

  memcpy( (void *)factor._perm, (void *)_factor->Perm,
          factor._n * sizeof( int ) );

  if ( _factor->xtype == CHOLMOD_REAL )
  {
    memcpy( (void *)factor._x, (void *)_factor->x,
            factor._xsize * sizeof( double ) );
  }
  else
  {
    // Symbolic factor only
    factor._x = NULL;
  }

  memcpy( (void *)factor._super, (void *)_factor->super,
          (factor._nsuper + 1) * sizeof( int ) );
  memcpy( (void *)factor._pi, (void *)_factor->pi,
          (factor._nsuper + 1) * sizeof( int ) );
  memcpy( (void *)factor._px, (void *)_factor->px,
          (factor._nsuper + 1) * sizeof( int ) );
  memcpy( (void *)factor._s, (void *)_factor->s,
          factor._ssize * sizeof( int ) );

  return true;
}

//////////////////////////////////////////////////////////////////////
// Solves the factor system using the given right hand side.
//////////////////////////////////////////////////////////////////////
void CHOLMOD_Environment::solveFactorSystem( const SPARSE_MATRIX &rhs,
                                             SPARSE_MATRIX &solution,
                                             bool permuted )
{
  // Put the input matrix in CHOLMOD's format
  cholmod_sparse    *rhs_cholmod;
  cholmod_triplet   *rhs_triplet;
  cholmod_sparse    *solution_cholmod;
  cholmod_sparse    *solution_permuted;
  cholmod_triplet   *solution_triplet;

  allocateSparse( rhs, rhs_cholmod, rhs_triplet );

  cholmod_free_triplet( &rhs_triplet, &_common );

  cout << "Running CHOLMOD spsolve" << endl;
  if ( _factor->xtype == CHOLMOD_PATTERN )
  {
    cout << "Factor only stores pattern!" << endl;
  }
  else
  {
    cout << "Factor stores numerical factor" << endl;
  }
  cout << SDUMP( _factor->is_ll ) << "    " << SDUMP( _factor->is_super ) << endl;
  _solveTimer.tick();
  solution_cholmod = cholmod_spsolve( CHOLMOD_L, _factor,
                                      rhs_cholmod, &_common );

  if ( permuted )
  {
    // Permute the solution
    solution_permuted = cholmod_spsolve( CHOLMOD_P, _factor,
                                         solution_cholmod, &_common );

    // Throw away the unpermuted solution
    cholmod_free_sparse( &solution_cholmod, &_common );

    solution_cholmod = solution_permuted;
  }

  _solveTimer.tock();
  cout << "Done" << endl << endl;

  // FIXME
  int *pinv = NULL;
#if 0
  if ( _factor->Perm )
  {
    pinv = new int[ _factor->n ];

    int *perm = (int *)_factor->Perm;

    for ( int i = 0; i < _factor->n; i++ )
    {
      pinv[ perm[ i ] ] = i;

      cout << "perm[ " << i << " ] = " << perm[ i ] << endl;
    }
  }
#endif

#if 0
  if ( !_factor->is_ll )
  {
    cout << "FACTOR NOT CHOLESKY!" << endl;
  }

  if ( !rhs_cholmod->packed )
  {
    cout << "RHS is NOT packed!" << endl;
  }

  if ( !solution_cholmod->sorted )
  {
    cout << "Solution columns not sorted!" << endl;
  }
#endif

  // Try looking at the supernodal factorization here...
  CHOLMOD_Supernodal *super = CHOLMOD_Supernodal::Build( _common, _factor );

  if ( super )
  {
#if 0
    cout << "Compressed supernodal columns " << endl;
    super->debugPrintCompressed();
    cout << endl << endl;
#endif
  }
  else
  {
    cout << "Supernodal construct failed to build!" << endl;
  }

  ThreadSpecificData<int *> xi_set( NULL );
  ThreadSpecificData<double *> x_set( NULL );
  ThreadSpecificData<int *> marker_list_set( NULL );
  ThreadSpecificData<int *> nodes_set( NULL );
  ThreadSpecificData<cholmod_dense *> E_set( NULL );

  const int *LCp = super->compressedP();

#if 0
  int *xi = new int[ 2 * super->compressedN() ];
  double *x = new double[ rhs.rows() ];
  int *marker_list = new int[ super->compressedN() + 1 ];
  int *nodes = new int[ rhs.rows() ];
  cholmod_dense *E = super->buildSolverWorkspace( 1 );

  memset( (void *)x, 0, rhs.rows() * sizeof( double ) );

  for ( int i = 0; i < super->compressedN() + 1; i++ )
  {
    marker_list[ i ] = LCp[ i ];
  }
#endif

  cout << "Running custom sparse solve" << endl;
  _testSolveTimer.tick();
//#pragma omp parallel for schedule(static) default(shared)
  for ( int i = 0; i < rhs.cols(); i++ )
  {
    int *xi = xi_set.get();
    double *x = x_set.get();
    int *marker_list = marker_list_set.get();
    int *nodes = nodes_set.get();
    cholmod_dense *E = E_set.get();

    if ( !xi )
    {
      xi = new int[ 2 * super->compressedN() ];
    }

    if ( !x )
    {
      x = new double[ rhs.rows() ];
      memset( (void *)x, 0, rhs.rows() * sizeof( double ) );
    }

    if ( !marker_list )
    {
      marker_list = new int[ super->compressedN() + 1 ];

      for ( int i = 0; i < super->compressedN() + 1; i++ )
      {
        marker_list[ i ] = LCp[ i ];
      }
    }

    if ( !nodes )
    {
      nodes = new int[ rhs.rows() ];
    }

    if ( !E )
    {
      E = super->buildSolverWorkspace( 1 );
    }

    super->spsolve( rhs_cholmod, E, marker_list, i, nodes, xi, x );
  }
  _testSolveTimer.tock();
  cout << "Done" << endl << endl;

  cout << SDUMP( super->_reachTimer.getTotalSecs() ) << endl;
  cout << SDUMP( super->_solveTimer.getTotalSecs() ) << endl;

  int *Sp = (int *)solution_cholmod->p;
  int *Si = (int *)solution_cholmod->i;
  double *Sx = (double *)solution_cholmod->x;

#if 0
  if ( solution_cholmod->packed )
  {
    for ( int i = Sp[29]; i < Sp[30]; i++ )
    {
      //cout << "Packed Cholmod solution[ " << Si[ i ] << " ] = " << Sx[ Si[ i ] ];
      cout << "Cholmod solution[ " << Si[ i ] << " ] = " << Sx[ i ];
      cout << "\t";
      cout << "\ts[ " << Si[ i ] << " ] = " << x[ Si[ i ] ] << endl;
    }
  }
  else
  {
  }
#endif

#if 0
  cout << SDUMP( _factor->n ) << endl;
  cout << SDUMP( rhs.rows() ) << "   " << SDUMP( rhs.cols() ) << endl;
  cholmod_factor *factorCopy = cholmod_copy_factor( _factor, &_common );
  cholmod_sparse *G_temp = cholmod_factor_to_sparse( factorCopy, &_common );
  CholmodWrapper G( G_temp );
  //CholmodWrapper G( _factor );
  cholmod_free_factor( &factorCopy, &_common );

  int *xi = new int[ 2 * rhs.rows() ];
  double *x = new double[ rhs.rows() ];
  int *marker_list = new int[ _factor->n + 1 ];
  //int *Lp = (int *)_factor->p;
  int *Lp = (int *)G_temp->p;
  for ( int i = 0; i < _factor->n + 1; i++ )
  {
    marker_list[ i ] = Lp[ i ];
  }

  _testSolveTimer.tick();
  int top;
  for ( int i = 0; i < rhs.cols(); i++ )
  {
    top = spsolve( &G, rhs_cholmod, marker_list, i, xi, x, pinv, 1 );
  }
  _testSolveTimer.tock();

  delete[] xi;
  delete[] x;
  delete[] marker_list;
#endif

#if 0
  for ( int i = top; i < _factor->n; i++ )
  {
    cout << "Solution[ " << xi[ i ] << " ] = " << x[ xi[ i ] ] << endl;
  }
  cout << endl;

  if ( Sp[1] - Sp[0] != _factor->n - top )
  {
    cout << "Solution size mismatch!" << endl;
  }

  if ( _factor->Perm )
  {
    cout << "Factor has a permutation" << endl;
  }

  double *Gp = (double *)_factor->x;

  cout << SDUMP( Gp[0] ) << endl;

  if ( solution_cholmod->packed )
  {
    int p = top;
    for ( int i = Sp[0]; i < Sp[1]; i++ )
    {
      //cout << "Packed Cholmod solution[ " << Si[ i ] << " ] = " << Sx[ Si[ i ] ];
      cout << "Packed Cholmod solution[ " << Si[ i ] << " ] = " << Sx[ i ];
      cout << "\t";
      cout << "\t[ " << xi[ p ] << " ] = " << x[ xi[ p ] ] << endl;
      p++;
    }
  }
  else
  {
  }
#endif

  // END FIXME

  cholmod_free_sparse( &rhs_cholmod, &_common );

  solution_triplet = cholmod_sparse_to_triplet( solution_cholmod, &_common );

  cholmod_free_sparse( &solution_cholmod, &_common );

  tripletToSparse( solution_triplet, solution );
}

//////////////////////////////////////////////////////////////////////
// Solves the matrix system (after factorization).
//////////////////////////////////////////////////////////////////////
void CHOLMOD_Environment::solveSystem( const VECTOR &rhs, VECTOR &solution )
{
  cholmod_dense           *rhsDense;
  cholmod_dense           *solutionDense;
  double                  *x;

  if ( !_factor ) return;

  solution.resizeAndWipe( rhs.size() );
  
  rhsDense = cholmod_allocate_dense( rhs.size(), 1, rhs.size(),
                                     CHOLMOD_REAL, &_common );

  x = (double *)rhsDense->x;

  memcpy( x, rhs.data(), rhs.size() * sizeof( double ) );

  solutionDense = cholmod_solve( CHOLMOD_A, _factor, rhsDense, &_common );

  x = (double *)solutionDense->x;

  memcpy( solution.data(), x, rhs.size() * sizeof( double ) );

  cholmod_free_dense( &rhsDense, &_common );
  cholmod_free_dense( &solutionDense, &_common );
}

//////////////////////////////////////////////////////////////////////
// Multiply the two given sparse matrices together.  The only
// reason this isn't static is because we need something to
// provide a cholmod_common environment.
// C = A * B (or, optionally use transposes)
//////////////////////////////////////////////////////////////////////
void CHOLMOD_Environment::sparseMultiply(
                                      const SPARSE_MATRIX &A, bool transposeA,
                                      const SPARSE_MATRIX &B, bool transposeB,
                                      SPARSE_MATRIX &C )
{
  cholmod_sparse          *sparseA;
  cholmod_sparse          *sparseB;
  cholmod_sparse          *sparseC;
  cholmod_triplet         *tripletA;
  cholmod_triplet         *tripletB;
  cholmod_triplet         *tripletC;

  allocateSparse( A, sparseA, tripletA, false /* not symmetric */, transposeA );
  allocateSparse( B, sparseB, tripletB, false /* not symmetric */, transposeB );

  // Don't need the triplet for anything
  cholmod_free_triplet( &tripletA, &_common );
  cholmod_free_triplet( &tripletB, &_common );

  sparseC = cholmod_ssmult( sparseA, sparseB,
                            0 /* not symmetric */,
                            1 /* calculate values, not just pattern */,
                            1 /* sort columns */,
                            &_common );

  cholmod_free_sparse( &sparseA, &_common );
  cholmod_free_sparse( &sparseB, &_common );

  tripletC = cholmod_sparse_to_triplet( sparseC, &_common );

  cholmod_free_sparse( &sparseC, &_common );

  // Convert this in to our object-oriented sparse matrix format.
  C.clearFull();
  tripletToSparse( tripletC, C );

  cholmod_free_triplet( &tripletC, &_common );
}

//////////////////////////////////////////////////////////////////////
// Writes the permutation used by the system to the given file
//////////////////////////////////////////////////////////////////////
void CHOLMOD_Environment::writePermutation( const char *filename )
{
  TRACE_ASSERT( _factor && _factor->Perm );

  // Copy to a vector
  int                       *perm = (int *)_factor->Perm;
  IntArray                   permutation( _factor->n );

  for ( int i = 0; i < _factor->n; i++ )
  {
    permutation[ i ] = perm[ i ];
  }

  writeVector( permutation, filename );

  cout << SDUMP( _common.method[ _common.selected ].prune_dense ) << endl;
  cout << SDUMP( _common.method[ _common.selected ].aggressive ) << endl;
  cout << SDUMP( _common.method[ _common.selected ].fl ) << endl;
  cout << SDUMP( _common.method[ _common.selected ].lnz ) << endl;
  cout << SDUMP( _common.postorder ) << endl;
}

//////////////////////////////////////////////////////////////////////
// Extracts and AMD ordering for the current system
//////////////////////////////////////////////////////////////////////
void CHOLMOD_Environment::approximateMinimumDegree( IntArray &permutation )
{
  int                        *perm = new int[ _matrix->nrow ];

  cholmod_amd( _matrix, NULL, 0 /* empty fset */, perm, &_common );

  permutation.resize( _matrix->nrow );

  for ( int i = 0; i < _matrix->nrow; i++ )
  {
    permutation[ i ] = perm[ i ];
  }

  delete[] perm;
}

//////////////////////////////////////////////////////////////////////
// Mark and flip routines for the cs_* function implementations
// copied from [Davis, 2006]
//////////////////////////////////////////////////////////////////////
#define CS_FLIP(i) ( -(i) - 2 )
#define CS_UNFLIP(i) ( ((i) < 0) ? CS_FLIP(i) : (i) )
#define CS_MARKED(w,j) ( w[ j ] < 0 )
#define CS_MARK(w,j) { w[j] = CS_FLIP( w[j] ); }

//////////////////////////////////////////////////////////////////////
// Solves a sparse triangular system with a sparse right hand side.
// xi is assumed to have size 2 * n, and x to have size n.
// This is based entirely on cs_spsolve from [Davis, 2006], pg. 34
//////////////////////////////////////////////////////////////////////
int CHOLMOD_Environment::spsolve( const CholmodWrapper *G,
                                  const cholmod_sparse *B,
                                  int *marker_list,
                                  int k, // The column of B
                                  int *xi,
                                  double *x,
                                  const int *pinv,
                                  int lo )
{
  int j, J, p, q, px, top, n;
  const int *Gp, *Gi, *Bp, *Bi;
  const double *Gx, *Bx;

  Gp = (const int *)G->p;
  Gi = (const int *)G->i;
  Gx = (const double *)G->x;
  n = G->n;

  Bp = (const int *)B->p;
  Bi = (const int *)B->i;
  Bx = (const double *)B->x;

  // x[top...n-1] = Reach(B(:,k))
  top = node_reach( G, B, marker_list, k, xi, pinv );

  // Clear x
  for ( p = top; p < n; p++ )
  {
    x[ xi[ p ] ] = 0.0;
  }

  // Scatter B
  for ( p = Bp[ k ]; p < Bp[ k+1 ]; p++ )
  {
    x[ Bi[ p ] ] = Bx[ p ];
  }

  for ( px = top; px < n; px++ )
  {
    j = xi[ px ];                                   // x(j) is nonzero
    J = pinv ? ( pinv[ j ] ) : j;                   // j maps to col J of G
    if ( J < 0 ) continue;                          // Column J is empty
    x[j] /= Gx[ lo ? (Gp[ J ]) : (Gp[ J+1 ] - 1) ]; // x(j) /= G(j,j)
    p = lo ? (Gp[ J ] + 1) : (Gp[ J ] );            // lo: L(j,j) 1st entry
    q = lo ? (Gp[ J+1 ]) : ( Gp[ J+1 ] - 1);        // up: U(j,j) last entry

    for ( ; p < q; p++ )
    {
      x[ Gi[ p ] ] -= Gx[ p ] * x[ j ];             // x(i) -= G(i,j) * x(j)
    }
  }

  return top;                                       // Return top of stack
}

//////////////////////////////////////////////////////////////////////
// Computes the reach of a set of nodes, defined to be the
// non-zero nodes in a given column of a sparse matrix
// This is based entirely on the cs_reach algorithm
// from Timothy Davis' "Direct Methods for Sparse Linear Systems"
// (page 33)
//////////////////////////////////////////////////////////////////////
int CHOLMOD_Environment::node_reach( const CholmodWrapper *G,
                                     const cholmod_sparse *B,
                                     int *marker_list,
                                     int k, // The column of B
                                     int *xi,
                                     const int *pinv )
{
  int p, n, top;
  const int *Bp, *Bi, *Gp;

  n = G->n;
  Bp = (const int*)B->p;
  Bi = (const int*)B->i;
  Gp = (const int*)G->p;

  top = n;

  for ( p = Bp[ k ]; p < Bp[ k+1 ]; p++ )
  {
    // Start a DFS at unmarked node
    if ( !CS_MARKED( marker_list, Bi[ p ] ) )
    {
      top = dfs_node_reach( Bi[ p ], G, marker_list, top, xi, xi+n, pinv );
    }
  }

  // Restore the marker list
  for ( p = top; p < n; p++ )
    CS_MARK( marker_list, xi[ p ] );

  return top;
}

//////////////////////////////////////////////////////////////////////
// Similar to the above, but with an explicitly provided
// list of starting nodes
//////////////////////////////////////////////////////////////////////
int CHOLMOD_Environment::node_reach( const CholmodWrapper *G,
                                     int *nodes, int numNodes,
                                     int *marker_list,
                                     int *xi,
                                     const int *pinv )
{
  int p, n, top;
  const int *Bp, *Bi, *Gp;

  n = G->n;
  Gp = (const int*)G->p;

  top = n;

  for ( p = 0; p < numNodes; p++ )
  {
    // Start a DFS at unmarked node
    if ( !CS_MARKED( marker_list, nodes[ p ] ) )
    {
      top = dfs_node_reach( nodes[ p ], G, marker_list, top, xi, xi+n, pinv );
    }
  }

  // Restore the marker list
  for ( p = top; p < n; p++ )
    CS_MARK( marker_list, xi[ p ] );

  return top;
}

//////////////////////////////////////////////////////////////////////
// Cleanup
//////////////////////////////////////////////////////////////////////
void CHOLMOD_Environment::clear()
{
  freeSparse( _matrix );
  freeTriplet( _matrixTriplet );
  freeFactor( _factor );
  freeSparse( _factorSparse );
  freeTriplet( _factorTriplet );
}

//////////////////////////////////////////////////////////////////////
// Converts our sparse matrix format to CHOLMOD's
//////////////////////////////////////////////////////////////////////
void CHOLMOD_Environment::allocateSparse( const SPARSE_MATRIX &M,
                                          cholmod_sparse * &A,
                                          cholmod_triplet * &T,
                                          bool symmetric,
                                          bool transpose )
{
  // This assumes that the diagonal has non-zeros on it
  int nnz;
  int idx = 0;

  if ( symmetric )
  {
    nnz = ( M.size() - M.cols() ) / 2 + M.cols();
  }
  else
  {
    nnz = M.size();
  }

  T = cholmod_allocate_triplet( transpose ? M.cols() : M.rows(),
                                transpose ? M.rows() : M.cols(),
                                nnz,
                                symmetric ? -1 : 0 /* symmetric, lower tri */,
                                CHOLMOD_REAL, &_common );

  int *row_indices = (int *)T->i;
  int *col_indices = (int *)T->j;
  double *matrix_entries = (double *)T->x;

  // Put the matrix information in to triplet form.  Pay attention
  // to only lower triangular entries.
  for ( SPARSE_MATRIX::const_iterator i = M.begin(); i != M.end(); i++ )
  {
    int row_idx = i->first.first;
    int col_idx = i->first.second;

    if ( row_idx < col_idx && symmetric )
      continue;

    row_indices[ idx ] = transpose ? col_idx : row_idx;
    col_indices[ idx ] = transpose ? row_idx : col_idx;
    matrix_entries[ idx ] = i->second;

    idx++;
  }

  T->nnz = nnz;

  A = cholmod_triplet_to_sparse( T, 0, &_common );
}

//////////////////////////////////////////////////////////////////////
// Reads in a sparse matrix from a file and puts it in CHOLMOD's
// internal format
//
// format = 0 (SPARSE_MATRIX object written to file)
// format = 1 (SPARSE_MATRIX::SparseColumnCopy object written to file)
//////////////////////////////////////////////////////////////////////
void CHOLMOD_Environment::allocateSparse( const char *filename,
                                          int format,
                                          cholmod_sparse * &A,
                                          cholmod_triplet * &T,
                                          bool symmetric,
                                          bool transpose )
{
  if ( format == 0 )
  {
  }
  else if ( format == 1 )
  {
    SPARSE_MATRIX::SparseColumnMatrix S;

    SPARSE_MATRIX::readFromBinary( S, filename );

    allocateSparseColumnCopy( S, A, T, symmetric, transpose );
  }
  else
  {
    cerr << "** WARNING ** Invalid format" << endl;
  }
}

//////////////////////////////////////////////////////////////////////
// Different versions of the above for different formats
//
// Assumes matrix written by SPARSE_MATRIX object
//////////////////////////////////////////////////////////////////////
void CHOLMOD_Environment::allocateSparseStandard( const char *filename,
                                                  cholmod_sparse * &A,
                                                  cholmod_triplet * &T,
                                                  bool symmetric,
                                                  bool transpose )
{
  TRACE_ASSERT( NULL, "Not implemented" );
}

//////////////////////////////////////////////////////////////////////
// Assumes matrix written by SparseColumnCopy object
//////////////////////////////////////////////////////////////////////
void CHOLMOD_Environment::allocateSparseColumnCopy(
                                    const SPARSE_MATRIX::SparseColumnMatrix &S,
                                    cholmod_sparse * &A,
                                    cholmod_triplet * &T,
                                    bool symmetric,
                                    bool transpose )
{
  // Allocate triplet version
  T = cholmod_allocate_triplet( transpose ? S._ncol : S._nrow,
                                transpose ? S._nrow : S._ncol,
                                S._nzmax,
                                symmetric ? 1 : 0 /* symmetric, lower tri */,
                                CHOLMOD_REAL, &_common );

  int *row_indices = (int *)T->i;
  int *col_indices = (int *)T->j;
  double *matrix_entries = (double *)T->x;
  long int idx = 0;

  size_t nzTotal = 0;

  for ( int col_idx = 0; col_idx < S._ncol; col_idx++ )
  {
    for ( int row_ptr = S._p[ col_idx ]; row_ptr < S._p[ col_idx + 1 ];
          row_ptr++ )
    {
      int row_idx = S._i[ row_ptr ];

      if ( symmetric && row_idx < col_idx )
      {
        continue;
      }

      // FIXME:
      if ( row_idx < 0 || row_idx >= S._nrow ) {
        cout << SDUMP( col_idx ) << SDUMP( row_ptr ) << endl;
        cout << SDUMP( S._p[ col_idx ] ) << SDUMP( S._p[ col_idx + 1 ] ) << endl;
        cout << SDUMP( row_idx ) << SDUMP( S._nrow ) << endl;
      }
      TRACE_ASSERT( row_idx >= 0 && row_idx < S._nrow, "Row out of range" );

      row_indices[ idx ] = transpose ? col_idx : row_idx;
      col_indices[ idx ] = transpose ? row_idx : col_idx;
      matrix_entries[ idx ] = S._x[ row_ptr ];

      idx++;

      nzTotal++;
    }
  }

  //T->nnz = S._nzmax;
  T->nnz = nzTotal;

  printf( "Found %zd non-zeros in matrix\n", nzTotal );

#if 0
  // WAY too slow
  printf( "Computing sparse fill-in\n" );
  nzTotal = SPARSE_MATRIX::computeFillIn( S );
  printf( "Unordered factor has %lld non-zeros\n\n", nzTotal );
#endif

  A = cholmod_triplet_to_sparse( T, 0, &_common );
}

//////////////////////////////////////////////////////////////////////
// Cleans up a sparse matrix
//////////////////////////////////////////////////////////////////////
void CHOLMOD_Environment::freeSparse( cholmod_sparse * &A )
{
  if ( A )
  {
    cholmod_free_sparse( &A, &_common );
  }
}

//////////////////////////////////////////////////////////////////////
// Cleans up a sparse matrix in triplet form
//////////////////////////////////////////////////////////////////////
void CHOLMOD_Environment::freeTriplet( cholmod_triplet * &T )
{
  if ( T )
  {
    cholmod_free_triplet( &T, &_common );
  }
}

//////////////////////////////////////////////////////////////////////
// Cleans up a sparse matrix factor
//////////////////////////////////////////////////////////////////////
void CHOLMOD_Environment::freeFactor( cholmod_factor * &G )
{
  if ( G )
  {
    cholmod_free_factor( &G, &_common );
  }
}

//////////////////////////////////////////////////////////////////////
// Given a nested dissection ordering, determines the level of
// each separator in the nested dissection tree
//////////////////////////////////////////////////////////////////////
void CHOLMOD_Environment::getSeparatorLevels(
                                    const IntArray &separatorParents,
                                    IntArray &separatorLevels )
{
  separatorLevels.clear();
  separatorLevels.resize( separatorParents.size(), -1 );

  for ( int sep_idx = 0; sep_idx < separatorParents.size(); sep_idx++ )
  {
    getSeparatorLevels( sep_idx, separatorParents, separatorLevels );
  }
}

//////////////////////////////////////////////////////////////////////
// Recursive helper function for the above
//////////////////////////////////////////////////////////////////////
void CHOLMOD_Environment::getSeparatorLevels(
                                    int sep_idx,
                                    const IntArray &separatorParents,
                                    IntArray &separatorLevels )
{
  if ( separatorLevels.at( sep_idx ) >= 0 )
  {
    return;
  }

  // If this has a parent set its level to 1 + it's parent's level.
  // Otherwise, it is at level 0
  if ( separatorParents.at( sep_idx ) < 0 )
  {
    separatorLevels.at( sep_idx ) = 0;
  }
  else
  {
    getSeparatorLevels( separatorParents[ sep_idx ],
                        separatorParents, separatorLevels );

    separatorLevels.at( sep_idx )
      = 1 + separatorLevels.at( separatorParents.at( sep_idx ) );
  }
}

//////////////////////////////////////////////////////////////////////
// Error handler for CHOLMOD
//////////////////////////////////////////////////////////////////////
void CHOLMOD_Environment::handleError( int status, const char *file,
                                       int line, const char *message )
{
  printf ("cholmod error: file: %s line: %d status: %d: %s\n",
          file, line, status, message) ;
  if (status < 0)
  {
    exit (0) ;
  }
}

//////////////////////////////////////////////////////////////////////
// Converts the CHOLMOD triplet format to our SPARSE_MATRIX
// format
//////////////////////////////////////////////////////////////////////
void CHOLMOD_Environment::tripletToSparse( const cholmod_triplet *T,
                                           SPARSE_MATRIX &A )
{
  A.resize( T->nrow, T->ncol );
  A.clear();

  int *row_indices = (int *)T->i;
  int *col_indices = (int *)T->j;
  double *matrix_entries = (double *)T->x;

  for ( int i = 0; i < T->nnz; i++ )
  {
    int row_idx = row_indices[ i ];
    int col_idx = col_indices[ i ];

    double value = matrix_entries[ i ];

    A( row_idx, col_idx ) = value;
  }
}

//////////////////////////////////////////////////////////////////////
// Computes the reach of a node in the graph of the given
// factor.  This is based entirely on the cs_dfs algorithm
// from Timothy Davis' "Direct Methods for Sparse Linear Systems"
//////////////////////////////////////////////////////////////////////
int CHOLMOD_Environment::dfs_node_reach( int j,
                                         const CholmodWrapper *G,
                                         int *marker_list,
                                         int top,
                                         int *xi,
                                         int *pstack,
                                         const int *pinv )
{
  int i, p, p2, done, jnew, head = 0;
  const int *Gp, *Gi;
  Gp = (const int *)G->p;
  Gi = (const int *)G->i;

  // Initialize the recursion stack
  xi[ 0 ] = j;

  while ( head >= 0 )
  {
    j = xi[ head ];
    jnew = pinv ? ( pinv[ j ] ) : j;

    if ( !CS_MARKED( marker_list, j ) )
    {
      CS_MARK( marker_list, j );
      pstack[ head ] = ( jnew < 0 ) ? 0 : CS_UNFLIP( marker_list[ jnew ] );
    }

    // Node j done if no unvisited neighbours
    done = 1;
    p2 = ( jnew < 0 ) ? 0 : CS_UNFLIP( marker_list[ jnew + 1 ] );

    // Examine all neighbours of j
    for ( p = pstack[ head ]; p < p2; p++ )
    {
      i = Gi[ p ];                  // Consider neighbour node i

      // Skip visited node i
      if ( CS_MARKED( marker_list, i ) )
        continue;

      pstack[ head ] = p;           // Pause depth-first search of node j
      xi [ ++head ] = i;            // Start DFS at node i
      done = 0;                     // Node j is not done
      break;                        // Break to start DFS(i)
    }

    if ( done )
    {
      head--;                       // Remove j from the recursion stack
      xi[ --top ] = j;              // and place in the output stack
    }
  }

  return top;
}

//////////////////////////////////////////////////////////////////////
// Long format definitions.  Yes, this is an ugly way of doing things
//////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////
// Constructor
//////////////////////////////////////////////////////////////////////
CHOLMOD_Environment_Long::CHOLMOD_Environment_Long( int ordering )
  : _matrix( NULL ),
    _matrixTriplet( NULL ),
    _factor( NULL ),
    _factorSparse( NULL ),
    _factorTriplet( NULL ),
    _factorTimer( "Factor timer" ),
    _solveTimer( "Solve timer" ),
    _testSolveTimer( "Test Solve Timer" )
{
  cholmod_l_start( &_common );

#if 0
  // Use nested dissection for now
  _common.nmethods = 1;
  _common.method[0].ordering = ordering;

#if 0
  _common.supernodal = CHOLMOD_SIMPLICIAL;
#endif
  _common.supernodal = CHOLMOD_SUPERNODAL;
#endif
#if 0
  _common.nmethods = 1;
  _common.method[0].ordering = CHOLMOD_NESDIS;
#endif

  _common.final_ll = 1;
  //_common.supernodal = CHOLMOD_SIMPLICIAL;

  _common.error_handler = CHOLMOD_Environment_Long::handleError;
}

//////////////////////////////////////////////////////////////////////
// Destructor
//////////////////////////////////////////////////////////////////////
CHOLMOD_Environment_Long::~CHOLMOD_Environment_Long()
{
  freeSparse( _matrix );
  freeSparse( _factorSparse );

  freeTriplet( _matrixTriplet );
  freeTriplet( _factorTriplet );

  freeFactor( _factor );

  cholmod_l_finish( &_common );
}

//////////////////////////////////////////////////////////////////////
// Sets a new matrix to factorize
//////////////////////////////////////////////////////////////////////
void CHOLMOD_Environment_Long::setMatrix( const SPARSE_MATRIX &M )
{
  clear();

  _factor = NULL;

  cout << SDUMP( M(0,0) ) << endl;

  allocateSparse( M, _matrix, _matrixTriplet, true );
}

//////////////////////////////////////////////////////////////////////
// Sets a matrix from a stored sparse matrix file
//
// format = 0 (SPARSE_MATRIX object written to file)
// format = 1 (SPARSE_MATRIX::SparseColumnCopy object written to file)
//////////////////////////////////////////////////////////////////////
void CHOLMOD_Environment_Long::setMatrix( const char *filename, int format )
{
  clear();

  _factor = NULL;

  allocateSparse( filename, format, _matrix, _matrixTriplet, true );
}

//////////////////////////////////////////////////////////////////////
// Factors the current matrix
//////////////////////////////////////////////////////////////////////
void CHOLMOD_Environment_Long::computeCholeskyFactor( bool numeric,
                                                      bool symbolic,
                                                      long *userPerm )
{
  if ( !_matrix )
    return;

  if ( symbolic )
  {
    if ( userPerm )
    {
      _factor = cholmod_l_analyze_p( _matrix, userPerm, NULL, 0, &_common );
    }
    else
    {
      cout << "Analyzing" << endl;
      _factor = cholmod_l_analyze( _matrix, &_common );
    }
  }

  if ( numeric )
  {
    _factorTimer.tick();
    cholmod_l_factorize( _matrix, _factor, &_common );
    _factorTimer.tock();
  }

  cholmod_l_print_sparse( _matrix, "Cholesky Factor Input", &_common );
  cholmod_l_print_factor( _factor, "Cholesky Factor Output", &_common );

  if ( _common.method[ _common.selected ].ordering == CHOLMOD_AMD )
  {
    cout << "AMD used!!!" << endl;
  }
  else if ( _common.method[ _common.selected ].ordering == CHOLMOD_METIS )
  {
    cout << "METIS used!!!" << endl;
  }
  else
  {
    cout << "Some other ordering was used :(" << endl;
  }

  cout << SDUMP( _factor->nsuper ) << endl;
  cout << SDUMP( _factor->nzmax ) << endl;

  // FIXME
  long int totalNz = 0;
  long *super = (long *)_factor->super;
  long *pi = (long *)_factor->pi;

  for ( int super_idx = 0; super_idx < _factor->nsuper; super_idx++ )
  {
    long nCols = super[ super_idx + 1 ] - super[ super_idx ];
    long nRows = pi[ super_idx + 1 ] - pi[ super_idx ];

    totalNz += nCols * nRows;
    //totalNz -= nCols * ( nCols - 1 ) / 2;
  }
  cout << SDUMP( totalNz ) << endl;
  cout << SDUMP( (Real)totalNz * 8.0 / 1024.0 / 1024.0 ) << endl;
}

//////////////////////////////////////////////////////////////////////
// Puts the factor in the given sparse matrix
//////////////////////////////////////////////////////////////////////
void CHOLMOD_Environment_Long::getFactor( SPARSE_MATRIX &G )
{
  if ( !_factor )
    return;

  freeSparse( _factorSparse );
  freeTriplet( _factorTriplet );

#if 0
  G.resize( _factor->n, _factor->n );
  G.clear();
#endif

  cholmod_factor *factorCopy = cholmod_l_copy_factor( _factor, &_common );

  _factorSparse = cholmod_l_factor_to_sparse( factorCopy, &_common );
  _factorTriplet = cholmod_l_sparse_to_triplet( _factorSparse, &_common );

  cholmod_l_free_factor( &factorCopy, &_common );

#if 0
  int *row_indices = (int *)_factorTriplet->i;
  int *col_indices = (int *)_factorTriplet->j;
  double *matrix_entries = (double *)_factorTriplet->x;

  for ( int i = 0; i < _factorTriplet->nnz; i++ )
  {
    int row_idx = row_indices[ i ];
    int col_idx = col_indices[ i ];

    double value = matrix_entries[ i ];

    G( row_idx, col_idx ) = value;
  }
#endif

  tripletToSparse( _factorTriplet, G );

  // If we have an LDL' factorization, convert it to a
  // standard Cholesky factor.
  if ( !_factor->is_ll )
  {
    VECTOR diagonals( _factor->n );

    for ( long i = 0; i < _factor->n; i++ )
    {
      diagonals(i) = sqrt( G(i, i) );
      G(i, i) = diagonals(i);
    }

    cout << "Converting LDL' factorization to Cholesky" << endl;

    for ( SPARSE_MATRIX::iterator i = G.begin(); i != G.end(); i++ )
    {
      int row_idx = i->first.first;
      int col_idx = i->first.second;

      if ( row_idx == col_idx )
        continue;

      i->second *= diagonals( col_idx );
    }
  }
  else
  {
    cout << "Factorization already in Cholesky format: "
            "no conversion needed" << endl;
  }
}

//////////////////////////////////////////////////////////////////////
// Cleanup
//////////////////////////////////////////////////////////////////////
void CHOLMOD_Environment_Long::clear()
{
  freeSparse( _matrix );
  freeTriplet( _matrixTriplet );
  freeFactor( _factor );
  freeSparse( _factorSparse );
  freeTriplet( _factorTriplet );
}

//////////////////////////////////////////////////////////////////////
// Converts our sparse matrix format to CHOLMOD's
//////////////////////////////////////////////////////////////////////
void CHOLMOD_Environment_Long::allocateSparse( const SPARSE_MATRIX &M,
                                          cholmod_sparse * &A,
                                          cholmod_triplet * &T,
                                          bool symmetric,
                                          bool transpose )
{
  // This assumes that the diagonal has non-zeros on it
  long nnz;
  long idx = 0;

  if ( symmetric )
  {
    nnz = ( M.size() - M.cols() ) / 2 + M.cols();
  }
  else
  {
    nnz = M.size();
  }

  T = cholmod_l_allocate_triplet( transpose ? M.cols() : M.rows(),
                                  transpose ? M.rows() : M.cols(),
                                  nnz,
                                  symmetric ? -1 : 0 /* symmetric, lower tri */,
                                  CHOLMOD_REAL, &_common );

  long *row_indices = (long *)T->i;
  long *col_indices = (long *)T->j;
  double *matrix_entries = (double *)T->x;

  // Put the matrix information in to triplet form.  Pay attention
  // to only lower triangular entries.
  for ( SPARSE_MATRIX::const_iterator i = M.begin(); i != M.end(); i++ )
  {
    int row_idx = i->first.first;
    int col_idx = i->first.second;

    if ( row_idx < col_idx && symmetric )
      continue;

    row_indices[ idx ] = transpose ? col_idx : row_idx;
    col_indices[ idx ] = transpose ? row_idx : col_idx;
    matrix_entries[ idx ] = i->second;

    idx++;
  }

  T->nnz = nnz;

  A = cholmod_l_triplet_to_sparse( T, 0, &_common );
}

//////////////////////////////////////////////////////////////////////
// Reads in a sparse matrix from a file and puts it in CHOLMOD's
// internal format
//
// format = 0 (SPARSE_MATRIX object written to file)
// format = 1 (SPARSE_MATRIX::SparseColumnMatrix object written to file)
//////////////////////////////////////////////////////////////////////
void CHOLMOD_Environment_Long::allocateSparse( const char *filename,
                                          int format,
                                          cholmod_sparse * &A,
                                          cholmod_triplet * &T,
                                          bool symmetric,
                                          bool transpose )
{
  if ( format == 0 )
  {
  }
  else if ( format == 1 )
  {
    allocateSparseColumnCopy( filename, A, T, symmetric, transpose );
  }
  else
  {
    cerr << "** WARNING ** Invalid format" << endl;
  }
}

//////////////////////////////////////////////////////////////////////
// Different versions of the above for different formats
//
// Assumes matrix written by SPARSE_MATRIX object
//////////////////////////////////////////////////////////////////////
void CHOLMOD_Environment_Long::allocateSparseStandard(
                                                  const char *filename,
                                                  cholmod_sparse * &A,
                                                  cholmod_triplet * &T,
                                                  bool symmetric,
                                                  bool transpose )
{
}

//////////////////////////////////////////////////////////////////////
// Assumes matrix written by SparseColumnCopy object
//////////////////////////////////////////////////////////////////////
void CHOLMOD_Environment_Long::allocateSparseColumnCopy(
                                                    const char *filename,
                                                    cholmod_sparse * &A,
                                                    cholmod_triplet * &T,
                                                    bool symmetric,
                                                    bool transpose )
{
  SPARSE_MATRIX::SparseColumnMatrix S;

  SPARSE_MATRIX::readFromBinary( S, filename );

  // Allocate triplet version
  T = cholmod_l_allocate_triplet( transpose ? S._ncol : S._nrow,
                                transpose ? S._nrow : S._ncol,
                                S._nzmax,
                                symmetric ? 1 : 0 /* symmetric, lower tri */,
                                CHOLMOD_REAL, &_common );

  long *row_indices = (long *)T->i;
  long *col_indices = (long *)T->j;
  double *matrix_entries = (double *)T->x;
  long int idx = 0;

  size_t nzTotal = 0;

  for ( int col_idx = 0; col_idx < S._ncol; col_idx++ )
  {
    for ( int row_ptr = S._p[ col_idx ]; row_ptr < S._p[ col_idx + 1 ];
          row_ptr++ )
    {
      int row_idx = S._i[ row_ptr ];

      if ( symmetric && row_idx < col_idx )
      {
        continue;
      }

      TRACE_ASSERT( row_idx >= 0 && row_idx < S._nrow, "Row out of range" );

      row_indices[ idx ] = transpose ? col_idx : row_idx;
      col_indices[ idx ] = transpose ? row_idx : col_idx;
      matrix_entries[ idx ] = S._x[ row_ptr ];

      idx++;

      nzTotal++;
    }
  }

  T->nnz = S._nzmax;

  printf( "Found %zd non-zeros in matrix %s\n", nzTotal, filename );

#if 0
  // WAY too slow
  printf( "Computing sparse fill-in\n" );
  nzTotal = SPARSE_MATRIX::computeFillIn( S );
  printf( "Unordered factor has %lld non-zeros\n\n", nzTotal );
#endif

  A = cholmod_l_triplet_to_sparse( T, 0, &_common );
}

//////////////////////////////////////////////////////////////////////
// Cleans up a sparse matrix
//////////////////////////////////////////////////////////////////////
void CHOLMOD_Environment_Long::freeSparse( cholmod_sparse * &A )
{
  if ( A )
  {
    cholmod_l_free_sparse( &A, &_common );
  }
}

//////////////////////////////////////////////////////////////////////
// Cleans up a sparse matrix in triplet form
//////////////////////////////////////////////////////////////////////
void CHOLMOD_Environment_Long::freeTriplet( cholmod_triplet * &T )
{
  if ( T )
  {
    cholmod_l_free_triplet( &T, &_common );
  }
}

//////////////////////////////////////////////////////////////////////
// Cleans up a sparse matrix factor
//////////////////////////////////////////////////////////////////////
void CHOLMOD_Environment_Long::freeFactor( cholmod_factor * &G )
{
  if ( G )
  {
    cholmod_l_free_factor( &G, &_common );
  }
}

//////////////////////////////////////////////////////////////////////
// Error handler for CHOLMOD
//////////////////////////////////////////////////////////////////////
void CHOLMOD_Environment_Long::handleError( int status, const char *file,
                                       int line, const char *message )
{
  printf ("cholmod error: file: %s line: %d status: %d: %s\n",
          file, line, status, message) ;
  if (status < 0)
  {
    exit (0) ;
  }
}

//////////////////////////////////////////////////////////////////////
// Converts the CHOLMOD triplet format to our SPARSE_MATRIX
// format
//////////////////////////////////////////////////////////////////////
void CHOLMOD_Environment_Long::tripletToSparse( const cholmod_triplet *T,
                                           SPARSE_MATRIX &A )
{
  A.resize( T->nrow, T->ncol );
  A.clear();

  long *row_indices = (long *)T->i;
  long *col_indices = (long *)T->j;
  double *matrix_entries = (double *)T->x;

  for ( int i = 0; i < T->nnz; i++ )
  {
    int row_idx = row_indices[ i ];
    int col_idx = col_indices[ i ];

    double value = matrix_entries[ i ];

    A( row_idx, col_idx ) = value;
  }
}
