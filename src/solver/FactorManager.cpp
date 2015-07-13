//////////////////////////////////////////////////////////////////////
// FactorManager.cpp: Implementation of the FactorManager class
//
//////////////////////////////////////////////////////////////////////

#include "FactorManager.h"

#include <node/DiagonalBlock.h>

#if 0
#include <ordering/MinimumDegree.h>
#endif

#include <util/MathUtil.h>
#include <util/StatsCounter.h>
#include <util/STLUtil.h>

#include <boost/bind.hpp>
#include <boost/function.hpp>

#include <map>

const int FactorManager::EMPTY = -1;

const Real FactorManager::RANK_CONSTANT = 0.4;
//const Real FactorManager::RANK_CONSTANT = 0.5;
//const Real FactorManager::RANK_CONSTANT = 1.0;
//const Real FactorManager::RANK_CONSTANT = 0.0025;

const Real FactorManager::DIAGONAL_RANK_CONSTANT = 0.5;
//const Real FactorManager::DIAGONAL_RANK_CONSTANT = 0.75;
//const Real FactorManager::DIAGONAL_RANK_CONSTANT = 1.0;
//const Real FactorManager::DIAGONAL_RANK_CONSTANT = 1.0;

const int FactorManager::OVER_SAMPLING = 5;
//const int FactorManager::OVER_SAMPLING = 2;

const int FactorManager::POWER_ITERATIONS = 1;

// FIXME: for implementing interior blocks
#define FINALIZE_INTERIOR_BLOCKS 1
//#undef FINALIZE_INTERIOR_BLOCKS

const char *FactorManager::TIMER_NAMES[] = 
{
  "Base matrix copy",
  "Block initialization",
  "Build interaction intersection",
  "Build relative map",
  "Apply standard update",
  "Apply ext. update (decompress)",
  "Apply ext. update (apply)",
  "Build low rank descendent list",
  "Orthogonal projection",
  "QR factorization",
  "Basis copy",
  "Random multiply number generation",
  "Random multiply initialization",
  "Random multiply application",
  "Update active blocks",
  "Error estimate",
  "Basis projection",
  "Basis assignment",
  "Compressed interaction expansion",
  "Find diagonal block ranges",
  "Diagonal orthogonal projection",
  "Diagonal QR factoriation",
  "Diagonal basis copy",
  "Diagonal random multiply",
  "Diagonal update active blocks",
  "Diagonal error estimate",
  "Diagonal basis projection",
  "Diagonal basis assignment",
  "Add diagonal contributions",
  "Factor node",
  "Init extended nodes",
  "Build extended schur complement",
  "Build extended schur (inversion)",
  "Build extended schur (multiply)",
  "Interior block: preamble",
  "Interior block: node list setup",
  "Interior block: workspace setup",
  "Interior block: multiply setup",
  "Interior block: sparse multiply",
  "Interior block: build multiply workspace",
  "Interior block: forward solve",
  "Interior block: backward solve",
  "Interior block: solve - triangular solve",
  "Interior block: solve - propagate solve result",
  "Interior block: copy result",
};

//////////////////////////////////////////////////////////////////////
// Constructor: Assumes that an ordering for the matrix
// has been generated and is stored in the nodes array
//////////////////////////////////////////////////////////////////////
FactorManager::FactorManager(
        const SPARSE_MATRIX::SparseColumnMatrix &A,
        //const std::vector<const NestedDissection::DissectionNode *> &nodes,
        const std::vector<Ordering::SupernodeSpecification> &nodes,
        const std::vector<std::vector<DenseBlock> > *diagonalBlocks,
        Real errorBoundMultiplier,
        Real rankConstant, Real diagonalRankConstant )
  : _A( A ),
    _ordering( nodes ),
    _diagonalBlocks( diagonalBlocks ),
    _factorData( NULL ),
    _factorDataSz( 0 ),
    _multWorkspace( NULL ),
    _updateWorkspace( NULL ),
    _solveWorkspace( NULL ),
    _extendedSolveWorkspace( NULL ),
    _extendedSolveWorkspaceSz( 0 ),
    _multWorkspaceSz( 0 ),
    _updateWorkspaceSz( 0 ),
    _solveWorkspaceSz( 0 ),
    _overSampling( OVER_SAMPLING ),
    _slackData( NULL ),
    _slackDataSz( 0 ),
    _availableSlackData( 0 ),
    _slackDataOffset( 0 ),
#if 0
    _subMatrix( NULL ),
#endif
    _decompWorkspace( NULL ),
    _decompWorkspaceSz( 0 ),
#if 0
    _decompMultWorkspace( NULL ),
    _decompMultWorkspaceSz( 0 ),
#endif
    _decompMultTransWorkspace( NULL ),
    _decompMultTransWorkspaceSz( 0 ),
#if 0
    _copyWorkspace( NULL ),
    _copyWorkspaceSz( 0 ),
#endif
    _qrExtraData( NULL ),
    _qrExtraDataSz( 0 ),
    _qrWorkspace( NULL ),
    _qrWorkspaceSz( 0 ),
    _errorBoundMultiplier( errorBoundMultiplier ),
#if 0
    _blockWorkspace( NULL ),
    _blockWorkspaceSz( 0 ),
#endif
    _vectorWork1( NULL ),
    _vectorWork2( NULL ),
    _vectorWorkSz( 0 ),
    _vectorWork1Sz( 0 ),
    _vectorWork2Sz( 0 ),
    _basisWorkspace( NULL ),
    _basisWorkspaceSz( 0 ),
#if 0
    _expansionWorkspace( NULL ),
    _expansionWorkspaceSz( 0 ),
#endif
    // Timers
    _totalTime( "Total" ),
    _useInteriorBlocks( false ),
    _offDiagCompressionType( COMPRESS_INDIVIDUAL_INTERACTIONS ),
    _offDiagonalRankEstimator( boost::bind( &FactorManager::blockRank,
                                            this, _1, _2 ) ),
    _diagonalRankEstimator( boost::bind( &FactorManager::diagonalBlockRank,
                                         this, _1, _2 ) ),
    _rankConstant( rankConstant ),
    _diagonalRankConstant( diagonalRankConstant )
{
  TRACE_ASSERT( A._nrow == A._ncol, "System is not square" );

  // Make sure we have valid rank constants
  _rankConstant = max( _rankConstant, RANK_CONSTANT );
  _diagonalRankConstant = max( _diagonalRankConstant, RANK_CONSTANT );

  initTimers();
}

//////////////////////////////////////////////////////////////////////
// Destructor
//////////////////////////////////////////////////////////////////////
FactorManager::~FactorManager()
{
#if 0
  if ( _factorDataSz > 0 ) delete[] _factorData;
  if ( _multWorkspaceSz > 0 ) delete[] _multWorkspace;
  if ( _updateWorkspaceSz > 0 ) delete[] _updateWorkspace;
  delete[] _solveWorkspace;
  if ( _extendedSolveWorkspaceSz > 0 ) delete[] _extendedSolveWorkspace;
#if 0
  delete[] _subMatrix;
#endif
  if ( _decompWorkspaceSz > 0 ) delete[] _decompWorkspace;
#if 0
  delete[] _decompMultWorkspace;
#endif
  if ( _decompMultTransWorkspaceSz > 0 ) delete[] _decompMultTransWorkspace;
#if 0
  delete[] _copyWorkspace;
#endif
  if ( _qrExtraDataSz > 0 ) delete[] _qrExtraData;
  if ( _qrWorkspaceSz > 0 ) delete[] _qrWorkspace;
#if 0
  delete[] _blockWorkspace;
#endif
  if ( _vectorWorkSz > 0 ) delete[] _vectorWork1;
  if ( _vectorWorkSz > 0 ) delete[] _vectorWork2;
  if ( _basisWorkspaceSz > 0 ) delete[] _basisWorkspace;
#if 0
  delete[] _expansionWorkspace;
#endif

  if ( _slackDataSz > 0 ) delete[] _slackData;
#endif

  clear();
}

#define FACTORMANAGER_FREE_WORKSPACE( workspaceName ) \
  if ( workspaceName ## Sz > 0 ) { \
    delete[] workspaceName; \
    workspaceName = NULL; \
    workspaceName ## Sz = 0; \
  }

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void FactorManager::clear()
{
  FACTORMANAGER_FREE_WORKSPACE( _factorData );
  FACTORMANAGER_FREE_WORKSPACE( _multWorkspace );
  FACTORMANAGER_FREE_WORKSPACE( _updateWorkspace );
  FACTORMANAGER_FREE_WORKSPACE( _solveWorkspace );
  FACTORMANAGER_FREE_WORKSPACE( _extendedSolveWorkspace );
  FACTORMANAGER_FREE_WORKSPACE( _decompWorkspace );
  FACTORMANAGER_FREE_WORKSPACE( _decompMultTransWorkspace );
  FACTORMANAGER_FREE_WORKSPACE( _qrExtraData );
  FACTORMANAGER_FREE_WORKSPACE( _qrWorkspace );
  FACTORMANAGER_FREE_WORKSPACE( _vectorWork1 );
  FACTORMANAGER_FREE_WORKSPACE( _vectorWork2 );
  FACTORMANAGER_FREE_WORKSPACE( _basisWorkspace );
  FACTORMANAGER_FREE_WORKSPACE( _slackData );
}

//////////////////////////////////////////////////////////////////////
// Builds a symbolic factor for this manager's system,
// optionally extending the system with slack variables.
//////////////////////////////////////////////////////////////////////
void FactorManager::factorSymbolic( int maxBlockSize,
                                    CompressionType compressType,
                                    bool useInteriorBlocks )
{
  vector<Supernode>          initialNodes;

  bool                       compressInPlace;

  compressInPlace = ( compressType == IN_PLACE );

  _factor.clear();

  cout << SDUMP( _diagonalBlocks->size() ) << endl;
  cout << SDUMP( _ordering.size() ) << endl;
  TRACE_ASSERT(
      _diagonalBlocks == NULL || _diagonalBlocks->size() == _ordering.size(),
      "Invalid diagonal block set" );

  Supernode::BuildBasicFactor( _ordering, _A, initialNodes );

  Supernode::BuildOffDiagonalInteractions( initialNodes, _A, compressInPlace );

  // Prior to system extension, set diagonal block ranges for each node
  if ( _diagonalBlocks != NULL )
  {
    for ( int node_idx = 0; node_idx < initialNodes.size(); node_idx++ )
    {
      initialNodes[ node_idx ].setDiagonalBlocks(
                                        //_diagonalBlocks->at( node_idx ) );
                                        _diagonalBlocks->at( node_idx ),
                                        compressInPlace,
                                        //-1 );
                                        4000 );
                                        //2000 );
    }
  }

  // Set up a block rank estimator for symbolic reordering of
  // extended variables
  Supernode::RankEstimator   estimator = boost::bind( &FactorManager::blockRank,
                                                      this, _1, _2 );

  Supernode::ExtendSystemAppend( initialNodes, _factor, maxBlockSize, _A,
                                 2.0 /* compression ratio */,
                                 //NULL );
                                 &estimator,
                                 compressInPlace );

  // Count non-zero columns in each block of the factor
  for ( int node_idx = 0; node_idx < _factor.size(); node_idx++ )
  {
    _factor[ node_idx ].countNonZeroInteractionColumns( _factor, _A );
    _factor[ node_idx ].countNonZeroDiagonalBlockColumns( _A );
  }

  _numExtendedNodes = _factor.size() - initialNodes.size();

  if ( useInteriorBlocks ) {
    identifyInteriorBlocks();
    printInteriorBlockDiagonalUsage();
  }

  calculateSlackUsage( maxBlockSize );

  initWorkspace();

  initDecompositionWorkspace();

  writeStandardMemoryUsage( maxBlockSize );
}

//////////////////////////////////////////////////////////////////////
// Computes the numerical factor for our system, assuming that
// the symbolic factor has already been computed.
// The pattern of A is assumed to match the pattern of the
// stored matrix.  No error checking is done for this, however.
//////////////////////////////////////////////////////////////////////
SolverErrorPtr FactorManager::factorNumeric(
                                  const SPARSE_MATRIX::SparseColumnMatrix &A,
                                  bool adaptive, Real tolerance,
                                  Real diagonalTolerance,
                                  int writeSize, bool writeFactors,
                                  ExtendedFactorType extFactorType,
                                  CompressionType compressType,
                                  DecompositionBlock blockType,
                                  DecompositionType decompType )
{
  int                        descendent_idx;
  int                        descendent_idx_next;
  int                        start_idx;
  int                        nCols, nRows1, nRows2;
  int                        parentNode;
  int                        descendentAncestor;
  int                        rank;
  int                        interaction_idx;

  int                        node_idx = 0;

  bool                       decomposeTranspose;
  bool                       useFactor;

  SolverErrorPtr             error;

  // FIXME: reseed the random number generator
  MathUtil::GENERATOR.seed( 152 );

  printf( "Factoring %d x %d sparse matrix\n", (int)A._nrow, (int)A._ncol );

  decomposeTranspose = ( blockType == DECOMPOSE_TRANSPOSE );
  useFactor = ( decompType == DECOMPOSE_FACTOR );

  clearWorkspace();

  _totalTime.tick();

  // TODO: refactoring

  // Factor all standard nodes
  for ( ; node_idx < _factor.size()
            && _factor[ node_idx ].type() == STANDARD_NODE; node_idx++ )
  {
    Supernode                 &node = _factor[ node_idx ];

    TRACE_ASSERT( node.type() == STANDARD_NODE );

    //printf( "Node %09d of %09d\r", node_idx + 1, (int)_factor.size() );

#ifdef DO_TIMING
    _timers[ BASE_MATRIX_COPY ].tick();
#endif

    // Copy data from the columns of A handled by this node
    node.constructScatteredMap( _scatteredMap, _lowRankScatteredMap, _factor );

    node.copyMatrixData( _scatteredMap, _lowRankScatteredMap, A, _slackData );

#ifdef DO_TIMING
    _timers[ BASE_MATRIX_COPY ].tock();
#endif

    //////////////////////////////////////////////////////////////////////
    // Prepare workspaces for low rank decomposition
    //////////////////////////////////////////////////////////////////////

#ifdef DO_TIMING
    _timers[ BLOCK_INIT ].tick();
#endif

    // Get the rank for off-diagonal components of this node.
    // FIXME: This is a hack.  For now, we are assigning all
    //        off-diagonal components the same rank (based on the
    //        number of columns.
    rank = maxBlockRank( 0, node.numColumns() );
    rank = min( rank, node.numColumns() );

    // Figure out which blocks are to be treated as low rank
    node.getLowRankBlocks( _lowRankBlocks, _lowRankBaseRows );

    setBlockRanks( node );

    // Set all blocks to active for now (since we are not done
    // factoring any of them yet)
    _lowRankBlockActive.clear();
    _lowRankBlockActive.resize( _lowRankBlocks.size(), true );

    _currentDescendents.clear();

    if ( _lowRankBlocks.size() > 0 )
    {
      // FIXME: better way of doing this
      for ( int i = 0; i < _lowRankDescendents.size(); i++ )
      {
        _lowRankDescendents[ i ].clear();
      }
    }

#ifdef DO_TIMING
    _timers[ BLOCK_INIT ].tock();
#endif

    //////////////////////////////////////////////////////////////////////
    // Iterate over each descendent of this node
    //////////////////////////////////////////////////////////////////////
    for ( descendent_idx = _head[ node_idx ]; descendent_idx != EMPTY;
          descendent_idx = descendent_idx_next )
    {
      Supernode               &descendent = _factor[ descendent_idx ];

      if ( descendent.compressOffDiagonal() ) {
        TRACE_ASSERT( node.compressOffDiagonal() );
      }

      // Keep track of all descendents for this node
      _currentDescendents.push_back( descendent_idx );

      descendent_idx_next = _next[ descendent_idx ];

      //////////////////////////////////////////////////////////////////////
      // Form Schur complement for this node
      //////////////////////////////////////////////////////////////////////
      start_idx = _nextInteraction[ descendent_idx ];

#ifdef DO_TIMING
      _timers[ INTERACTION_INTERSECTION ].tick();
#endif

      node.findInteractionIntersection( descendent, start_idx,
                                        _interactionIndices,
                                        _lowRankInteractionIndices,
                                        _extendedInteractionIndices );

#ifdef DO_TIMING
      _timers[ INTERACTION_INTERSECTION ].tock();
#endif

#ifdef DO_TIMING
      _timers[ RELATIVE_MAP ].tick();
#endif

      if ( !descendent.compressOffDiagonal() ) {
        node.buildRelativeMap( descendent, start_idx, _interactionIndices,
                               _relativeMap );
      }

#ifdef DO_TIMING
      _timers[ RELATIVE_MAP ].tock();
#endif

#ifdef DO_TIMING
      _timers[ STANDARD_UPDATE ].tick();
#endif
      FLOP_COUNT_START;
      TIMING_START( "Standard update" );

      // Copy the relevant interaction data from descendent into
      // another work space, and multiply to form an update matrix.
      // Skip this step for extended nodes - it is not needed
      if ( node.type() == STANDARD_NODE ) {
        if ( !descendent.compressOffDiagonal() ) {
          descendent.copyInteractionData( descendent, start_idx,
                                          _interactionIndices,
                                          _multWorkspace,
                                          nCols, nRows1, nRows2 );

          node.buildUpdateMatrix( _multWorkspace,
                                  nCols, nRows1, nRows2,
                                  _relativeMap,
                                  _updateWorkspace,
                                  descendent_idx );

          // Subtract the update matrix from this node's data
          // (requires the relative map to figure out correct indices)
          node.applyUpdate( _relativeMap, nRows1, nRows2, _updateWorkspace );
        }

        // Apply an updates resulting from compressed interactions
        // in the descendent.  Only necessary if node doesn't
        // build large diagonal explicitly.
        if ( !node.nodeDiagonal().hasExplicitBlocks() ) {
          node.compressedDiagonalUpdate( descendent, start_idx );
        }
      }

      TIMING_STOP( "Standard update" );
      FLOP_COUNT_END( "Standard update" );
#ifdef DO_TIMING
      _timers[ STANDARD_UPDATE ].tock();
#endif

      //_timers[ EXTENDED_UPDATE ].tick();

      // Form update matrices for interactions with extended nodes
      FLOP_COUNT_START;
      TIMING_START( "Extended update (standard)" );
      node.applyExtendedUpdate( descendent, _extendedInteractionIndices,
                                _slackData, start_idx,
                                _realWorkspaceManager,
#ifdef DO_TIMING
                                &_timers[ EXTENDED_UPDATE_STANDARD ],
                                &_timers[ EXTENDED_UPDATE_EXTENDED ]
#else
                                NULL, NULL
#endif
                                );
      TIMING_STOP( "Extended update (standard)" );
      FLOP_COUNT_END( "Extended update (standard)" );

      //_timers[ EXTENDED_UPDATE ].tock();

      //////////////////////////////////////////////////////////////////////
      // First pass for building low rank decompositions
      //////////////////////////////////////////////////////////////////////

#ifdef DO_TIMING
      _timers[ LOW_RANK_DESC_LIST ].tick();
#endif

      if ( _lowRankBlocks.size() > 0 ) {
        // For each low rank block, add to its list of descendents, if
        // this descendent interacts with that block
        node.appendLowRankDescendentList( descendent, _lowRankInteractionIndices,
                                          _lowRankBlocks, _lowRankDescendents );
      }

#ifdef DO_TIMING
      _timers[ LOW_RANK_DESC_LIST ].tock();
#endif

      //////////////////////////////////////////////////////////////////////
      // Update descendent relationships
      //////////////////////////////////////////////////////////////////////
      _nextInteractionCache[ descendent_idx ] = start_idx;

      TRACE_ASSERT(
          descendent.offDiagonal()[ start_idx ]._nodeID == node.nodeID(),
          "Invalid ancestor relationship" );

      start_idx++;
      _nextInteraction[ descendent_idx ] = start_idx;

      // Increment the start index until we find an active node.
      // ie. skip any inactive nodes, since they no longer matter
      // for the factorization.
      while ( start_idx < descendent.offDiagonal().size()
           && !descendent.offDiagonal()[ start_idx ]._active )
      {
        start_idx++;
        _nextInteraction[ descendent_idx ] = start_idx;
      }

      if ( start_idx < descendent.offDiagonal().size() ) {
        descendentAncestor = descendent.offDiagonal()[ start_idx ]._nodeID;

        TRACE_ASSERT( descendentAncestor > node_idx
                   && descendentAncestor < _factor.size(),
                   "Invalid descendent ancestor" );

        _next[ descendent_idx ] = _head[ descendentAncestor ];
        _head[ descendentAncestor ] = descendent_idx;
      }
    }

    // Next, we need to figure out which interior blocks are descendents
    // of this node (but only if this node itself is not in an interior
    // block)
    if ( _useInteriorBlocks && !isInInteriorBlock( node_idx ) ) {
      findInteriorBlockDescendents();

#ifdef FINALIZE_INTERIOR_BLOCKS
      // Apply diagonal contributions from interior blocks.  Only necessary
      // if this node doesn't form large diagonal blocks explicitly.
      if ( !node.nodeDiagonal().hasExplicitBlocks() ) {
        applyInteriorBlockDiagonalUpdate( A, node );
      }
#endif
    }

    // If we are doing in-place compression of diagonal blocks, then
    // we should do the initial diagonal factorization here, since
    // we may want to use the diagonal in the process of low-rank
    // decomposition.
    //
    // FIXME: This doesn't currently account for diagonal compression
    if ( compressType == IN_PLACE ) {
      error = factorCompressedDiagonal( A, node, POWER_ITERATIONS, tolerance,
                                        adaptive, decomposeTranspose,
                                        useFactor );

      if ( error != SolverErrorPtr() ) {
        // Something has gone wrong, so return immediately
        return error;
      }
#if 0
      node.factor( _slackData, writeSize );
#endif
    }

    //////////////////////////////////////////////////////////////////////
    // First pass for low rank decomposition
    //////////////////////////////////////////////////////////////////////
    FLOP_COUNT_START;
    TIMING_START( "Off-diagonal decomposition" );
    if ( _lowRankBlocks.size() > 0 )
    {
      // Decompose based on the given tolerance
      offDiagonalDecomposition( A, node, POWER_ITERATIONS, tolerance,
                                adaptive, decomposeTranspose,
                                useFactor,
                                // FIXME
                                COMPRESS_FULL_OFF_DIAGONAL );
                                //COMPRESS_INDIVIDUAL_INTERACTIONS );
    }
    TIMING_STOP( "Off-diagonal decomposition" );
    FLOP_COUNT_END( "Off-diagonal decomposition" );

    //////////////////////////////////////////////////////////////////////
    // Second pass for low rank decomposition
    //////////////////////////////////////////////////////////////////////

    FLOP_COUNT_START;
    TIMING_START( "Off-diagonal basis projection" );
    if ( _lowRankBlocks.size() > 0 )
    {
      // Form projections in to each basis, assuming a fixed rank
      // for each block.
      //
      // Copies new extended data to the node.
      //fullBasisProject( A, node, writeFactors, adaptive ? 0 : rank );
      fullBasisProject( A, node, writeFactors,
                        0, true,
                        decomposeTranspose,
                        useFactor,
                        // FIXME
                        COMPRESS_FULL_OFF_DIAGONAL );
                        //COMPRESS_INDIVIDUAL_INTERACTIONS );
    }
    TIMING_STOP( "Off-diagonal basis projection" );
    FLOP_COUNT_END( "Off-diagonal basis projection" );

    //////////////////////////////////////////////////////////////////////
    // Now that we have expanded diagonal contributions, we can handle
    // low rank blocks in the main diagonal
    //////////////////////////////////////////////////////////////////////
    FLOP_COUNT_START;
    TIMING_START( "Diagonal decomposition" );
    if ( compressType == EXTENDED_VARIABLE ) {
      diagonalDecomposition( A, node, POWER_ITERATIONS, diagonalTolerance,
                             adaptive, false, false /* FIXME */ );
    }
    TIMING_STOP( "Diagonal decomposition" );
    FLOP_COUNT_END( "Diagonal decomposition" );

    //////////////////////////////////////////////////////////////////////
    // Add diagonal contributions resulting from slack variables to the
    // diagonal
    //
    //////////////////////////////////////////////////////////////////////
#ifdef DO_TIMING
    _timers[ ADD_DIAGONAL_CONTRIBUTIONS ].tick();
#endif
    FLOP_COUNT_START;
    TIMING_START( "Add low rank diag. contrib." );
    if ( compressType == EXTENDED_VARIABLE ) {
      node.addLowRankDiagonalContributions( _slackData, _realWorkspaceManager );
    }
    TIMING_STOP( "Add low rank diag. contrib." );
    FLOP_COUNT_END( "Add low rank diag. contrib." );
#ifdef DO_TIMING
    _timers[ ADD_DIAGONAL_CONTRIBUTIONS ].tock();
#endif

    //////////////////////////////////////////////////////////////////////
    // Factor, then triangular solve
    //////////////////////////////////////////////////////////////////////
#ifdef DO_TIMING
    _timers[ FACTOR_NODE ].tick();
#endif
    FLOP_COUNT_START;
    TIMING_START( "Factor node (standard)" );
    // If we are doing an extended variable factorization, we have
    // to compute the final factor here, since it may include diagonal
    // contributions from prior low-rank decompositions
    if ( compressType == EXTENDED_VARIABLE ) {
      node.factor( _slackData, writeSize );
    }
    node.offDiagonalSolve( // Whether or not to solve compressed blocks
                           // We need to do this if we only decomposed
                           // Schur complements
                           decompType == DECOMPOSE_SCHUR,
                           _slackData, writeSize );
    TIMING_STOP( "Factor node (standard)" );
    FLOP_COUNT_END( "Factor node (standard)" );
#ifdef DO_TIMING
    _timers[ FACTOR_NODE ].tock();
#endif

    // Find the first interaction which is not low rank
    start_idx = 0;

    while ( start_idx < node.offDiagonal().size()
         && !node.offDiagonal()[ start_idx ]._active )
    {
      start_idx++;
    }

    // Update descendent relationships
    if ( start_idx < node.offDiagonal().size() )
    {
      _nextInteraction[ node_idx ] = start_idx;

      parentNode = node.offDiagonal()[ start_idx ]._nodeID;

      _next[ node_idx ] = _head[ parentNode ];
      _head[ parentNode ] = node_idx;
    }

    // Link list for supernode s no longer needed (?)
    _head[ node_idx ] = EMPTY;

    node.clearScatteredMap( _scatteredMap, _lowRankScatteredMap, _factor );

#if 0
    // FIXME: debugging
    if ( node.numColumns() >= 200 ) {
      for ( int interaction_idx = 0; interaction_idx < node.offDiagonal().size();
            interaction_idx++ )
      {
        const SupernodeInteraction &interaction
                                      = node.offDiagonal()[ interaction_idx ];

        cout << "Writing a matrix?" << endl;

        const Real *data = node.interactionMatrix( interaction_idx );
        char buf[ 1024 ];
        sprintf( buf, "node_%d_interaction_%d.matrix", node.nodeID(),
                 interaction_idx );
        MATRIX::write( data, interaction._rowList.size(), node.numColumns(),
                       buf );

        if ( node.numColumns() >= 3000 )
          abort();
      }
    }
#endif
  }

  // Initialize slack variable memory for extended nodes
#ifdef DO_TIMING
  _timers[ INIT_EXTENDED ].tick();
#endif
  FLOP_COUNT_START;
  TIMING_START( "Init extended nodes" );
  initExtendedNodes( node_idx /* First extended node index */ );
  TIMING_STOP( "Init extended nodes" );
  FLOP_COUNT_END( "Init extended nodes" );
#ifdef DO_TIMING
  _timers[ INIT_EXTENDED ].tock();
#endif

  // Build the extended variable schur complement
  //_timers[ EXTENDED_SCHUR_FORMATION ].tick();
  FLOP_COUNT_START;
  TIMING_START( "Build extended Schur complement" );
  buildExtendedSchurComplement( node_idx /* Number of standard nodes */ );
  TIMING_STOP( "Build extended Schur complement" );
  FLOP_COUNT_END( "Build extended Schur complement" );
  //_timers[ EXTENDED_SCHUR_FORMATION ].tock();

#if 0
  //
  // Write out this block
  writeExtendedBlock( "extended_schur_complement.bcsm" );
#endif

  // Factor all extended nodes
  Timer extendedTimer( "extendedTimer" );
  extendedTimer.tick();
  for ( ; node_idx < _factor.size(); node_idx++ )
  {
    Supernode                 &node = _factor[ node_idx ];

    TRACE_ASSERT( node.type() == EXTENDED_NODE );

    // Set whether or not to use an LDL factorization
    if ( extFactorType == LDL ) {
      node.setLDL( true );
    } else {
      node.setLDL( false );
    }

    //printf( "Node %09d of %09d\r", node_idx + 1, (int)_factor.size() );

#ifdef DO_TIMING
    _timers[ BASE_MATRIX_COPY ].tick();
#endif

    // Copy data from the columns of A handled by this node
    node.constructScatteredMap( _scatteredMap, _lowRankScatteredMap, _factor );
    node.copyMatrixData( _scatteredMap, _lowRankScatteredMap, A, _slackData );

#ifdef DO_TIMING
    _timers[ BASE_MATRIX_COPY ].tock();
#endif

    //////////////////////////////////////////////////////////////////////
    // Prepare workspaces for low rank decomposition
    //////////////////////////////////////////////////////////////////////

#ifdef DO_TIMING
    _timers[ BLOCK_INIT ].tick();
#endif

    _currentDescendents.clear();

#ifdef DO_TIMING
    _timers[ BLOCK_INIT ].tock();
#endif

    //////////////////////////////////////////////////////////////////////
    // Iterate over each descendent of this node
    //////////////////////////////////////////////////////////////////////
    for ( descendent_idx = _head[ node_idx ]; descendent_idx != EMPTY;
          descendent_idx = descendent_idx_next )
    {
      Supernode               &descendent = _factor[ descendent_idx ];

      // Keep track of all descendents for this node
      _currentDescendents.push_back( descendent_idx );

      descendent_idx_next = _next[ descendent_idx ];

      //////////////////////////////////////////////////////////////////////
      // Form Schur complement for this node
      //////////////////////////////////////////////////////////////////////
      start_idx = _nextInteraction[ descendent_idx ];

#ifdef DO_TIMING
      _timers[ INTERACTION_INTERSECTION ].tick();
#endif

      node.findInteractionIntersection( descendent, start_idx,
                                        _interactionIndices,
                                        _lowRankInteractionIndices,
                                        _extendedInteractionIndices );

#ifdef DO_TIMING
      _timers[ INTERACTION_INTERSECTION ].tock();
#endif

#ifdef DO_TIMING
      _timers[ RELATIVE_MAP ].tick();
#endif

      node.buildRelativeMap( descendent, start_idx, _interactionIndices,
                             _relativeMap );

#ifdef DO_TIMING
      _timers[ RELATIVE_MAP ].tock();
#endif

#ifdef DO_TIMING
      _timers[ STANDARD_UPDATE ].tick();
#endif

#ifdef DO_TIMING
      _timers[ STANDARD_UPDATE ].tock();
#endif

      //_timers[ EXTENDED_UPDATE ].tick();

      // Form update matrices for interactions with extended nodes
#if 0
      node.applyExtendedUpdate( descendent, _extendedInteractionIndices,
                                _copyWorkspace, _slackData, start_idx,
                                _expansionWorkspace,
                                &_timers[ EXTENDED_UPDATE_STANDARD ],
                                &_timers[ EXTENDED_UPDATE_EXTENDED ] );
#endif
      FLOP_COUNT_START;
      TIMING_START( "Extended update (extended)" );
      node.applyExtendedUpdate( descendent, _extendedInteractionIndices,
                                _slackData, start_idx,
                                _realWorkspaceManager,
#ifdef DO_TIMING
                                &_timers[ EXTENDED_UPDATE_STANDARD ],
                                &_timers[ EXTENDED_UPDATE_EXTENDED ]
#else
                                NULL, NULL
#endif
                                );
      TIMING_STOP( "Extended update (extended)" );
      FLOP_COUNT_END( "Extended update (extended)" );

      //_timers[ EXTENDED_UPDATE ].tock();

      //////////////////////////////////////////////////////////////////////
      // Update descendent relationships
      //////////////////////////////////////////////////////////////////////
      _nextInteractionCache[ descendent_idx ] = start_idx;

      TRACE_ASSERT(
          descendent.offDiagonal()[ start_idx ]._nodeID == node.nodeID(),
          "Invalid ancestor relationship" );

      start_idx++;
      _nextInteraction[ descendent_idx ] = start_idx;

      // Increment the start index until we find an active node.
      // ie. skip any inactive nodes, since they no longer matter
      // for the factorization.
      while ( start_idx < descendent.offDiagonal().size()
           && !descendent.offDiagonal()[ start_idx ]._active )
      {
        start_idx++;
        _nextInteraction[ descendent_idx ] = start_idx;
      }

      if ( start_idx < descendent.offDiagonal().size() )
      {
        descendentAncestor = descendent.offDiagonal()[ start_idx ]._nodeID;

        TRACE_ASSERT( descendentAncestor > node_idx
                   && descendentAncestor < _factor.size(),
                   "Invalid descendent ancestor" );

        _next[ descendent_idx ] = _head[ descendentAncestor ];
        _head[ descendentAncestor ] = descendent_idx;
      }
    }

    //////////////////////////////////////////////////////////////////////
    // Factor, then triangular solve
    //////////////////////////////////////////////////////////////////////
#ifdef DO_TIMING
    _timers[ FACTOR_NODE ].tick();
#endif
    FLOP_COUNT_START;
    TIMING_START( "Factor node (extended)" );
    node.factor( _slackData, writeSize, &_realWorkspaceManager );
    node.offDiagonalSolve( // Whether or not to solve compressed blocks
                           // We need to do this if we only decomposed
                           // Schur complements
                           decompType == DECOMPOSE_SCHUR,
                           _slackData, writeSize, &_realWorkspaceManager );
    TIMING_STOP( "Factor node (extended)" );
    FLOP_COUNT_END( "Factor node (extended)" );
#ifdef DO_TIMING
    _timers[ FACTOR_NODE ].tock();
#endif
    //////////////////////////////////////////////////////////////////////

    // We must "un-permute" descendent interactions here, if we are
    // doing an LDL factorization.
    if ( extFactorType == LDL ) {
      for ( int i = 0; i < _currentDescendents.size(); i++ )
      {
        node.permuteDescendentInteraction(
                        _factor.at( _currentDescendents[ i ] ),
                        _nextInteractionCache.at( _currentDescendents[ i ] ),
                        _slackData );
      }
    }

    // Find the first interaction which is not low rank
    start_idx = 0;

    while ( start_idx < node.offDiagonal().size()
         && !node.offDiagonal()[ start_idx ]._active )
    {
      start_idx++;
    }

    // Update descendent relationships
    if ( start_idx < node.offDiagonal().size() )
    {
      _nextInteraction[ node_idx ] = start_idx;

      parentNode = node.offDiagonal()[ start_idx ]._nodeID;

      _next[ node_idx ] = _head[ parentNode ];
      _head[ parentNode ] = node_idx;
    }

    // Link list for supernode s no longer needed (?)
    _head[ node_idx ] = EMPTY;

    node.clearScatteredMap( _scatteredMap, _lowRankScatteredMap, _factor );
  }
  extendedTimer.tock();
  printf( "\n\nFactoring extended nodes took %f seconds\n\n",
          extendedTimer.getTotalSecs() );

  _totalTime.tock();

  printf("\n\nFactorization required %f MB of extra space\n",
         (Real)slackUsage() * 8.0 / 1024.0 / 1024.0 );

  printf( "\n" );

  return error;
}

//////////////////////////////////////////////////////////////////////
// Copies the off-diagonal slack variable data to a sparse matrix
//////////////////////////////////////////////////////////////////////
void FactorManager::copySlackOffDiagonal( SPARSE_MATRIX &S )
{
  IndexRange                 rowRange;
  IndexRange                 colRange;

  S.clear();

  // Resize S according to the number of extended variables
  S.resize( factorSystemSize() - _A._nrow, _A._nrow );

  rowRange.first = _A._nrow;
  rowRange.second = factorSystemSize() - 1;

  colRange.first = 0;
  colRange.second = _A._ncol - 1;

  for ( int node_idx = 0; node_idx < _factor.size(); node_idx++ )
  {
    _factor[ node_idx ].copySubMatrix( rowRange, colRange, _slackData,
                                       _factor, S );
  }
}

//////////////////////////////////////////////////////////////////////
// Check for any strange stuff happening in the factor
//////////////////////////////////////////////////////////////////////
bool FactorManager::checkFactor( int maxNode )
{
  int                        lastColumn;
  bool                       okay = true;

  lastColumn = _factor[ _factor.size() - _numExtendedNodes - 1 ].endColumn();

  if ( maxNode < 0 )
  {
    maxNode = _factor.size() - 1;
  }

  for ( int node_idx = 0; node_idx <= maxNode; node_idx++ )
  {
    const Supernode         &node = _factor[ node_idx ];

    for ( int interaction_idx = 0; interaction_idx < node.offDiagonal().size();
          interaction_idx++ )
    {
      const SupernodeInteraction
                            &interaction = node.offDiagonal()[ interaction_idx ];

      if ( interaction._active && interaction._type == STANDARD_NODE )
      {
        const IntArray      &rowList = interaction._rowList;

        for ( int row_idx = 0; row_idx < rowList.size(); row_idx++ )
        {
          int                real_row_idx = rowList[ row_idx ];

          TRACE_ASSERT( real_row_idx >= 0 && real_row_idx <= lastColumn,
                        "Bad row entry" );

          okay = false;
        }
      }
    }
  }

  return okay;
}

//////////////////////////////////////////////////////////////////////
// Counts the columns in the given list of nodes
//////////////////////////////////////////////////////////////////////
int FactorManager::nodeColumns( const IntArray &nodeList )
{
  int                        numColumns = 0;

  for ( int i = 0; i < nodeList.size(); i++ )
  {
    numColumns += _factor[ nodeList[ i ] ].numColumns();
  }

  return numColumns;
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void FactorManager::writeTimings()
{
  Real                       totalTime = _totalTime.getTotalSecs();
  Real                       otherTime = totalTime;
  FILE                      *timingFile;
  
  timingFile = fopen( "super_numeric/timings.dat", "w" );

  if ( !timingFile )
  {
    cerr << "Could not write timing info" << endl;
    return;
  }

  for ( int timer_idx = 0; timer_idx < NUM_TIMERS; timer_idx++ )
  {
    Timer                   &timer = _timers[ timer_idx ];
    Real                     timing = timer.getTotalSecs();
    Real                     percent = 100.0 * timing / totalTime;

    otherTime -= timing;

    fprintf( timingFile,
             "%38s:   %04.2f seconds,   %02.2f %%\n", timer.name().c_str(),
             timing, percent );
  }

  fprintf( timingFile, "%38s:   %04.2f seconds\n", "Other", otherTime );
  fprintf( timingFile, "%38s:   %04.2f seconds\n", "Total", totalTime );

  fclose( timingFile );
}

//////////////////////////////////////////////////////////////////////
// Once all standard nodes have been handled, call this to
// build a symbolic representation of the Schur complement formed
// by eliminating just the standard nodes in the system.
//
// Write this sparse matrix to the given file.
//////////////////////////////////////////////////////////////////////
void FactorManager::writeExtendedSymbolicSchurComplement( const char *filename )
{
  // Build an interaction list for each extended node
  vector<set<int> >          extendedInteractions;
  int                        numStandardNodes;
  int                        startColumn, endColumn;
  
  numStandardNodes = _factor.size() - _numExtendedNodes;
  startColumn = _factor[ numStandardNodes ].startColumn();
  endColumn = _factor.back().endColumn();

  TRACE_ASSERT( startColumn >= 0, "Not all standard nodes have been factored" );
  TRACE_ASSERT( endColumn >= 0, "Not all standard nodes have been factored" );

  buildExtendedSchurComplement( extendedInteractions );

  // Compute row lists for each column in the schur complement
  SPARSE_MATRIX::SparseColumnMatrix S;

  S._nrow = endColumn - startColumn + 1;
  S._ncol = endColumn - startColumn + 1;
  S._nzmax = 0;

  // Count non-zeros
  for ( int node_idx = 0; node_idx < _numExtendedNodes; node_idx++ )
  {
    int        nodeID = node_idx + numStandardNodes;
    int        nCols = _factor[ nodeID ].numColumns();

    TRACE_ASSERT( nodeID < _factor.size(), "Bad node ID" );

    set<int>  &nodeInteractions = extendedInteractions[ node_idx ];

    S._nzmax += nCols * nCols;

    for ( set<int>::const_iterator iter = nodeInteractions.begin();
          iter != nodeInteractions.end(); iter++ )
    {
      int      nextNodeID = *iter + numStandardNodes;

      TRACE_ASSERT( nextNodeID < _factor.size(), "Bad node ID" );

      S._nzmax += _factor[ nextNodeID ].numColumns() * nCols;
    }
  }

  // Initialize
  //S._p = new int[ endColumn - startColumn + 2 ];
  S._p = (int *)malloc( ( endColumn - startColumn + 2 ) * sizeof( int ) );
  //S._i = new int[ S._nzmax ];
  S._i = (int *)malloc( S._nzmax * sizeof( int ) );
  //S._x = new Real[ S._nzmax ];
  S._x = (Real *)malloc( S._nzmax * sizeof( Real ) );

  int          currentIndex = 0;

  memset( (void *)S._x, 0, S._nzmax * sizeof( Real ) );

  // Figure out actual data locations in matrix
  for ( int node_idx = 0; node_idx < _numExtendedNodes; node_idx++ )
  {
    int        nodeID = node_idx + numStandardNodes;
    int        nCols = _factor[ nodeID ].numColumns();
    int        columnEntries = 0;
    int        numColumns;
    int        start_idx = currentIndex;

    TRACE_ASSERT( nodeID < _factor.size(), "Bad node ID" );

    IndexRange columnRange( _factor[ nodeID ].startColumn() - startColumn,
                            _factor[ nodeID ].endColumn() - startColumn );

    TRACE_ASSERT( columnRange.first >= 0 && columnRange.second < S._nrow,
                  "Invalid range" );

    numColumns = range_size( columnRange );

    S._p[ columnRange.first ] = currentIndex;

    S._x[ currentIndex ] = 1.0;

    // Figure out first column for this block.  All other
    // columns will follow.
    for ( int i = columnRange.first; i <= columnRange.second; i++ )
    {
      //printf( "Column %d gets row %d\n", columnRange.first, i );
      S._i[ currentIndex ] = i;
      currentIndex++;
      columnEntries++;
    }

    set<int>  &nodeInteractions = extendedInteractions[ node_idx ];

    for ( set<int>::const_iterator iter = nodeInteractions.begin();
          iter != nodeInteractions.end(); iter++ )
    {
      int      nextNodeID = *iter + numStandardNodes;

      TRACE_ASSERT( nextNodeID < _factor.size(), "Bad node ID" );

      IndexRange nextColumnRange(
                      _factor[ nextNodeID ].startColumn() - startColumn,
                      _factor[ nextNodeID ].endColumn() - startColumn );

      TRACE_ASSERT(
            nextColumnRange.first >= 0 && nextColumnRange.second < S._nrow,
            "Invalid range" );

      for ( int i = nextColumnRange.first; i <= nextColumnRange.second; i++ )
      {
        S._i[ currentIndex ] = i;
        //printf( "Column %d gets row %d\n", columnRange.first, i );
        currentIndex++;
        columnEntries++;
      }
    }

    // Fill in the remaining columns
    for ( int col_idx = 1; col_idx < numColumns; col_idx++ )
    {
      columnEntries -= 1;
      start_idx += 1;

      // Copy column entries starting from the start index
      // to the next column
      S._p[ columnRange.first + col_idx ] = currentIndex;

      TRACE_ASSERT( start_idx + columnEntries <= currentIndex, "Bad memcpy" );

      memcpy( (void *)( S._i + currentIndex ), (void *)( S._i + start_idx ),
              columnEntries * sizeof( int ) );

      S._x[ currentIndex ] = 1.0;

      currentIndex += columnEntries;

      TRACE_ASSERT( currentIndex <= S._nzmax, "Ran out of space" );
    }
  }

  S._p[ endColumn - startColumn + 1 ] = currentIndex;

  printf( "\nWriting extended schur complement: %ld non-zeros\n",
          (long int)currentIndex );

  SPARSE_MATRIX::writeToBinary( S, filename );

  // Now continue the fill-in process inside the Schur complement
  printf( "Computing Schur complement fill-in\n" );
  long int schurFillIn = 0;

  for ( int node_idx = 0; node_idx < _numExtendedNodes; node_idx++ )
  {
    set<int> &columnInteractions = extendedInteractions.at( node_idx );

    //printf( "Node %d of %d\r", node_idx + 1, _numExtendedNodes );

    int nodeID = node_idx + numStandardNodes;
    TRACE_ASSERT( nodeID < _factor.size(), "Invalid node ID" );

    int nodeColumns = _factor[ nodeID ].numColumns();

    //schurFillIn += nodeColumns * ( nodeColumns + 1 ) / 2;
    schurFillIn += nodeColumns * nodeColumns;

    for ( set<int>::const_iterator i = columnInteractions.begin();
          i != columnInteractions.end(); i++ )
    {
      set<int> &ancestorInteractions = extendedInteractions.at( *i );

      int ancestorID = *i + numStandardNodes;
      TRACE_ASSERT( ancestorID < _factor.size(), "Invalid node ID ");

      int ancestorColumns = _factor[ ancestorID ].numColumns();

      schurFillIn += ancestorColumns * nodeColumns;

      set<int>::const_iterator j( i );
      j++;

      for ( ; j != columnInteractions.end(); j++ )
      {
        ancestorInteractions.insert( *j );
      }
    }
  }
  printf( "\n" );
  printf( "Filled schur complement has %lld non-zeros\n\n",
          (long long int)schurFillIn );
}

//////////////////////////////////////////////////////////////////////
// Writes the current contents of the extended variable block
// to a sparse matrix
//////////////////////////////////////////////////////////////////////
void FactorManager::writeExtendedBlock( const char *filename )
{
  if ( _numExtendedNodes == 0 )
    return;

  int              nCols;
  int              nodeColumns;
  const Supernode &firstNode = _factor[ _factor.size() - _numExtendedNodes ];
  int              startColumn = firstNode.startColumn();
  int              nodeStartColumn;
  int              startRow;

  const Real      *baseData;

  // Also write start and end columns for each block
  IntArray         startColumns;
  IntArray         endColumns;

  char             buf[ 1024 ];

  nCols = factorSystemSize() - startColumn;

  SPARSE_MATRIX    dataMatrix( nCols, nCols );

  for ( int node_idx = _factor.size() - _numExtendedNodes;
        node_idx < _factor.size(); node_idx++ )
  {
    const Supernode &node = _factor[ node_idx ];

    nodeColumns = node.numColumns();
    nodeStartColumn = node.startColumn();
    nodeStartColumn -= startColumn;

    startColumns.push_back( nodeStartColumn );
    endColumns.push_back( node.endColumn() - startColumn );

    baseData = _slackData + node.extendedDataOffset();

    for ( int i = 0; i < nodeColumns; i++ )
    for ( int j = 0; j <= i; j++ )
    {
      dataMatrix( nodeStartColumn + i, nodeStartColumn + j )
        = MATRIX::access( baseData, nodeColumns, nodeColumns, i, j );
    }

    for ( int interaction_idx = 0; interaction_idx < node.offDiagonal().size();
          interaction_idx++ )
    {
      const SupernodeInteraction &interaction
                                  = node.offDiagonal()[ interaction_idx ];
  
      startRow = _factor[ interaction._nodeID ].startColumn();
      startRow -= startColumn;

      baseData = _slackData + interaction._extendedDataOffset;

      for ( int i = 0; i < interaction._numExtendedRows; i++ )
      for ( int j = 0; j < nodeColumns; j++ )
      {
        dataMatrix( startRow + i, nodeStartColumn + j )
          = MATRIX::access( baseData, interaction._numExtendedRows, nodeColumns,
                            i, j );
      }
    }
  }

  printf( "Writing sparse matrix with size %d by %d and name %s\n",
          dataMatrix.rows(), dataMatrix.cols(), filename );
  dataMatrix.writeToBinary( filename );

  sprintf( buf, "%s_startColumns.vector", filename );
  writeVector( startColumns, buf );

  sprintf( buf, "%s_endColumns.vector", filename );
  writeVector( endColumns, buf );
}

//////////////////////////////////////////////////////////////////////
// For an extended block factored via an LDL factorization, writes
// the diagonal and permutation matrices to disk
//////////////////////////////////////////////////////////////////////
void FactorManager::writeExtendedLDLData( const char *filename )
{
  if ( _numExtendedNodes == 0 )
    return;

  int              nCols;
  int              nodeColumns;
  const Supernode &firstNode = _factor[ _factor.size() - _numExtendedNodes ];
  int              startColumn = firstNode.startColumn();
  int              nodeStartColumn;
  int              startRow;

  char             buf[ 1024 ];

  nCols = factorSystemSize() - startColumn;

  SPARSE_MATRIX    diagonalMatrix( nCols, nCols );
  SPARSE_MATRIX    permutationMatrix( nCols, nCols );

  for ( int node_idx = _factor.size() - _numExtendedNodes;
        node_idx < _factor.size(); node_idx++ )
  {
    const Supernode &node = _factor[ node_idx ];

    nodeColumns = node.numColumns();
    nodeStartColumn = node.startColumn();
    nodeStartColumn -= startColumn;

    const LDL_Data          &ldlData = node.ldlData();

    // Get the permutation
    {
      MATRIX                 P( nodeColumns, nodeColumns );
      IntArray               permutation;

      MATRIX::buildLDLPermutation( ldlData._pivotData, permutation, P.data() );

      // Copy to the sparse matrix
      for ( int i = 0; i < nodeColumns; i++ )
      for ( int j = 0; j < nodeColumns; j++ )
      {
        if ( P( i, j ) > 0.0 ) {
          permutationMatrix( nodeStartColumn + i,
                             nodeStartColumn + j ) = P( i, j );
        }
      }
    }

    // Add to the diagonal matrix
    {
      const IntArray        &pivotData = ldlData._pivotData;

      int                    blockIndex_1x1 = 0;
      int                    blockIndex_2x2 = 0;

      for ( int block_start = 0; block_start < nodeColumns; block_start++ )
      {
        if ( pivotData[ block_start ] < 0 ) {
          // Copy a 2x2 block
          const LDL_Data::Block_2x2 &block
                            = ldlData._diagonalBlockEntries[ blockIndex_2x2 ];

          for ( int i = 0; i < 2; i++ )
          for ( int j = 0; j < 2; j++ )
          {
            diagonalMatrix( nodeStartColumn + block_start + i,
                            nodeStartColumn + block_start + j )
                                  = MATRIX::access( block.data, 2, 2, i, j );
          }

          block_start += 1;
          blockIndex_2x2 += 1;
        }
        else {
          // Copy a 1x1 block
          const LDL_Data::Block_1x1 &block
                            = ldlData._diagonalEntries[ blockIndex_1x1 ];

          diagonalMatrix( nodeStartColumn + block_start,
                          nodeStartColumn + block_start ) = block;

          blockIndex_1x1 += 1;
        }
      }
    }
  }

  sprintf( buf, "D_%s", filename );
  diagonalMatrix.writeToBinary( buf );

  sprintf( buf, "P_%s", filename );
  permutationMatrix.writeToBinary( buf );
}

//////////////////////////////////////////////////////////////////////
// Prints out usage statistics for both off diagonal extended
// interactions occuring in the main part of the factorization,
// and those occuring in the extended variable space
//////////////////////////////////////////////////////////////////////
void FactorManager::printExtendedUsage()
{
  map<string, long int>      extendedUsage;

  for ( int node_idx = 0; node_idx < _factor.size(); node_idx++ )
  {
    const Supernode         &node = _factor[ node_idx ];
    int                      nCols = node.numColumns();

    if ( node.type() == STANDARD_NODE )
    {
      int                    blockSize = (int)log2( nCols );
      char                   blockName[ 1024 ];

      blockSize = (int)pow( 2.0, (Real)blockSize );
      sprintf( blockName, "Compression usage [%d, %d)",
               blockSize, 2 * blockSize );

      for ( int interaction_idx = node.firstExtendedInteraction();
            interaction_idx < node.offDiagonal().size(); interaction_idx++ )
      {
        const SupernodeInteraction  &interaction
                                      = node.offDiagonal()[ interaction_idx ];

        if ( interaction._compressed )
        {
          extendedUsage[ "Standard usage" ]
            += interaction._compressedColumnList.size()
               * interaction._numExtendedRows * sizeof( Real );
        }
        else
        {
          extendedUsage[ "Standard usage" ]
            += nCols * interaction._numExtendedRows * sizeof( Real );
        }

        // Determine compression cost
        if ( interaction._extendedInteraction >= 0
          && interaction._firstInteraction
          && !interaction._lowRankDiagonalInteraction )
        {
          const SupernodeInteraction &originalInteraction
                  = node.offDiagonal().at( interaction._extendedInteraction );
          const SupernodeInteraction &nextInteraction
            = _factor[ originalInteraction._nodeID ].
                offDiagonal()[ interaction._forwardInteractions[ 0 ].second ];

          if ( interaction._compressed )
          {
            extendedUsage[ blockName ]
              += interaction._numExtendedRows
                 * ( nCols + originalInteraction._rowList.size() )
                   * sizeof( Real );

            TRACE_ASSERT( NULL, "Should never get here!" );
          }
          else
          {
            if ( nextInteraction._compressed )
            {
              extendedUsage[ blockName ]
                += interaction._numExtendedRows
                   * ( nextInteraction._compressedColumnList.size()
                       + nCols ) * sizeof( Real );
            }
            else
            {
              extendedUsage[ blockName ]
                += interaction._numExtendedRows
                   * ( _factor[ originalInteraction._nodeID ].numColumns()
                       + nCols ) * sizeof( Real );
            }
          }
        }

        if ( interaction._lowRankDiagonalInteraction )
        {
          extendedUsage[ "Diagonal compression" ]
            += interaction._numExtendedRows
               * range_size( interaction._compressedColumnRange )
                 * sizeof( Real );
        }
      }
    }
    else
    {
      extendedUsage[ "Extended usage" ] += nCols * nCols * sizeof( Real );

      for ( int interaction_idx = 0;
            interaction_idx < node.offDiagonal().size(); interaction_idx++ )
      {
        const SupernodeInteraction  &interaction
                                      = node.offDiagonal()[ interaction_idx ];

        extendedUsage[ "Extended usage" ]
          += nCols * interaction._numExtendedRows * sizeof( Real );
      }
    }
  }

  printUsageMap( extendedUsage, "Extended variable memory usage" );
}

//////////////////////////////////////////////////////////////////////
// Estimates the storage used to just store the diagonal
// parts of the factor between nodes with sparse blocks
//////////////////////////////////////////////////////////////////////
void FactorManager::printUncompressedDiagonalStorage()
{
  int                        start_node = 0;
  int                        end_node = 0;

  long int                   totalSize = 0;

  int                        firstInteraction;

  int                        nCols;

  for ( ; start_node < _factor.size(); start_node = end_node )
  {
    // Find the next starting block of an uncompressed range
    for ( ; start_node < _factor.size(); start_node++ )
    {
      firstInteraction = _factor[ start_node ].firstExtendedInteraction();

      if ( firstInteraction >= _factor[ start_node ].offDiagonal().size() )
      {
        // Must be a normal block
        break;
      }
    }

    // Find the ending block
    end_node = start_node;
    for ( ; end_node < _factor.size(); end_node++ )
    {
      firstInteraction = _factor[ end_node ].firstExtendedInteraction();

      if ( firstInteraction < _factor[ end_node ].offDiagonal().size() )
      {
        // Must have extended variables
        break;
      }
    }

    // Count the sizes in this range
    for ( int node_idx = start_node; node_idx < end_node; node_idx++ )
    {
      const Supernode &node = _factor[ node_idx ];

      TRACE_ASSERT(
        node.firstExtendedInteraction() >= node.offDiagonal().size() );

      nCols = node.numColumns();

      totalSize += nCols * nCols;

      for ( int interaction_idx = 0;
            interaction_idx < node.offDiagonal().size();
            interaction_idx++ )
      {
        const SupernodeInteraction &interaction
                                      = node.offDiagonal()[ interaction_idx ];

        if ( interaction._nodeID >= start_node
          && interaction._nodeID < end_node )
        {
          totalSize += interaction._rowList.size() * nCols;
        }
      }
    }
  }

  printf( "Diagonal blocks require %lld non-zeros = %f MB\n",
          (long long int)totalSize, (Real)totalSize * 8.0 / 1024.0 / 1024.0 );
}
 
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void FactorManager::printInteriorBlocks() const
{
  for ( int block_idx = 0; block_idx < _interiorBlocks.size(); block_idx++ )
  {
    const InteriorBlock     &block = _interiorBlocks[ block_idx ];

    printf( "Interior block %d: Node range [ %d, %d ]\n",
            block_idx, block._nodeRange.first, block._nodeRange.second );
  }
}

//////////////////////////////////////////////////////////////////////
// Determines the storage used to just store the diagonal parts
// of interior blocks
//////////////////////////////////////////////////////////////////////
void FactorManager::printInteriorBlockDiagonalUsage() const
{
  size_t                     totalSize = 0;
  int                        nCols;
  int                        nRows;
  int                        start_node;
  int                        end_node;

  for ( int block_idx = 0; block_idx < _interiorBlocks.size(); block_idx++ ) {
    start_node = _interiorBlocks[ block_idx ]._nodeRange.first;
    end_node = _interiorBlocks[ block_idx ]._nodeRange.second;

    printf( "Interior block range %d to %d\n", start_node, end_node );

    for ( int node_idx = start_node; node_idx <= end_node; node_idx++ ) {
      const Supernode       &node = _factor[ node_idx ];

      nCols = node.numColumns();

      totalSize += nCols * nCols;

      for ( int interaction_idx = 0;
            interaction_idx < node.offDiagonal().size(); interaction_idx++ )
      {
        const SupernodeInteraction &interaction
                                      = node.offDiagonal()[ interaction_idx ];

        if ( interaction._nodeID < start_node
          || interaction._nodeID > end_node )
        {
          continue;
        }

        nRows = interaction._rowList.size();

        totalSize += nCols * nRows;
      }
    }
  }

  printf( "Interior block diagonals require %f MB\n",
          (Real)totalSize * sizeof( Real ) / 1024.0 / 1024.0 );
}

//////////////////////////////////////////////////////////////////////
// Extract a sub-matrix of the factor
//////////////////////////////////////////////////////////////////////
void FactorManager::extractSubMatrix( const IntArray &nodeList,
                                      MATRIX &subMatrix ) const
{
  int                        totalColumns = 0;

  for ( int node_idx = 0; node_idx < nodeList.size(); node_idx++ )
  {
    totalColumns += nodeColumns( nodeList[ node_idx ] );
  }

  subMatrix.resizeAndWipe( totalColumns, totalColumns );

  for ( int node_idx = 0; node_idx < nodeList.size(); node_idx++ )
  {
    _factor.at( nodeList[ node_idx ] ).extractSubMatrix( nodeList,
                                                         _factor,
                                                         subMatrix.data(),
                                                         totalColumns );
  }
}

//////////////////////////////////////////////////////////////////////
// Extract a sub-matrix of the factor in sparse format
//////////////////////////////////////////////////////////////////////
void FactorManager::extractSubMatrix( const IntArray &nodeList,
                                      SPARSE_MATRIX &subMatrix ) const
{
  int                        totalColumns = 0;

  for ( int node_idx = 0; node_idx < nodeList.size(); node_idx++ )
  {
    totalColumns += nodeColumns( nodeList[ node_idx ] );
  }

  subMatrix.resize( totalColumns, totalColumns );
  subMatrix.clear();


}

//////////////////////////////////////////////////////////////////////
// Extracts an interior block interaction matrix
//////////////////////////////////////////////////////////////////////
void FactorManager::extractInteriorBlockInteraction( int interior_block_idx,
                                                     int interaction_idx,
                                                     MATRIX &subMatrix )
{
  const InteriorBlock       &block = _interiorBlocks[ interior_block_idx ];

  const InteriorBlock::Interaction
                          &interaction = block._interactions[ interaction_idx ];

  int                        startColumn;
  int                        fullColumnCount;
  int                        numRows;
  int                        interNodeID;

  const Real                *interactionData;
  Real                      *baseData;

  fullColumnCount = _factor[ block._nodeRange.second ].endColumn()
                  - _factor[ block._nodeRange.first ].startColumn() + 1;

  TRACE_ASSERT( interaction.size() > 0 );
  {
    const Supernode         &baseNode = _factor[ interaction[ 0 ].first ];

    interNodeID = baseNode.offDiagonal()[ interaction[ 0 ].second ]._nodeID;

    numRows = _factor[ interNodeID ].numColumns();
  }

  subMatrix.resizeAndWipe( numRows, fullColumnCount );

  for ( int interaction_idx = 0; interaction_idx < interaction.size();
        interaction_idx++ )
  {
    int nodeID = interaction[ interaction_idx ].first;
    int interID = interaction[ interaction_idx ].second;
    int nCols = _factor[ nodeID ].numColumns();
    int nodeStartCol = _factor[ nodeID ].startColumn() - startColumn;

    TRACE_ASSERT(
      _factor[ nodeID ].offDiagonal()[ interID ]._nodeID == interNodeID );

    const IntArray &rowList
                  = _factor[ nodeID ].offDiagonal()[ interID ]._rowList;

    // Move to the right column
    baseData = subMatrix.data() + nodeStartCol;

    interactionData = _factor[ nodeID ].interactionMatrix( interID );

    for ( int row_idx = 0; row_idx < rowList.size(); row_idx++ )
    {
      MATRIX::copyRow( interactionData, baseData,
                       // Copy row row_idx to rowList[ row_idx ]
                       row_idx, rowList[ row_idx ],
                       // Input has nCols columns
                       nCols,
                       // Output has fullColumnCount columns
                       fullColumnCount,
                       // Only copy nCols columns
                       nCols );
    }
  }
}

//////////////////////////////////////////////////////////////////////
// Helper function for factorNumeric
//
// Initializes slack variable data for all extended nodes
//////////////////////////////////////////////////////////////////////
void FactorManager::initExtendedNodes( int start_node )
{
  for ( int node_idx = start_node; node_idx < _factor.size(); node_idx++ )
  {
    Supernode               &node = _factor[ node_idx ];

    TRACE_ASSERT( node.type() == EXTENDED_NODE );

    node.initializeExtendedNode( _factor, node_idx,
                                 _slackData, _slackDataOffset,
                                 _availableSlackData );

    TRACE_ASSERT( _availableSlackData >= 0, "Ran out of slack variable space" );
  }
}

//////////////////////////////////////////////////////////////////////
// Helper function for factorNumeric
//
// Builds the schur complement in extended variable space
//////////////////////////////////////////////////////////////////////
void FactorManager::buildExtendedSchurComplement( int numStandardNodes )
{
  for ( int node_idx = 0; node_idx < numStandardNodes; node_idx++ )
  {
    Supernode               &node = _factor[ node_idx ];

    TRACE_ASSERT( node.type() == STANDARD_NODE );

#if 0
    node.addExtendedSchurComplementContribution(
                    _factor, _expansionWorkspace, _copyWorkspace, _slackData,
                    &_timers[ EXTENDED_SCHUR_INVERSION ],
                    &_timers[ EXTENDED_SCHUR_MULTIPLY ] );
#endif
    node.addExtendedSchurComplementContribution(
                    _factor, _realWorkspaceManager, _slackData,
#ifdef DO_TIMING
                    &_timers[ EXTENDED_SCHUR_INVERSION ],
                    &_timers[ EXTENDED_SCHUR_MULTIPLY ]
#else
                    NULL, NULL
#endif
                    );
  }
}

//////////////////////////////////////////////////////////////////////
// Initializes workspace data according to the size of the system
// to be factored.
//////////////////////////////////////////////////////////////////////
void FactorManager::initWorkspace()
{
  int                        rowsMax = 0;
  int                        interactionMax = 0;
  int                        nCols, nRows;

  size_t                     blockSizeMax = 0;
  size_t                     blockSize;

  long int                   floatWorkTotalSize = 0;
  long int                   intWorkTotalSize = 0;

  map<string, long int>      floatWorkUsage;

  // Figure out which supernode has the largest number of rows, and
  // the largest block size (columns * rows)
  for ( int node_idx = 0; node_idx < _factor.size(); node_idx++ )
  {
    nCols = _factor[ node_idx ].numColumns();
    nRows = _factor[ node_idx ].numRows();
    
    // Note: nRows only refers to the number of off-diagonal rows

    blockSize = _factor[ node_idx ].uncompressedEntries();

    rowsMax = max( rowsMax, nCols + nRows );

    //blockSizeMax = max( blockSizeMax, nCols * ( nCols + nRows ) );
    blockSizeMax = max( blockSizeMax, blockSize );

    interactionMax = max( interactionMax,
                          (int)_factor[ node_idx ].offDiagonal().size() );
  }

  // Both _multWorkspace and _updateWorkspace can have size blockSizeMax,
  // though this is actually potentially larger than necessary for
  // _updateWorkspace
  _multWorkspace = new Real[ blockSizeMax ];
  _updateWorkspace = new Real[ blockSizeMax ];

  floatWorkTotalSize += blockSizeMax;
  floatWorkTotalSize += blockSizeMax;

  floatWorkUsage[ "mult. workspace" ] += blockSizeMax * sizeof( Real );
  floatWorkUsage[ "update workspace" ] += blockSizeMax * sizeof( Real );

  // FIXME: debugging
  Supernode::multWorkspaceSz = blockSizeMax;
  Supernode::workspaceSz = blockSizeMax;

  _multWorkspaceSz = blockSizeMax;
  _updateWorkspaceSz = blockSizeMax;

  // The scattered map size should match the number of unknowns in the
  // original system

  // FIXME: I don't know why I have to do this, but if I don't
  // I get a crash
  _scatteredMap.clear();
  for ( int i = 0; i < _A._nrow; i++ )
  {
    _scatteredMap.push_back( EMPTY );
  }
  _lowRankScatteredMap.resize( _A._nrow, EMPTY );

  intWorkTotalSize += 2 * _A._nrow;

  // The relative map should just have enough space allocated
  // so that it never needs to re-allocate
  _relativeMap.reserve( rowsMax );

  intWorkTotalSize += rowsMax;

  _interactionIndices.reserve( interactionMax );
  _lowRankInteractionIndices.reserve( interactionMax );
  _extendedInteractionIndices.reserve( interactionMax );

  intWorkTotalSize += 6 * interactionMax;

  _solveWorkspaceSz = _A._nrow;
  _solveWorkspace = new Real[ _A._nrow ];

  floatWorkTotalSize += _A._nrow;

  floatWorkUsage[ "solve workspace" ] += _A._nrow * sizeof( Real );

  _nextInteraction.resize( _factor.size(), EMPTY );
  _nextInteractionCache.resize( _factor.size(), EMPTY );
  _head.resize( _factor.size(), EMPTY );
  _next.resize( _factor.size(), EMPTY );

  intWorkTotalSize += 4 * _factor.size();

  printf( "Standard float workspace: %f MB\n",
          (Real)( floatWorkTotalSize * sizeof( Real ) ) / ( 1024.0 * 1024.0 ) );
  printf( "Standard integer workspace: %f MB\n\n",
          (Real)( intWorkTotalSize * sizeof( int ) ) / ( 1024.0 * 1024.0 ) );

  // Finally, allocate memory for the numerical factor itself
  _factorData = Supernode::AllocateFixedData( _factor, _factorDataSz );

  printUsageMap( floatWorkUsage, "Standard float workspace usage" );
}

//////////////////////////////////////////////////////////////////////
// Initialize workspaces for forming low-rank decompositions.
//
// TODO: For now, consider the fixed rank problem; that is,
//       we will determine the desired rank for a block using
//       the blockRank function (instead of actually trying to
//       approximate to some tolerance)
//////////////////////////////////////////////////////////////////////
void FactorManager::initDecompositionWorkspace()
{
  // Figure out the maximum rank and maximum number of
  // rows of any low-rank block to decompose.
  // This is needed for the QR workspace.
  // Also need the maximum possible size for a decomposition
  // workspace.
  int                        maxRank = 0;
  int                        maxRows = 0;
  int                        maxBlockSize = 0;
  int                        maxTransBlockSize = 0;

  int                        maxNodeSize = 0;

  int                        nRows, rank;

  IntArray                   expansionWorkSizes( _factor.size(), 0 );

  // For diagonal block decomposition workspaces
  IntArray                   descendentCounts;
  IntArray                   lowRankDiagonalBlockCounts;

  int                        maxDescendents;
  int                        maxDiagonalLowRankBlocks;

  int                        maxDiagonalRows = 0;
  int                        maxDiagonalColumns = 0;
  int                        maxDiagonalRank = 0;
  int                        maxVectorSize = 0;

  long int                   intWorkTotalSize = 0;
  long int                   floatWorkTotalSize = 0;

  map<string, long int>      floatWorkUsage;

  _lowRankBlocks.reserve( _A._nrow );
  _blockRanks.reserve( _A._nrow );
  _blockRanksCache.reserve( _A._nrow );
  _lowRankBaseRows.reserve( _A._nrow );
  _lowRankBlockActive.reserve( _A._nrow );

  intWorkTotalSize += 4 * _A._nrow * sizeof( int );
  intWorkTotalSize += _A._nrow * sizeof( bool );

  for ( int node_idx = 0; node_idx < _factor.size(); node_idx++ )
  {
    const Supernode         &node = _factor[ node_idx ];

    maxNodeSize = max( maxNodeSize, node.numColumns() );

    diagonalBlockSizes( node, maxDiagonalRows, maxDiagonalColumns,
                        maxDiagonalRank );

    if ( node.type() == EXTENDED_NODE || node.offDiagonal().size() == 0 )
    {
      continue;
    }

    node.getLowRankBlocks( _lowRankBlocks, _lowRankBaseRows );

    if ( _lowRankDescendents.size() < _lowRankBlocks.size() )
    {
      _lowRankDescendents.resize( _lowRankBlocks.size() );
    }

    if ( _lowRankStartRows.size() < _lowRankBlocks.size() )
    {
      _lowRankStartRows.resize( _lowRankBlocks.size() );
    }

    maxVectorSize = max( node.numColumns() * (int)_lowRankBlocks.size(),
                         maxVectorSize );

    nRows = node.countLowRankRows();

    maxRows = max( maxRows, nRows );

    rank = maxBlockRank( nRows, node.numColumns() ) + _overSampling;
    maxRank = max( maxRank, rank );

    maxBlockSize = max( maxBlockSize, maxRows * rank );

    maxTransBlockSize = max(
      maxTransBlockSize, (int)_lowRankBlocks.size() * node.numColumns() * rank );

    // Add to expansion workspace totals resulting from this node
    node.estimateCompressedExpansionStorage( _factor, rank,
                                             expansionWorkSizes );
  }

  printf( "Reserving %d entries in %d low rank descendent lists.\n",
          (int)_factor.size(), (int)_lowRankDescendents.size() );

  _extendedSolveWorkspaceSz = maxRank * _numExtendedNodes;
  if (_extendedSolveWorkspaceSz > 0)
    _extendedSolveWorkspace = new Real[ _extendedSolveWorkspaceSz ];

  floatWorkTotalSize += _extendedSolveWorkspaceSz * sizeof( Real );

  floatWorkUsage[ "extended solve" ]
    += _extendedSolveWorkspaceSz * sizeof( Real );

  for ( int i = 0; i < _lowRankDescendents.size(); i++ )
  {
    _lowRankDescendents[ i ].reserve( _factor.size() );
    _lowRankStartRows[ i ].reserve( maxNodeSize );

    intWorkTotalSize += _factor.size() * sizeof( IndexPair );
    intWorkTotalSize += maxNodeSize * sizeof( int );
  }

  _lowRankBlocks.clear();
  _lowRankBaseRows.clear();

  //////////////////////////////////////////////////////////////////////
  // Include diagonal blocks in all workspace size calculations
  //////////////////////////////////////////////////////////////////////
  maxRank = max( maxRank, maxDiagonalRank );
  maxTransBlockSize = max( maxTransBlockSize,
                           maxDiagonalColumns * maxDiagonalRank );
  maxRows = max( maxRows, maxDiagonalRows );
  maxRows = max( maxRows, maxDiagonalColumns );

#if 0
  // Initialize a maxNodeSize x maxRank random matrix
  _testMatrix = MathUtil::randomGaussianMatrix( maxNodeSize, maxRank );

  floatWorkTotalSize += maxNodeSize * maxRank * sizeof( Real );
  floatWorkUsage[ "test matrix" ] += maxNodeSize * maxRank * sizeof( Real );
#endif

#if 0
  _subMatrix = new Real[ maxNodeSize * maxRank ];

  floatWorkTotalSize += maxNodeSize * maxRank * sizeof( Real );
  floatWorkUsage[ "sub matrix" ] += maxNodeSize * maxRank * sizeof( Real );
#endif

  maxBlockSize = max( maxBlockSize, maxNodeSize * maxRank );

  // Decomposition work space should be as big as the largest
  // block size we wish to handle
  // FIXME: size change
  _decompWorkspaceSz = maxBlockSize;
  //_decompWorkspaceSz = maxBlockSize * 2;
  _decompWorkspace = new Real[ _decompWorkspaceSz ];

  floatWorkTotalSize += _decompWorkspaceSz * sizeof( Real );
  floatWorkUsage[ "decomp. workspace" ] += _decompWorkspaceSz * sizeof( Real );

#if 0
  _decompMultWorkspaceSz = max( maxBlockSize, maxNodeSize * maxRank );
  _decompMultWorkspace = new Real[ _decompMultWorkspaceSz ];

  floatWorkTotalSize += _decompMultWorkspaceSz * sizeof( Real );
  floatWorkUsage[ "decomp. mult. workspace" ]
    += _decompMultWorkspaceSz * sizeof( Real );
#endif

  //_decompMultTransWorkspaceSz = maxNodeSize * maxRank;
  _decompMultTransWorkspaceSz = maxTransBlockSize;
  if (_decompMultTransWorkspaceSz > 0)
    _decompMultTransWorkspace = new Real[ _decompMultTransWorkspaceSz ];

  floatWorkTotalSize += _decompMultTransWorkspaceSz * sizeof( Real );
  floatWorkUsage[ "decomp. mult. trans. workspace" ]
    += _decompMultTransWorkspaceSz * sizeof( Real );

#if 0
  // I *think* this is big enough for the copy workspace
  _copyWorkspaceSz = maxNodeSize * maxRank;
  _copyWorkspace = new Real[ _copyWorkspaceSz ];

  floatWorkTotalSize += _copyWorkspaceSz * sizeof( Real );
  floatWorkUsage[ "copy workspace" ] += _copyWorkspaceSz * sizeof( Real );
#endif

  _qrExtraDataSz = max( maxRows, maxRank );
  _qrExtraData = new Real[ _qrExtraDataSz ];

  floatWorkTotalSize += _qrExtraDataSz * sizeof( Real );
  floatWorkUsage[ "QR extra data" ] += _qrExtraDataSz * sizeof( Real );

  // I *think* this is enough space
  _qrWorkspaceSz = maxRank * 64;
  _qrWorkspace = new Real[ _qrWorkspaceSz ];

  floatWorkTotalSize += _qrWorkspaceSz * sizeof( Real );
  floatWorkUsage[ "QR workspace" ] += _qrWorkspaceSz * sizeof( Real );

  //////////////////////////////////////////////////////////////////////
  // Initialize workspaces for adaptive low rank decomposition
  //////////////////////////////////////////////////////////////////////

  //_blockWorkspaceSz = MAX_BLOCK_SIZE * maxRank;
#if 0
  _blockWorkspaceSz = maxRank * maxRank;
  _blockWorkspace = new Real[ _blockWorkspaceSz ];

  floatWorkTotalSize += _blockWorkspaceSz * sizeof( Real );
  floatWorkUsage[ "block workspace" ] += _blockWorkspaceSz * sizeof( Real );
#endif

  _vectorWorkSz = max( maxRows, maxRank * (int)_lowRankDescendents.size() );
  _vectorWorkSz = max( _vectorWorkSz, maxVectorSize );
  _vectorWork1Sz = _vectorWorkSz;
  _vectorWork2Sz = _vectorWorkSz;
  if (_vectorWorkSz > 0) {
    _vectorWork1 = new Real[ _vectorWorkSz ];
    _vectorWork2 = new Real[ _vectorWorkSz ];
  }

  floatWorkTotalSize += 2 * _vectorWorkSz * sizeof( Real );
  floatWorkUsage[ "vector workspace" ] += 2 * _vectorWorkSz * sizeof( Real );

  //_basisWorkspaceSz = maxRows * MAX_BLOCK_SIZE;
  // FIXME: size change
  _basisWorkspaceSz = maxRows * maxRank;
  //_basisWorkspaceSz = maxRows * maxRank * 2;
  if ( _basisWorkspace != NULL )
    delete[] _basisWorkspace;
  if ( _basisWorkspaceSz > 0 )
    _basisWorkspace = new Real[ _basisWorkspaceSz ];

  floatWorkTotalSize += _basisWorkspaceSz * sizeof( Real );
  floatWorkUsage[ "basis workspace" ] += _basisWorkspaceSz * sizeof( Real );

  // Expansion workspace should be large enough to store the fully
  // expanded versions of all low rank bases in a node.
#if 0
  _expansionWorkspaceSz = 0;
  for ( int node_idx = 0; node_idx < _factor.size(); node_idx++ )
  {
    _expansionWorkspaceSz = max( _expansionWorkspaceSz,
                                 (long int)expansionWorkSizes[ node_idx ] );
  }
  _expansionWorkspaceSz = max( 2 * _copyWorkspaceSz, _expansionWorkspaceSz );
  _expansionWorkspace = new Real[ _expansionWorkspaceSz ];

  floatWorkTotalSize += _expansionWorkspaceSz * sizeof( Real );
  floatWorkUsage[ "expansion workspace" ]
    += _expansionWorkspaceSz * sizeof( Real );
#endif

  // FIXME: size change
  _basisStorage._dataSz = maxBlockSize;
  //_basisStorage._dataSz = maxBlockSize * 1.5;
  _basisStorage._data = new Real[ _basisStorage._dataSz ];
  _basisStorage._dataPtrs.reserve( _lowRankDescendents.size() );

  floatWorkTotalSize += _basisStorage._dataSz * sizeof( Real );
  floatWorkUsage[ "basis storage" ] += _basisStorage._dataSz * sizeof( Real );
  intWorkTotalSize += _lowRankDescendents.size() * sizeof( BlockStorageData );

  _blockErrors.reserve( _A._nrow );

  floatWorkTotalSize += _A._nrow * sizeof( VEC3F );
  floatWorkUsage[ "block errors" ] += _A._nrow * sizeof( VEC3F );

  //////////////////////////////////////////////////////////////////////
  // Initialize workspaces needed for decomposition of low rank blocks
  // in the diagonal
  //////////////////////////////////////////////////////////////////////
  countDescendents( descendentCounts );
  countLowRankDiagonalBlocks( lowRankDiagonalBlockCounts );

  maxDescendents = maxEntry( descendentCounts );
  maxDiagonalLowRankBlocks = maxEntry( lowRankDiagonalBlockCounts );

  _currentDescendents.reserve( maxDescendents );
  _inverseRowMap.reserve( _A._nrow );
  _diagonalLowRankBlockActive.reserve( maxDescendents );

  intWorkTotalSize += maxDescendents * sizeof( int );

  for ( int i = 0; i < maxDescendents; i++ )
  {
    _diagonalBlockRanges.push_back( vector<DenseBlock>() );
    _diagonalBlockRanges.back().reserve( maxDiagonalLowRankBlocks );

    intWorkTotalSize += maxDiagonalLowRankBlocks * sizeof( DenseBlock );
  }

  //////////////////////////////////////////////////////////////////////
  // Clear everything
  //////////////////////////////////////////////////////////////////////

#if 0
  memset( (void *)_subMatrix, 0, maxNodeSize * maxRank * sizeof( Real ) );
#endif
  memset( (void *)_decompWorkspace, 0, _decompWorkspaceSz * sizeof( Real ) );
#if 0
  memset( (void *)_decompMultWorkspace, 0,
          _decompMultWorkspaceSz * sizeof( Real ) );
#endif
  memset( (void *)_decompMultTransWorkspace, 0,
          _decompMultTransWorkspaceSz * sizeof( Real ) );
#if 0
  memset( (void *)_copyWorkspace, 0, _copyWorkspaceSz * sizeof( Real ) );
#endif
  memset( (void *)_qrExtraData, 0, _qrExtraDataSz * sizeof( Real ) );
  memset( (void *)_qrWorkspace, 0, _qrWorkspaceSz * sizeof( Real ) );

#if 0
  memset( (void *)_blockWorkspace, 0, _blockWorkspaceSz * sizeof( Real ) );
#endif
  memset( (void *)_vectorWork1, 0, _vectorWorkSz * sizeof( Real ) );
  memset( (void *)_vectorWork2, 0, _vectorWorkSz * sizeof( Real ) );
  memset( (void *)_basisWorkspace, 0, _basisWorkspaceSz * sizeof( Real ) );
  memset( (void *)_basisStorage._data, 0,
          _basisStorage._dataSz * sizeof( Real ) );

  // Get a size for our slack variable data
  _slackDataSz = estimateExtraStorage();

  Real dataSizeMB = (Real)( _slackDataSz * 8 ) / 1024.0 / 1024.0;

  printf( "\nRequesting %lld slack variable nonzeros, %f MB of data\n\n",
          (long long int)_slackDataSz, dataSizeMB );
  printf( "Extended float workspace: %f MB\n",
          (Real)( floatWorkTotalSize ) / ( 1024.0 * 1024.0 ) );
  printf( "Extended int workspace: %f MB\n\n",
          (Real)( intWorkTotalSize ) / ( 1024.0 * 1024.0 ) );

  printUsageMap( floatWorkUsage, "Extended float workspace usage" );
}

//////////////////////////////////////////////////////////////////////
// Walks through the structure of a numerical factorization to
// count the number of descendents each node expects to encounter
//////////////////////////////////////////////////////////////////////
void FactorManager::countDescendents( IntArray &descendentCounts ) const
{
  // Temporary arrays to mirror those used for the actual numerical
  // factorization
  IntArray                   head( _factor.size(), EMPTY );
  IntArray                   next( _factor.size(), EMPTY );
  IntArray                   nextInteraction( _factor.size(), EMPTY );

  int                        descendent_idx;
  int                        descendent_idx_next;
  int                        start_idx;
  int                        descendentAncestor;
  int                        parentNode;

  descendentCounts.clear();
  descendentCounts.resize( _factor.size(), 0 );

  for ( int node_idx = 0; node_idx < _factor.size(); node_idx++ )
  {
    const Supernode         &node = _factor[ node_idx ];

    for ( descendent_idx = head[ node_idx ]; descendent_idx != EMPTY;
          descendent_idx = descendent_idx_next )
    {
      const Supernode       &descendent = _factor[ descendent_idx ];

      descendentCounts[ node_idx ]++;

      descendent_idx_next = next[ descendent_idx ];

      start_idx = nextInteraction[ descendent_idx ];

      start_idx++;
      nextInteraction[ descendent_idx ] = start_idx;

      // Increment the start index until we find an active node.
      // ie. skip any inactive nodes, since they no longer matter
      // for the factorization.
      while ( start_idx < descendent.offDiagonal().size()
           && !descendent.offDiagonal()[ start_idx ]._active )
      {
        start_idx++;
        nextInteraction[ descendent_idx ] = start_idx;
      }

      if ( start_idx < descendent.offDiagonal().size() )
      {
        descendentAncestor = descendent.offDiagonal()[ start_idx ]._nodeID;

        TRACE_ASSERT( descendentAncestor > node_idx
                   && descendentAncestor < _factor.size(),
                   "Invalid descendent ancestor" );

        next[ descendent_idx ] = head[ descendentAncestor ];
        head[ descendentAncestor ] = descendent_idx;
      }
    }

    // Find the first interaction which is not low rank
    start_idx = 0;

    while ( start_idx < node.offDiagonal().size()
         && !node.offDiagonal()[ start_idx ]._active )
    {
      start_idx++;
    }

    // Update descendent relationships
    if ( start_idx < node.offDiagonal().size() )
    {
      nextInteraction[ node_idx ] = start_idx;

      parentNode = node.offDiagonal()[ start_idx ]._nodeID;

      next[ node_idx ] = head[ parentNode ];
      head[ parentNode ] = node_idx;
    }

    // Link list for supernode s no longer needed (?)
    head[ node_idx ] = EMPTY;
  }
}

//////////////////////////////////////////////////////////////////////
// Counts the number of low rank block stored in the main diagonal
// of each node
//////////////////////////////////////////////////////////////////////
void FactorManager::countLowRankDiagonalBlocks( IntArray &blockCounts ) const
{
  blockCounts.clear();
  blockCounts.resize( _factor.size(), 0 );

  for ( int node_idx = 0; node_idx < _factor.size(); node_idx++ )
  {
    const Supernode         &node = _factor[ node_idx ];

    const vector<DenseBlock> &lowRankBlocks = node.diagonalLowRankBlocks();

    blockCounts[ node_idx ] = lowRankBlocks.size();
  }
}

#if 0
//////////////////////////////////////////////////////////////////////
// Estimates the total extra data needed for decomposing
// low-rank blocks in the factor
//////////////////////////////////////////////////////////////////////
long int FactorManager::estimateExtraStorage()
{
  long int                   totalStorage = 0;
  int                        rank;
  int                        nCols;

  // For each node
  //    1) Look at each low rank block
  //    2) Estimate its maximum rank
  //    3) Based on this rank, figure out how much additional
  //       fill-in will be generated amongst extended nodes.
  //
  for ( int node_idx = 0; node_idx < _factor.size(); node_idx++ )
  {
    Supernode               &node = _factor[ node_idx ];

    if ( node.type() == EXTENDED_NODE )
    {
      continue;
    }

    for ( int interaction_idx = 0; interaction_idx < node.offDiagonal().size();
          interaction_idx++ )
    {
      const SupernodeInteraction &interaction
                            = node.offDiagonal()[ interaction_idx ];

      if ( interaction._active )
      {
        continue;
      }

      const SupernodeInteraction &extendedInteraction
                      = node.offDiagonal()[ interaction._extendedInteraction ];

      rank = blockRank( interaction._rowList.size(), node.numColumns() );

      // This will generate a new node, which will have diagonal
      // size rank x rank.  It will also introduce off-diagonal
      // contributions to this node of size rank x numColumns
      if ( extendedInteraction._compressed )
      {
        totalStorage
          += rank * ( rank + extendedInteraction._compressedColumnList.size() );
      }
      else
      {
        totalStorage += rank * ( rank + node.numColumns() );
      }

      // Temporarily assign a column range to the node introduced
      // by this interaction
      Supernode &newNode = _factor[ node_idx + interaction._extendedOffset ];
      Supernode &prevNode = _factor[ node_idx + interaction._extendedOffset-1 ];

      TRACE_ASSERT( newNode.type() == EXTENDED_NODE, "Extended node expected" );
      TRACE_ASSERT( newNode.numColumns() == 0, "Empty node expected" );

      newNode.setColumnRange( prevNode.endColumn() + 1,
                              prevNode.endColumn() + rank );

      TRACE_ASSERT( interaction._extendedInteraction >= 0,
                    "No extended interaction index" );

      // Fill-in
      TRACE_ASSERT( extendedInteraction._firstInteraction,
                    "Not the first interaction" );

      const PairArray &forwardInteractions
        = node.offDiagonal()[
                    interaction._extendedInteraction ]._forwardInteractions;

      // We will have fill-in contributions at each node in this list
      for ( int fill_idx = 0; fill_idx < forwardInteractions.size();
            fill_idx++ )
      {
        const Supernode &nextNode
          = _factor[ forwardInteractions[ fill_idx ].first ];
        const SupernodeInteraction &nextInteraction
          = nextNode.offDiagonal()[ forwardInteractions[ fill_idx ].second ];

        if ( nextInteraction._compressed )
        {
          nCols = nextInteraction._compressedColumnList.size();
          TRACE_ASSERT( nCols > 0, "No column range set for this node" );
        }
        else
        {
          nCols = nextNode.numColumns();
          TRACE_ASSERT( nCols > 0, "Node has zero columns" );
        }

        totalStorage += rank * nCols;
      }
    }

    // Estimate storage required for compressed blocks in the diagonal
    for ( int block_idx = 0; block_idx < node.diagonalLowRankBlocks().size();
          block_idx++ )
    {
      int                    interaction_idx;
      
      interaction_idx = node.lowRankDiagonalInteraction( block_idx );

      const DenseBlock      &block = node.diagonalLowRankBlocks()[ block_idx ];

      const SupernodeInteraction &interaction
                            = node.offDiagonal()[ interaction_idx ];

      // Estimate the rank
      rank = blockRank( block.numRows(), block.numColumns() );

      TRACE_ASSERT( interaction._compressed,
                    "Diagonal compression should produce "
                    "compressed interactions" );

      nCols = interaction._compressedColumnList.size();

      TRACE_ASSERT( nCols > 0, "Empty column list" );

      totalStorage += rank * ( rank + nCols );

#if 0
      printf( "Estimating extra storage resulting from block %d of "
              "node %d\nNew node ID = %d\n\n", block_idx, node_idx,
              interaction._nodeID );
#endif

      // Temporarily assign a column range to the node introduced
      // by this interaction
      Supernode &newNode = _factor[ interaction._nodeID ];
      Supernode &prevNode = _factor[ interaction._nodeID - 1 ];

      TRACE_ASSERT( newNode.type() == EXTENDED_NODE, "Extended node expected" );
      TRACE_ASSERT( newNode.numColumns() == 0, "Empty node expected" );

      newNode.setColumnRange( prevNode.endColumn() + 1,
                              prevNode.endColumn() + rank );

      // Fill-in
      const PairArray &forwardInteractions = interaction._forwardInteractions;

      // We will have fill-in contributions at each node in this list
      for ( int fill_idx = 0; fill_idx < forwardInteractions.size();
            fill_idx++ )
      {
        const Supernode &nextNode
          = _factor[ forwardInteractions[ fill_idx ].first ];
        const SupernodeInteraction &nextInteraction
          = nextNode.offDiagonal()[ forwardInteractions[ fill_idx ].second ];

        if ( nextInteraction._compressed )
        {
          nCols = nextInteraction._compressedColumnList.size();
        }
        else
        {
          nCols = nextNode.numColumns();
        }

        TRACE_ASSERT( nCols > 0, "No column range set for this node" );

        totalStorage += rank * nCols;
      }
    }
  }

  // Test out extended system ordering here...
  IntArray permutation;
  IntArray inversePermuation;
  buildExtendedPermutation( permutation, inversePermuation );

#if 0
  writeExtendedSchurComplement( "super_numeric/extended_schur_initial.bcsm" );
#endif

  // Undo any changes to the factor
  for ( int node_idx = 0; node_idx < _factor.size(); node_idx++ )
  {
    if ( _factor[ node_idx ].type() == EXTENDED_NODE )
    {
      _factor[ node_idx ].setColumnRange( -1, -2 ); // So numColumns() == 0
    }
  }

  // FIXME
  //totalStorage /= 7;

  return totalStorage;
}
#endif

//////////////////////////////////////////////////////////////////////
// Estimates the total extra data needed for decomposing
// low-rank blocks in the factor
//////////////////////////////////////////////////////////////////////
long int FactorManager::estimateExtraStorage()
{
  long int                   totalStorage = 0;
  int                        rank;
  int                        nCols;
  int                        nRows;

  long int                   storageEstimate = 0;

  // Estimate sizes for all slack variable blocks
  estimateSlackVariableSizes();

  // Add up the total size of all extended interactions and nodes
  // given the estimated slack variable sizes
  for ( int node_idx = 0; node_idx < _factor.size(); node_idx++ )
  {
    const Supernode         &node = _factor[ node_idx ];

    if ( node.type() == STANDARD_NODE )
    {
      nCols = node.numColumns();
    }
    else
    {
      nCols = node.numExtendedColumns();

      TRACE_ASSERT( nCols > 0, "Empty extended node" );

      // Add the diagonal contribution
      totalStorage += nCols * nCols;
    }

    // Accumulate storage for diagonal blocks that are compressed
    // in-place
    if ( node.lowRankDiagonal() && node.inPlaceDiagonal() ) {
      for ( int block_idx = 0; block_idx < node.diagonalLowRankBlocks().size();
            block_idx++ )
      {
        const DenseBlock    &block = node.diagonalLowRankBlocks()[ block_idx ];

        rank = diagonalBlockRank( block.numRows(), block.numColumns() );

        totalStorage += rank * ( block.numRows() + block.numColumns() );
      }
    }

    // FIXME: debugging an idea here
    int                      maxRank = 0;
    int                      totalRank = 0;

    for ( int interaction_idx = 0; interaction_idx < node.offDiagonal().size();
          interaction_idx++ )
    {
      const SupernodeInteraction &interaction
                                    = node.offDiagonal()[ interaction_idx ];

      if ( interaction._type != COMPRESSED_NODE ) {
        continue;
      }

      rank = blockRank( interaction._rowList.size(), nCols );

      totalRank += rank;
      maxRank = max( rank, maxRank );
    }

    storageEstimate += totalRank * nCols;

    for ( int interaction_idx = 0; interaction_idx < node.offDiagonal().size();
          interaction_idx++ )
    {
      const SupernodeInteraction &interaction
                                    = node.offDiagonal()[ interaction_idx ];

      if ( interaction._type == STANDARD_NODE
        || interaction._type == IMPLICIT_NODE ) {
        continue;
      }

      if ( interaction._type == EXTENDED_NODE ) {
        const Supernode       &ancestor = _factor[ interaction._nodeID ];

        nRows = ancestor.numExtendedColumns();

        TRACE_ASSERT( nRows > 0, "Empty extended node" );

        if ( interaction._compressed )
        {
          // Add space only for the columns we store
          totalStorage += nRows * interaction._compressedColumnList.size();
        }
        else
        {
          // Add space for the full column set
          totalStorage += nRows * nCols;
        }
      } else if ( interaction._type == COMPRESSED_NODE ) {
        rank = blockRank( interaction._rowList.size(), nCols );

#if 0
        totalStorage += rank * ( interaction._rowList.size() + nCols );
#endif
#if 0
        // FIXME: try using the maximum rank here for now
        totalStorage += maxRank * ( interaction._rowList.size() + nCols );
#endif
        // FIXME: try using the total rank here for now
        totalStorage += totalRank * ( interaction._rowList.size() + nCols );

        storageEstimate += totalRank * interaction._rowList.size();
      } else {
        TRACE_ASSERT( NULL, "Invalid interaction type" );
      }
    }
  }

  // Reverse any changes we've made
  clearSlackVariableEstimates();

  printf( "Estimated extra storage is actually %f MB\n",
          (Real)storageEstimate * 8.0 / 1024.0 / 1024.0 );
  //abort();

  return totalStorage;
}

//////////////////////////////////////////////////////////////////////
// Sets size estimates on all extended nodes
//////////////////////////////////////////////////////////////////////
void FactorManager::estimateSlackVariableSizes()
{
  int                        rank;
  int                        nCols;

  // For each node
  //    1) Look at each low rank block
  //    2) Estimate its maximum rank
  //    3) Based on this rank, figure out how much additional
  //       fill-in will be generated amongst extended nodes.
  //
  for ( int node_idx = 0; node_idx < _factor.size(); node_idx++ )
  {
    Supernode               &node = _factor[ node_idx ];

    nCols = node.numColumns();

    if ( node.type() == EXTENDED_NODE )
    {
      continue;
    }

    for ( int interaction_idx = 0; interaction_idx < node.offDiagonal().size();
          interaction_idx++ )
    {
      const SupernodeInteraction &interaction
                            = node.offDiagonal()[ interaction_idx ];

      if ( interaction._active )
      {
        continue;
      }

      const SupernodeInteraction &extendedInteraction
                      = node.offDiagonal()[ interaction._extendedInteraction ];

      rank = blockRank( interaction._rowList.size(), nCols );

      // Set a size estimate on the slack variable space for this interaction
      _factor[ extendedInteraction._nodeID ].numExtendedColumns() = rank;
    }

    // Next, handle the node's diagonal
    if ( !node.inPlaceDiagonal() ) {
      for ( int block_idx = 0; block_idx < node.diagonalLowRankBlocks().size();
            block_idx++ )
      {
        // The interaction corresponding to this block
        int interaction_idx = node.lowRankDiagonalInteraction( block_idx );

        const SupernodeInteraction &interaction
                              = node.offDiagonal()[ interaction_idx ];

        TRACE_ASSERT( interaction._type == EXTENDED_NODE, "Invalid interaction" );
        TRACE_ASSERT( _factor[ interaction._nodeID ].numExtendedColumns() == 0,
                      "Node already initialized" );

        const DenseBlock &block = node.diagonalLowRankBlocks()[ block_idx ];

        rank = diagonalBlockRank( block.numRows(), block.numColumns() );

        // Set a size estimate on the slack variable space for this interaction
        _factor[ interaction._nodeID ].numExtendedColumns() = rank;
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void FactorManager::clearSlackVariableEstimates()
{
  for ( int node_idx = 0; node_idx < _factor.size(); node_idx++ )
  {
    _factor[ node_idx ].numExtendedColumns() = 0;
  }
}

//////////////////////////////////////////////////////////////////////
// Clears all workspaces, preparing them for numerical factorization
//////////////////////////////////////////////////////////////////////
void FactorManager::clearWorkspace()
{
  memset( (void *)_multWorkspace, 0, _multWorkspaceSz * sizeof( Real ) );
  memset( (void *)_updateWorkspace, 0, _updateWorkspaceSz * sizeof( Real ) );
  
  memset( (void *)_factorData, 0, _factorDataSz * sizeof( Real ) );

  for ( int i = 0; i < _scatteredMap.size(); i++ )
  {
    _scatteredMap[ i ] = EMPTY;
  }

  for ( int i = 0; i < _lowRankScatteredMap.size(); i++ )
  {
    _lowRankScatteredMap[ i ] = EMPTY;
  }

  _relativeMap.clear();
  _lowRankRelativeMap.clear();

  for ( int i = 0; i < _factor.size(); i++ )
  {
    _nextInteraction[ i ] = EMPTY;
    _nextInteractionCache[ i ] = EMPTY;
    _head[ i ] = EMPTY;
    _next[ i ] = EMPTY;
  }

  _slackDataOffset = 0;
  _availableSlackData = _slackDataSz;

  if ( _slackData == NULL && _slackDataSz > 0 )
  {
    _slackData = new Real[ _slackDataSz ];
  }

  memset( (void *)_slackData, 0, _slackDataSz * sizeof( Real ) );
}

//////////////////////////////////////////////////////////////////////
// Gets ranks for each low rank block in the given node
//////////////////////////////////////////////////////////////////////
void FactorManager::setBlockRanks( const Supernode &node )
{
  int                        nCols = node.numColumns();
  int                        nRows;
  int                        interaction_idx;

  _blockRanks.clear();

  for ( int block_idx = 0; block_idx < _lowRankBlocks.size(); block_idx++ )
  {
    interaction_idx = _lowRankBlocks[ block_idx ];

    nRows = node.offDiagonal()[ interaction_idx ]._rowList.size();

    _blockRanks.push_back( fixedBlockRank( nRows, nCols ) );
  }
}

//////////////////////////////////////////////////////////////////////
// Applies diagonal updates to the given node using all current
// interior block descendents (assumes that findInteriorBlockDescendents
// has been called)
//////////////////////////////////////////////////////////////////////
void FactorManager::applyInteriorBlockDiagonalUpdate(
                                  const SPARSE_MATRIX::SparseColumnMatrix &A,
                                  Supernode &node )
{
  for ( int i = 0; i < _currentBlockDescendents.size(); i++ ) {
    applyInteriorBlockDiagonalUpdate( A, node, _currentBlockDescendents[ i ] );
  }
}

//////////////////////////////////////////////////////////////////////
// Applies diagonal updates to the given node using an interior block
//////////////////////////////////////////////////////////////////////
void FactorManager::applyInteriorBlockDiagonalUpdate(
                                  const SPARSE_MATRIX::SparseColumnMatrix &A,
                                  Supernode &node, int block_idx )
{
  const InteriorBlock       &block = _interiorBlocks[ block_idx ];

  int                        node_interaction = 0;
  size_t                     totalSize = 0;

  int                        nodeColumns = node.numColumns();

  int                        numNodes = range_size( block._nodeRange );

  _inverseRowMap.clear();
  _inverseRowMap.resize( node.numColumns(), EMPTY );

  node_interaction = findInteriorBlockAncestorInteraction( block, node );

  TRACE_ASSERT( node_interaction < block._interactions.size(),
                "Block does not interact with given node" );

  const InteriorBlock::Interaction  &blockInteraction
                                      = block._interactions[ node_interaction ];

  // Build a workspace
  totalSize = getInteriorBlockInteractionSize( block, blockInteraction );

  // Build a workspace big enough to expand the contents of this
  // block interaction
  PairArray                  dataSizes( 2 );

  dataSizes[ 0 ] = IndexPair( totalSize, 1 );
  dataSizes[ 1 ] = IndexPair( totalSize, 1 );

  RealWorkspace              workspace( _realWorkspaceManager, dataSizes );

  // Pointers to each workspace location for each node
  vector<Real *>             workspacePointers( numNodes, NULL );
  IntArray                   interactionMap( numNodes, EMPTY );

  MATRIX::clear( workspace.workspaceData( 0 ), totalSize, 1 );

  // Initialize workspaces for each block that we're going to deal
  // with when building this interaction.
  initializeInteriorBlockInteractionWorkspace( block, blockInteraction,
                                               totalSize,
                                               workspace.workspaceData( 0 ),
                                               workspacePointers,
                                               interactionMap );

  // Process contributions from each sub-interaction in this block
  // interaction
  for ( int i = 0; i < blockInteraction.size(); i++ ) {
    int                      node_idx = blockInteraction[ i ].first;
    int                      interaction_idx = blockInteraction[ i ].second;

    const Supernode             &blockNode = _factor[ node_idx ];
    const SupernodeInteraction  &interaction
                                  = blockNode.offDiagonal()[ interaction_idx ];

    int                          nCols = blockNode.numColumns();

    Real                        *data;

    data = workspacePointers[ node_idx - block._nodeRange.first ];

    // Get the parts of the sparse matrix that we need
    copyInteractionSparseMatrix( A, node, blockNode, interaction, data );

    // Solve for the contents in this block
    blockNode.diagonalSolve( data, interaction._rowList.size(), _slackData,
                             // Transposed right side solve
                             true, false );

    // Process fill-in from this node through the rest of the interaction
    // block row
    pushInteriorBlockInteractionFill( block, blockInteraction, node,
                                      blockNode, interaction,
                                      workspacePointers, interactionMap,
                                      workspace.workspaceData( 1 ) );

    // Finally, subtract this contribution from node's diagonal
    node.applyDiagonalUpdate( data, workspace.workspaceData( 1 ),
                              nCols, interaction._rowList );
  }
}

//////////////////////////////////////////////////////////////////////
// Same as the function above, but applies this for a specific block
// in this node's diagonal.  Also, places the result in a designated
// workspace, rather than the node's own data.
//////////////////////////////////////////////////////////////////////
void FactorManager::constructInteriorBlockDiagonalUpdate(
                              const SPARSE_MATRIX::SparseColumnMatrix &A,
                              Supernode &node,
                              const DenseBlock &factorBlock,
                              // Interior block index
                              int block_idx,
                              Real *schurComplement )
{
  const InteriorBlock       &block = _interiorBlocks[ block_idx ];

  int                        node_interaction = 0;
  size_t                     totalSize = 0;

  int                        nodeColumns = node.numColumns();

  int                        numNodes = range_size( block._nodeRange );

  RangeArray                 nodeRowRanges;

  _inverseRowMap.clear();
  _inverseRowMap.resize( node.numColumns(), EMPTY );

  // Find the interaction between 'block' and 'node'
  node_interaction = findInteriorBlockAncestorInteraction( block, node );

  TRACE_ASSERT( node_interaction < block._interactions.size(),
                "Block does not interact with given node" );

  const InteriorBlock::Interaction  &blockInteraction
                                      = block._interactions[ node_interaction ];

  // Get the size of the workspace we need to build.  Also, figure
  // out the row range over which each sub-interaction intersects
  // with the row range of factorBlock.
  totalSize = getInteriorBlockInteractionSize( block, blockInteraction,
                                               factorBlock._rowRange,
                                               &nodeRowRanges );

  // Build a workspace big enough to expand the contents of this
  // block interaction
  PairArray                  dataSizes( 2 );

  dataSizes[ 0 ] = IndexPair( totalSize, 1 );
  dataSizes[ 1 ] = IndexPair( totalSize, 1 );

  RealWorkspace              workspace( _realWorkspaceManager, dataSizes );

  // Pointers to each workspace location for each node
  vector<Real *>             workspacePointers( numNodes, NULL );
  IntArray                   interactionMap( numNodes, EMPTY );

  MATRIX::clear( workspace.workspaceData( 0 ), totalSize, 1 );

  // Initialize workspaces for each block that we're going to deal
  // with when building this interaction.
  initializeInteriorBlockInteractionWorkspace( block, blockInteraction,
                                               totalSize,
                                               workspace.workspaceData( 0 ),
                                               workspacePointers,
                                               interactionMap,
                                               &nodeRowRanges );

  // Process contributions from each sub-interaction in this block
  // interaction
  for ( int i = 0; i < blockInteraction.size(); i++ ) {
    int                      node_idx = blockInteraction[ i ].first;
    int                      interaction_idx = blockInteraction[ i ].second;

    const Supernode             &blockNode = _factor[ node_idx ];
    const SupernodeInteraction  &interaction
                                  = blockNode.offDiagonal()[ interaction_idx ];

    int                          nCols = blockNode.numColumns();

    IndexRange                   rowRange;

    Real                        *data;

    rowRange = nodeRowRanges[ node_idx - block._nodeRange.first ];

    // Skip if this interaction's row range doesn't intersect the row
    // range for this factor block
    if ( rowRange.first < 0 || rowRange.second < rowRange.first ) {
      continue;
    }

    data = workspacePointers[ block.localNodeIdx(node_idx) ];

    // Get the parts of the sparse matrix that we need
    copyInteractionSparseMatrix( A, node, blockNode, interaction, data,
                                 factorBlock._rowRange, rowRange );

    // Solve for the contents in this block
    blockNode.diagonalSolve( data, range_size( rowRange ), _slackData,
                             // Transposed right side solve
                             true, false );

    // Process fill-in from this node through the rest of the interaction
    // block row
    pushInteriorBlockInteractionFill( block, blockInteraction, node,
                                      blockNode, interaction,
                                      workspacePointers, interactionMap,
                                      workspace.workspaceData( 1 ),
                                      factorBlock._rowRange,
                                      &nodeRowRanges );

    // Finally, subtract this contribution from node's diagonal
    node.applyDiagonalUpdate( factorBlock,
                              data, workspace.workspaceData( 1 ), nCols,
                              interaction._rowList,
                              schurComplement,
                              rowRange );
  }
}

//////////////////////////////////////////////////////////////////////
// Finds the ancestor interaction between this interior block and
// the given node
//////////////////////////////////////////////////////////////////////
int FactorManager::findInteriorBlockAncestorInteraction(
                                                  const InteriorBlock &block,
                                                  const Supernode &node )
{
  int                        node_interaction = 0;

  // Find the interaction between this block and the given node
  for ( ; node_interaction < block._interactions.size(); node_interaction++ ) {
    int                      node_idx;
    int                      interaction_idx;

#if 0
    node_idx = block._interactions[ 0 ][ node_interaction ].first;
    interaction_idx = block._interactions[ 0 ][ node_interaction ].second;
#endif
    node_idx = block._interactions[ node_interaction ][ 0 ].first;
    interaction_idx = block._interactions[ node_interaction ][ 0 ].second;

    const Supernode         &blockNode = _factor[ node_idx ];

    if ( blockNode.offDiagonal()[ interaction_idx ]._nodeID == node.nodeID() ) {
      break;
    }
  }

  return node_interaction;
}

//////////////////////////////////////////////////////////////////////
// The number of non-zeros stored in an interior block interaction,
// assuming it is explicitly constructed
//////////////////////////////////////////////////////////////////////
size_t FactorManager::getInteriorBlockInteractionSize(
                            const InteriorBlock &block,
                            const InteriorBlock::Interaction &blockInteraction,
                            const IndexRange &rowRange,
                            RangeArray *nodeRowRanges )
{
  size_t                     totalSize = 0;

  if ( nodeRowRanges ) {
    nodeRowRanges->resize( block.numNodes() );
  }

  for ( int i = 0; i < blockInteraction.size(); i++ ) {
    int                      node_idx = blockInteraction[ i ].first;
    int                      interaction_idx = blockInteraction[ i ].second;

    const Supernode             &blockNode = _factor[ node_idx ];
    const SupernodeInteraction  &interaction
                                  = blockNode.offDiagonal()[ interaction_idx ];

    if ( rowRange.first < 0 ) {
      // Ignore the row range
      totalSize += blockNode.numColumns() * interaction._rowList.size();
    } else {
      int                    rowCount = 0;
      int                    startRow = 0;

      // Count rows inside the range
      for ( int row_idx = 0; row_idx < interaction._rowList.size(); row_idx++ ) {
        if ( interaction._rowList[ row_idx ] < rowRange.first ) {
          startRow += 1;
          continue;
        } else if ( interaction._rowList[ row_idx ] > rowRange.second ) {
          // The row list should be sorted, so this is safe
          break;
        } else {
          rowCount += 1;
        }
      }

      totalSize += blockNode.numColumns() * rowCount;

      if ( nodeRowRanges ) {
        nodeRowRanges->at( block.localNodeIdx(node_idx) )
          = IndexRange( startRow, startRow + rowCount - 1 );
      }
    }
  }

  return totalSize;
}

//////////////////////////////////////////////////////////////////////
// Initializes a workspace for explicitly expanding the contents
// of an interior block interaction.  Optionally provide local
// node row ranges if we are restricting ourselves to a particular
// block row of this interaction.
//////////////////////////////////////////////////////////////////////
void FactorManager::initializeInteriorBlockInteractionWorkspace(
                      const InteriorBlock &block,
                      const InteriorBlock::Interaction &blockInteraction,
                      size_t interactionSize,
                      Real *workspaceBase,
                      vector<Real *> &workspacePointers,
                      IntArray &interactionMap,
                      RangeArray *nodeRanges )
{
  // Pointers to each workspace location for each node
  workspacePointers.resize( block.numNodes() );
  interactionMap.resize( block.numNodes() );

  MATRIX::clear( workspaceBase, interactionSize, 1 );

  TRACE_ASSERT( nodeRanges == NULL
                  || nodeRanges->size() == block.numNodes() );

  for ( int i = 0; i < blockInteraction.size(); i++ ) {
    int                      node_idx = blockInteraction[ i ].first;
    int                      interaction_idx = blockInteraction[ i ].second;
    int                      localNodeIdx = block.localNodeIdx(node_idx);

    const Supernode             &blockNode = _factor[ node_idx ];
    const SupernodeInteraction  &interaction
                                = blockNode.offDiagonal()[ interaction_idx ];

    workspacePointers[ localNodeIdx ] = workspaceBase;
    interactionMap[ localNodeIdx ] = i;

    if ( nodeRanges != NULL ) {
      workspaceBase
        += blockNode.numColumns() * range_size(nodeRanges->at(localNodeIdx));
      interactionSize
        -= blockNode.numColumns() * range_size(nodeRanges->at(localNodeIdx));
    } else {
      workspaceBase += blockNode.numColumns() * interaction._rowList.size();
      interactionSize -= blockNode.numColumns() * interaction._rowList.size();
    }
  }

  TRACE_ASSERT( interactionSize >= 0,
                "Not enough space in interaction workspace" );
}

//////////////////////////////////////////////////////////////////////
// Copies sparse matrix entries from the given 'interaction' of 'blockNode'
// (which interacts with 'node') over the given row range in to the
// provided workspace.
//
// fullRowRange is the range of rows we wish to consider from the
// sparse matrix (relative to the starting column of node).
//
// localRowRange is the range of rows from 'interaction' that we
// need to worry about.
//////////////////////////////////////////////////////////////////////
void FactorManager::copyInteractionSparseMatrix(
                              const SPARSE_MATRIX::SparseColumnMatrix &A,
                              const Supernode &node,
                              const Supernode &blockNode,
                              const SupernodeInteraction &interaction,
                              Real *workspace,
                              IndexRange fullRowRange,
                              IndexRange localRowRange )
{
  // Error checking on inputs
  if ( fullRowRange.first < 0 || fullRowRange.second < fullRowRange.first
        || fullRowRange.second >= node.numColumns() ) {
    fullRowRange.first = 0;
    fullRowRange.second = node.numColumns() - 1;
  }

  if ( localRowRange.first < 0 || localRowRange.second < localRowRange.first
        || localRowRange.second >= interaction._rowList.size() ) {
    localRowRange.first = 0;
    localRowRange.second = interaction._rowList.size() - 1;
  }

  // Initialize an inverse row list
  invertPartialIntArray( interaction._rowList, _inverseRowMap, localRowRange );
  
  // Fill in sparse matrix entries in the workspace
  for ( int col_idx = blockNode.startColumn();
        col_idx <= blockNode.endColumn(); col_idx++ )
  {
    for ( int row_ptr = A._p[ col_idx ]; row_ptr < A._p[ col_idx + 1 ];
          row_ptr++ ) 
    {
      int                  row_idx = A._i[ row_ptr ] - node.startColumn();

      if ( row_idx < fullRowRange.first || row_idx > fullRowRange.second ) {
        continue;
      }

      row_idx = _inverseRowMap[ row_idx ];

      TRACE_ASSERT( row_idx >= 0 && row_idx < range_size(localRowRange) );

      MATRIX::access( workspace,
                      range_size(localRowRange), blockNode.numColumns(),
                      row_idx, col_idx - blockNode.startColumn() )
                            += A._x[ row_ptr ];
    }
  }

  clearPartialInverseArray( interaction._rowList, _inverseRowMap, EMPTY,
                            localRowRange );
}

//////////////////////////////////////////////////////////////////////
// Pushes forward fill-in when explicitly constructing the contents of
// an interior block interaction.
//
// Specifically, pushes fill from 'blockNode' to all ancestors.
// The data for each node is assumed to be stored in the provided
// workspace pointers.  'node' is the node that this particular
// block interaction is interacting with
//
// Optionally, fullRowRange specifies the block row of this interior
// block interaction that we are trying to form.
//
// If fullRowRange is specified then you should also specify a list of local
// row ranges for each sub-interaction in the interior block interaction.
//////////////////////////////////////////////////////////////////////
void FactorManager::pushInteriorBlockInteractionFill(
                              const InteriorBlock &block,
                              const InteriorBlock::Interaction &blockInteraction,
                              const Supernode &node,
                              const Supernode &blockNode,
                              const SupernodeInteraction &nodeInteraction,
                              const vector<Real *> &workspaces,
                              const IntArray &interactionMap,
                              Real *multWorkspace,
                              IndexRange fullRowRange,
                              RangeArray *localRowRanges )
{
  const Real                *data;
  IndexRange                 localRowRange;
  int                        localNodeID;
  int                        numNodes = range_size( block._nodeRange );

  localNodeID = block.localNodeIdx(blockNode.nodeID());

  // Error checking on inputs
  if ( fullRowRange.first < 0 || fullRowRange.second < fullRowRange.first
        || fullRowRange.second >= node.numColumns() ) {
    fullRowRange.first = 0;
    fullRowRange.second = node.numColumns() - 1;
  }

  // The workspace data for blockNode
  data = workspaces[ localNodeID ];

  if ( localRowRanges == NULL ) {
    localRowRange.first = 0;
    localRowRange.second = nodeInteraction._rowList.size() - 1;
  } else {
    localRowRange = localRowRanges->at( localNodeID );
  }

  // Push forward fill from each relevant interaction of blockNode
  for ( int fill_int_idx = 0; fill_int_idx < blockNode.offDiagonal().size();
        fill_int_idx++ )
  {
    int                          fill_node_id;
    const SupernodeInteraction  &fillInteraction
                                  = blockNode.offDiagonal()[ fill_int_idx ];

    Real                        *fillData;
    const Real                  *baseBlock;

    int                          recipient_id;
    int                          recipient_interaction_idx;

    IndexRange                   fillRowRange;

    if ( fillInteraction._type != STANDARD_NODE || !fillInteraction._active ) {
      continue;
    }

    fill_node_id = fillInteraction._nodeID - block._nodeRange.first;

    TRACE_ASSERT( fillInteraction._nodeID > blockNode.nodeID()
                    && fill_node_id < numNodes
                    && workspaces[ fill_node_id ] != NULL );

    fillData = workspaces[ fill_node_id ];

    baseBlock = blockNode.interactionMatrix( fill_int_idx );

    // Multiply the current block data with base block, and subtract
    // from fillData
    MATRIX::gemm( data, baseBlock, multWorkspace,
                  range_size( localRowRange ), blockNode.numColumns(),
                  fillInteraction._rowList.size(), blockNode.numColumns(),
                  // Transpose the second argument
                  false, true );

    // Find the node and interaction receiving fill-in here
    recipient_id = blockInteraction[ interactionMap[ fill_node_id ] ].first;
    recipient_interaction_idx
                 = blockInteraction[ interactionMap[ fill_node_id ] ].second;

    const Supernode &recipient = _factor[ recipient_id ];
    const IntArray  &fillRows
            = recipient.offDiagonal()[ recipient_interaction_idx ]._rowList;
    const IntArray  &fillColumns = fillInteraction._rowList;

    TRACE_ASSERT( recipient.nodeID() == fillInteraction._nodeID );
    TRACE_ASSERT( recipient.offDiagonal()[ recipient_interaction_idx ]._nodeID
                    == nodeInteraction._nodeID );

    if ( localRowRanges == NULL ) {
      fillRowRange.first = 0;
      fillRowRange.second = fillRows.size() - 1;
    } else {
      fillRowRange = localRowRanges->at( fill_node_id );
    }

    // We're filling in, so there must be rows here
    TRACE_ASSERT( range_size(fillRowRange) >= range_size(localRowRange) );

    // Build an inverse row map for the rows in the block we are
    // filling in
    invertPartialIntArray( fillRows, _inverseRowMap, fillRowRange );

    for ( int row_idx = localRowRange.first; row_idx <= localRowRange.second;
          row_idx++ )
    {
      int fillRow = _inverseRowMap[ nodeInteraction._rowList[ row_idx ] ];

      TRACE_ASSERT( fillRow >= 0 && fillRow < range_size(fillRowRange) );

      for ( int col_idx = 0; col_idx < fillColumns.size(); col_idx++ ) {
        int fillCol = fillColumns[ col_idx ];

        MATRIX::access( fillData,
                        range_size( fillRowRange ), recipient.numColumns(),
                        fillRow, fillCol )
          -= MATRIX::access( multWorkspace,
                             range_size( localRowRange ), fillColumns.size(),
                             row_idx - localRowRange.first, col_idx );
      }
    }

    clearPartialInverseArray( fillRows, _inverseRowMap, EMPTY, fillRowRange );
  }
}

//////////////////////////////////////////////////////////////////////
// Constructs a diagonal update matrix from 'descendent' to 'node' over
// the given sub-block of 'node's diagonal and subtracts it from the
// given Schur complement.
//
// This will be used for explicit diagonal block formation/compression
//////////////////////////////////////////////////////////////////////
void FactorManager::constructNodeDiagonalUpdate( const Supernode &node,
                                                 const Supernode &descendent,
                                                 const DenseBlock &factorBlock,
                                                 Real *schurComplement )
{
  int                        ancestor_idx;
  IndexRange                 rowRange;
  int                        nRows;
  
  ancestor_idx = _nextInteractionCache[ descendent.nodeID() ];

  // The interaction between descendent and node
  const SupernodeInteraction &interaction
                                    = descendent.offDiagonal()[ ancestor_idx ];

  TRACE_ASSERT( interaction._nodeID == node.nodeID(),
                "Incorrect ancestor interaction" );

  // Figure out the row subset of this interaction that we need to use
  findEntryRangeIntersection( interaction._rowList,
                              factorBlock._rowRange, rowRange );

  // If no range is found, skip this node
  if ( rowRange.first < 0 || rowRange.second < rowRange.first ) {
    return;
  }

  nRows = range_size( rowRange );

  if ( interaction._type == COMPRESSED_NODE ) {
    const Real              *V = interaction._compressedV; 
    const Real              *U = interaction._compressedU; 
    int                      rank = interaction._numExtendedRows;

    // Make a workspace for storing matrix products
    //
    // We need a workspace of size (rank * rank) to form U' * U,
    // then a workspace of size rank * nRows to form (U' * U) * V',
    // then a workspace of size nRows * nRows to form V * (U' * U) * V'
    RealWorkspace            workspace( _realWorkspaceManager,
                                        rank * (rank + nRows) + nRows * nRows );

    node.compressedDiagonalUpdate( factorBlock, V, U, rank,
                                   workspace.workspaceData( 0 ),
                                   descendent.numColumns(),
                                   interaction._rowList,
                                   schurComplement, rowRange );
  } else {
    const Real              *interactionData;

    interactionData = descendent.interactionMatrix( ancestor_idx );
    interactionData += rowRange.first * descendent.numColumns();

    // Get a workspace capable of storing the matrix product
    RealWorkspace              workspace( _realWorkspaceManager, nRows * nRows );

    node.applyDiagonalUpdate( factorBlock, interactionData,
                              workspace.workspaceData( 0 ),
                              descendent.numColumns(),
                              interaction._rowList, schurComplement,
                              rowRange );
  }
}

#if 0
//////////////////////////////////////////////////////////////////////
// Initializes workspaces for interaction multiplication.
// This basically just means stacking a bunch of random matrices
// on top of each other in a work space.
//////////////////////////////////////////////////////////////////////
void FactorManager::initMultiplyWorkspace( const Supernode &node,
                                           int fixedRank )
{
  int                        nCols = node.numColumns();
  Real                      *baseData = _decompMultTransWorkspace;
  int                        rank;

  for ( int block_idx = 0; block_idx < _lowRankBlocks.size(); block_idx++ )
  {
    if ( !_lowRankBlockActive[ block_idx ] )
    {
      continue;
    }

    rank = ( fixedRank > 0 ) ? fixedRank : _blockRanks[ block_idx ];

    MATRIX::copy( baseData, _testMatrix.data(), nCols, rank );

    baseData += nCols * rank;
  }
}
#endif

//////////////////////////////////////////////////////////////////////
// Multiplies data in _decompMultTransWorkspace on the left
// with each low rank interaction matrix for the current node.
// Puts the result in _decompWorkspace.
//
// (ie. half a step of power iteration)
//////////////////////////////////////////////////////////////////////
void FactorManager::interactionMultiply(
                              const SPARSE_MATRIX::SparseColumnMatrix &A,
                              Supernode &node,
                              bool useFactor,
                              Real *inputData,
                              Real *outputData,
                              int fixedRank )
{
  int                        nCols = node.numColumns();
  int                        interaction_idx;
  int                        rank;

  MATRIX                     solveWorkspace;
  const Real                *multInput;

  IntArray                   nodeList;

  // Workspace stuff
  PairArray                  dataSizes( 3 );

  if ( !inputData ) {
    inputData = _decompMultTransWorkspace;
  }

  if ( !outputData ) {
    outputData = _decompWorkspace;
  }

  // Make space for storing node lists if we're doing interior block multiplies
  if ( _useInteriorBlocks ) {
    nodeList.reserve( _factor.size() );
  }

  for ( int block_idx = 0; block_idx < _lowRankBlocks.size(); block_idx++ )
  {
    if ( !_lowRankBlockActive[ block_idx ] )
    {
      continue;
    }

    interaction_idx = _lowRankBlocks[ block_idx ];

    rank = ( fixedRank > 0 ) ? fixedRank : _blockRanks[ block_idx ];

    const SupernodeInteraction &interaction
                          = node.offDiagonal()[ interaction_idx ];

    // If we are decomposing the factored block, apply the inverse
    // transpose of the factor here
    //
    // NOTE: This changes the input
    TIMING_START( "interactionMultiply: diagonalSolve" );
    if ( useFactor ) {
      // We need to make a copy of the input here, since applying the
      // solver directly will change it.
      solveWorkspace.resizeAndWipe( nCols, rank );

      MATRIX::copy( solveWorkspace.data(), inputData, nCols, rank );

      node.diagonalSolve( solveWorkspace.data(), rank, _slackData, true );

      multInput = solveWorkspace.data();
    } else {
      multInput = inputData;
    }
    TIMING_STOP( "interactionMultiply: diagonalSolve" );

    TIMING_START( "interactionMultiply: baseSystemMultiply" );
    node.baseSystemMultiply( A, _factor, interaction_idx,
                             multInput, rank, outputData );
    TIMING_STOP( "interactionMultiply: baseSystemMultiply" );

    PairArray &descendents = _lowRankDescendents[ block_idx ];

    // Interaction multiplies with each descendent
    interactionMultiply_nodeDescendents( node, interaction, descendents,
                                         dataSizes, rank,
                                         multInput, outputData );

    if ( _useInteriorBlocks ) {
      interactionMultiply_blockDescendents( A, node, interaction,
                                            rank, multInput, outputData,
                                            nodeList,
                                            false /* no transpose */ );
    }

    inputData += nCols * rank;
    outputData += interaction._rowList.size() * rank;
  }
}

//////////////////////////////////////////////////////////////////////
// Same as the above, but multiplies the entire off-diagonal block
// all at once
//////////////////////////////////////////////////////////////////////
void FactorManager::interactionMultiply_allRows(
                          const SPARSE_MATRIX::SparseColumnMatrix &A,
                          Supernode &node,
                          bool useFactor,
                          Real *inputData,
                          Real *outputData,
                          int fixedRank )
{
  int                        nCols = node.numColumns();
  int                        rank; 

  MATRIX                     solveWorkspace;
  const Real                *multInput;

  IntArray                   nodeList;

  // Workspace stuff
  PairArray                  dataSizes( 3 );

  if ( !inputData ) {
    inputData = _decompMultTransWorkspace;
  }

  if ( !outputData ) {
    outputData = _decompWorkspace;
  }

  // Make space for storing node lists if we're doing interior block
  // multiplies
  if ( _useInteriorBlocks ) {
    nodeList.reserve( _factor.size() );
  }

  TRACE_ASSERT( fixedRank > 0 || _blockRanks.size() > 0 );
  rank = ( fixedRank > 0 ) ? fixedRank : _blockRanks[ 0 ];

  TIMING_START( "interactionMultiply_allRows: diagonalSolve" );
  if ( useFactor ) {
    // We need to make a copy of the input here, since applying
    // the solver directly will change it
    solveWorkspace.resizeAndWipe( nCols, rank );

    MATRIX::copy( solveWorkspace.data(), inputData, nCols, rank );

    node.diagonalSolve( solveWorkspace.data(), rank, _slackData, true );

    multInput = solveWorkspace.data();
  } else {
    multInput = inputData;
  }
  TIMING_STOP( "interactionMultiply_allRows: diagonalSolve" );

  TIMING_START( "interactionMultiply_allRows: baseSystemMultiply" );
  node.baseSystemMultiply( A, _factor, multInput, rank, outputData );
  TIMING_STOP( "interactionMultiply_allRows: baseSystemMultiply" );

  interactionMultiply_nodeDescendents( node, dataSizes, rank,
                                       multInput, outputData,
                                       false /* No transpose */ );

  if ( _useInteriorBlocks ) {
    interactionMultiply_blockDescendents( A, node, rank,
                                          multInput, outputData,
                                          nodeList,
                                          false /* No transpose */ );
  }
}

//////////////////////////////////////////////////////////////////////
// Processes normal descendent multiplications when forming a low-rank
// interaction
//////////////////////////////////////////////////////////////////////
void FactorManager::interactionMultiply_nodeDescendents(
                                    const Supernode &node,
                                    const SupernodeInteraction &interaction,
                                    const PairArray &descendents,
                                    PairArray &dataSizes,
                                    int rank,
                                    const Real *multInput,
                                    Real *outputData )
{
  int                        descendent_idx;
  int                        ancestor_idx;
  int                        desc_interaction_idx;

  // Interaction multiplies with each descendent
  for ( int i = 0; i < descendents.size(); i++ )
  {
    descendent_idx = descendents[ i ].first;

    if ( isInInteriorBlock( descendent_idx ) ) {
      continue;
    }

#if 0
    cout << "Applying node descendent " << descendent_idx << endl;
#endif

    desc_interaction_idx = descendents[ i ].second;
    ancestor_idx = _nextInteractionCache[ descendent_idx ];

    Supernode         &descendent = _factor[ descendent_idx ];

    int                descendentRows;


    // Construct workspaces
    //
    // Sub-matrix workspace
    descendentRows = descendent.offDiagonal()[ ancestor_idx ]._rowList.size();
    dataSizes[ 0 ] = IndexPair( descendentRows, rank );
    // Initial interaction multiply workspace
    dataSizes[ 1 ] = IndexPair( descendent.numColumns(), rank );
    // Final interaction multiply workspace
    descendentRows
      = descendent.offDiagonal()[ desc_interaction_idx ]._rowList.size();
    dataSizes[ 2 ] = IndexPair( descendentRows, rank );

    RealWorkspace          workspace( _realWorkspaceManager, dataSizes );

    Real                  *subMatrix = workspace.workspaceData( 0 );
    Real                  *multWorkspaceInitial = workspace.workspaceData( 1 );
    Real                  *multWorkspaceFinal = workspace.workspaceData( 2 );

    // Form a sub-matrix of our precomputed random matrix based on
    // the row set in the interaction between descendent and node
    TIMING_START( "interactionMultiply: BuildInteractionSubMatrix" );
    Supernode::BuildInteractionSubMatrix(
                                  descendent.offDiagonal()[ ancestor_idx ],
                                  multInput, subMatrix, rank );
    TIMING_STOP( "interactionMultiply: BuildInteractionSubMatrix" );

    TIMING_START( "interactionMultiply: interactionMultiply" );
    descendent.interactionMultiply( desc_interaction_idx, ancestor_idx,
                                    node, subMatrix, rank,
                                    multWorkspaceInitial,
                                    multWorkspaceFinal );
    TIMING_STOP( "interactionMultiply: interactionMultiply" );

    // Scatter back to the appropriate workspace
    TIMING_START( "interactionMultiply: ScatterLowRankUpdate" );
    Supernode::ScatterLowRankUpdate(
                          multWorkspaceFinal, outputData,
                          descendent.offDiagonal()[ desc_interaction_idx ],
                          interaction._rowList, rank );
    TIMING_STOP( "interactionMultiply: ScatterLowRankUpdate" );
  }
}

//////////////////////////////////////////////////////////////////////
// Overloaded version of the above function which is assumes that
// we are simultaneously decomposing the entire off-diagonal block
// of node
//////////////////////////////////////////////////////////////////////
void FactorManager::interactionMultiply_nodeDescendents( const Supernode &node,
                                                         PairArray &dataSizes,
                                                         int rank,
                                                         const Real *multInput,
                                                         Real *outputData,
                                                         bool transpose )
{
  int                        descendent_idx;
  int                        ancestor_idx;

  // Interaction multiplies with each descendent
  for ( int i = 0; i < _currentDescendents.size(); i++ )
  {
    descendent_idx = _currentDescendents[ i ];

    if ( isInInteriorBlock( descendent_idx ) ) {
      continue;
    }

    ancestor_idx = _nextInteractionCache[ descendent_idx ];

    Supernode               &descendent = _factor[ descendent_idx ];

    int                      descendentRows;

    if ( ancestor_idx == descendent.offDiagonal().size() - 1 ) {
      continue;
    }

    // Construct workspaces
    //
    // Sub-matrix workspaces
    descendentRows = descendent.offDiagonal()[ ancestor_idx ]._rowList.size();
    dataSizes[ 0 ] = IndexPair( descendentRows, rank );
    // Initial interaction multiply workspace
    dataSizes[ 1 ] = IndexPair( descendent.numColumns(), rank );
    // Final interaction multiply workspace
    descendentRows = descendent.countInteractionRows(
                        ancestor_idx + 1, descendent.offDiagonal().size() - 1 );
    dataSizes[ 2 ] = IndexPair( descendentRows, rank );

    RealWorkspace          workspace( _realWorkspaceManager, dataSizes );

    Real                  *subMatrix;
    Real                  *multWorkspaceInitial;
    Real                  *multWorkspaceFinal;

    if ( transpose ) {
      subMatrix = workspace.workspaceData( 2 );
      multWorkspaceInitial = workspace.workspaceData( 1 );
      multWorkspaceFinal = workspace.workspaceData( 0 );
    } else {
      subMatrix = workspace.workspaceData( 0 );
      multWorkspaceInitial = workspace.workspaceData( 1 );
      multWorkspaceFinal = workspace.workspaceData( 2 );
    }

    // Form a sub-matrix of our precomputed random matrix based on the row
    // set in the interaction between descendent and node
    if ( transpose ) {
      TIMING_START( "interactionTransMultiply: BuildMultiInteractionSubMatrix" );
      Supernode::BuildMultiInteractionSubMatrix(
                                    descendent, ancestor_idx, node,
                                    multInput, subMatrix, rank );
      TIMING_STOP( "interactionTransMultiply: BuildMultiInteractionSubMatrix" );
    } else {
      TIMING_START( "interactionMultiply: BuildInteractionSubMatrix" );
      Supernode::BuildInteractionSubMatrix(
                                    descendent.offDiagonal()[ ancestor_idx ],
                                    multInput, subMatrix, rank );
      TIMING_STOP( "interactionMultiply: BuildInteractionSubMatrix" );
    }

    // Next, multiply with all interactions after ancestor_idx in
    // descendent
    if ( transpose ) {
      TIMING_START( "interactionTransMultiply: interactionMultiply" );
      descendent.interactionTransMultiply_left( ancestor_idx,
                                                node, subMatrix, rank,
                                                multWorkspaceInitial,
                                                multWorkspaceFinal );
      TIMING_STOP( "interactionTransMultiply: interactionMultiply" );
    } else {
      TIMING_START( "interactionMultiply: interactionMultiply" );
      descendent.interactionMultiply( ancestor_idx,
                                      node, subMatrix, rank,
                                      multWorkspaceInitial,
                                      multWorkspaceFinal );
      TIMING_STOP( "interactionMultiply: interactionMultiply" );
    }

    // Scatter back to the appropriate workspace
    if ( transpose ) {
      // Scatter back to the desired row space
      Supernode::ScatterLowRankTransUpdate_row( ancestor_idx,
                                                descendent, node,
                                                multWorkspaceFinal,
                                                outputData, rank );
    } else {
      TIMING_START( "interactionMultiply: ScatterMultiLowRankUpdate" );
      Supernode::ScatterMultiLowRankUpdate( multWorkspaceFinal, outputData,
                                            descendent, ancestor_idx,
                                            node, rank );
      TIMING_STOP( "interactionMultiply: ScatterMultiLowRankUpdate" );
    }
  }
}

//////////////////////////////////////////////////////////////////////
// Processes interior block descendent multiplications when forming
// a low-rank interaction
//////////////////////////////////////////////////////////////////////
void FactorManager::interactionMultiply_blockDescendents(
                                    const SPARSE_MATRIX::SparseColumnMatrix &A,
                                    const Supernode &node,
                                    const SupernodeInteraction &interaction,
                                    int rank,
                                    const Real *multInput,
                                    Real *outputData,
                                    IntArray &nodeList,
                                    bool transpose )
{
  int                        descendent_idx;
  int                        base_block_interaction;
  int                        block_interaction;
  int                        nCols = node.numColumns();

  // If we are doing interior block multiplies, then an inverse of
  // this interaction's row list will come in handy
  _inverseRowMap.clear();
  _inverseRowMap.resize( _factor[ interaction._nodeID ].numColumns(),
                         EMPTY );

  for ( int row_idx = 0; row_idx < interaction._rowList.size(); row_idx++ )
  {
    _inverseRowMap[ interaction._rowList[ row_idx ] ] = row_idx;
  }

#if 0
  cout << "Interaction size " << interaction._rowList.size();
#endif

  // Interaction multiply with each interior block descendent
  for ( int i = 0; i < _currentBlockDescendents.size(); i++ )
  {
#ifdef DO_TIMING
    _timers[ INTERIOR_BLOCK_PREAMBLE ].tick();
#endif
    TIMING_START( "Interior block preamble" );

    descendent_idx = _currentBlockDescendents[ i ];

    const InteriorBlock   &block = _interiorBlocks[ descendent_idx ];

#if 0
    cout << "Applying block descendent " << descendent_idx;
    cout << " with size " << range_size( block._nodeRange ) << endl;
#endif

    int                    totalColumns;

#if 0
    printf( "Interaction multiply with block %d which has %d interactions\n",
            descendent_idx, (int)block._interactions.size() );
#endif

    // Determine the interaction index in this block which refers to node
    // (it must exist) and the one which applies to the current low rank
    // interaction being processed (doesn't necessarily have to exist)
    base_block_interaction = block._interactionMap[ node.nodeID() ];
    block_interaction = block._interactionMap[ interaction._nodeID ];

    TRACE_ASSERT( base_block_interaction >= 0,
                  "No interaction found in interior block" );

    TIMING_STOP( "Interior block preamble" );
#ifdef DO_TIMING
    _timers[ INTERIOR_BLOCK_PREAMBLE ].tock();
#endif

    if ( block_interaction < 0 ) {
      continue;
    }

#ifdef DO_TIMING
    _timers[ INTERIOR_BLOCK_NODELIST ].tick();
#endif
    TIMING_START( "Interior block node list" );

    nodeList.clear();

    // Build the list of nodes in this block which contribute to this
    // low rank interaction
    buildInteriorBlockNodeList( descendent_idx, base_block_interaction,
                                block_interaction, nodeList );

#if 0
    {
      nodeList.clear();
      for ( int node_idx = block._nodeRange.first;
            node_idx <= block._nodeRange.second; node_idx++ )
      {
        nodeList.push_back( node_idx );
      }
    }
    cout << SDUMP( nodeList.size() ) << endl;
#endif

    TIMING_STOP( "Interior block node list" );
#ifdef DO_TIMING
    _timers[ INTERIOR_BLOCK_NODELIST ].tock();
#endif

    if ( nodeList.size() == 0 )
    {
      continue;
    }

#ifdef DO_TIMING
    _timers[ INTERIOR_BLOCK_NODELIST ].tick();
#endif
    totalColumns = nodeColumns( nodeList );
#ifdef DO_TIMING
    _timers[ INTERIOR_BLOCK_NODELIST ].tock();
#endif

#ifdef DO_TIMING
    _timers[ INTERIOR_BLOCK_WORKSPACE ].tick();
#endif

    // Generate workspaces for interior block multiplication
#if 0
    RealWorkspace            workspace( _realWorkspaceManager,
                                        IndexPair( totalColumns, rank ) );
#endif
    RealWorkspace            workspace( _realWorkspaceManager,
                                        IndexPair( totalColumns, 2 * rank ) );

#ifdef DO_TIMING
      _timers[ INTERIOR_BLOCK_WORKSPACE ].tock();
#endif

    // Perform the actual multiplication, subtracting from the
    // output workspace
#if 0
    interiorBlockSchurMultiply( A, nodeList, descendent_idx,
                                base_block_interaction, block_interaction,
                                // Input matrix and its size
                                multInput,
                                transpose ?  interaction._rowList.size() : nCols,
                                rank,
                                // Map of matrix rows back to rows in the
                                // interaction and the number of rows in
                                // the interaction
                                _inverseRowMap, interaction._rowList.size(),
                                workspace.workspaceData( 0 ),
                                outputData,
                                // Dummy row/column ranges
                                IndexRange( -1, -1 ), IndexRange( -1, -1 ),
                                transpose );
#endif
    interiorBlockSchurMultiply( A, nodeList, descendent_idx,
                                base_block_interaction, block_interaction,
                                // Input matrix and its size
                                multInput,
                                transpose ?  interaction._rowList.size() : nCols,
                                rank,
                                // Map of matrix rows back to rows in the
                                // interaction and the number of rows in
                                // the interaction
                                _inverseRowMap, interaction._rowList.size(),
                                workspace.workspaceData( 0 ),
                                workspace.workspaceData( 0 )
                                  + totalColumns * rank,
                                outputData,
                                // Dummy row/column ranges
                                IndexRange( -1, -1 ), IndexRange( -1, -1 ),
                                transpose );
  }
}

//////////////////////////////////////////////////////////////////////
// Same as the above, but assumes we are forming an off-diagonal
// decomposition for node all at once
//////////////////////////////////////////////////////////////////////
void FactorManager::interactionMultiply_blockDescendents(
                                const SPARSE_MATRIX::SparseColumnMatrix &A,
                                const Supernode &node,
                                int rank,
                                const Real *multInput,
                                Real *outputData,
                                IntArray &nodeList,
                                bool transpose )
{
  int                        descendent_idx;
  int                        base_block_interaction;
  int                        startRow = node.endColumn() + 1;
  int                        nRows = A._nrow - startRow;
  int                        nRowsCompressed;
  int                        nCols = node.numColumns();
  int                        currentRow = 0;

  // If we are doing interior block multiplication, then an inverse
  // of the off-diagonal row set for node will come in handy
  _inverseRowMap.clear();
  _inverseRowMap.resize( nRows, EMPTY );

  for ( int interaction_idx = 0; interaction_idx < node.offDiagonal().size();
        interaction_idx++ )
  {
    const SupernodeInteraction  &interaction
                                      = node.offDiagonal()[ interaction_idx ];
    const Supernode             &ancestor = _factor[ interaction._nodeID ];

    int                          row;

    for ( int row_idx = 0; row_idx < interaction._rowList.size(); row_idx++ ) {
      row = ancestor.startColumn() + interaction._rowList[ row_idx ] - startRow;

      _inverseRowMap[ row ] = currentRow;
      currentRow += 1;
    }
  }
  nRowsCompressed = currentRow;

  // Interaction multiply with each interior block descendent
  for ( int i = 0; i < _currentBlockDescendents.size(); i++ )
  {
#ifdef DO_TIMING
    _timers[ INTERIOR_BLOCK_PREAMBLE ].tick();
#endif
    TIMING_START( "Interior block preamble" );

    descendent_idx = _currentBlockDescendents[ i ];

    const InteriorBlock   &block = _interiorBlocks[ descendent_idx ];

    int                    totalColumns;

    // Determine the interaction index in this block which refers to node
    // (it must exist) and the one which applies to the current low rank
    // interaction being processed (doesn't necessarily have to exist)
    base_block_interaction = block._interactionMap[ node.nodeID() ];

    TIMING_STOP( "Interior block preamble" );
#ifdef DO_TIMING
    _timers[ INTERIOR_BLOCK_PREAMBLE ].tock();
#endif

   // Skip blocks for which this is the last interaction
    if ( base_block_interaction >= block._interactions.size() - 1 ) {
      continue;
    }

    TRACE_ASSERT( base_block_interaction >= 0,
                  "No interaction found in interior block" );

#ifdef DO_TIMING
    _timers[ INTERIOR_BLOCK_NODELIST ].tick();
#endif
    TIMING_START( "Interior block node list" );

    nodeList.clear();

    // Build the list of nodes in this block which contribute
    // to the off-diagonal low rank interaction
    buildInteriorBlockNodeList( descendent_idx, base_block_interaction,
                                nodeList );

    TIMING_STOP( "Interior block node list" );
#ifdef DO_TIMING
    _timers[ INTERIOR_BLOCK_NODELIST ].tock();
#endif

    if ( nodeList.size() == 0 ) {
      continue;
    }

#ifdef DO_TIMING
    _timers[ INTERIOR_BLOCK_NODELIST ].tick();
#endif
    totalColumns = nodeColumns( nodeList );
#ifdef DO_TIMING
    _timers[ INTERIOR_BLOCK_NODELIST ].tock();
#endif

#ifdef DO_TIMING
    _timers[ INTERIOR_BLOCK_WORKSPACE ].tick();
#endif
    RealWorkspace            workspace( _realWorkspaceManager,
                                        IndexPair( totalColumns, 2 * rank ) );

#ifdef DO_TIMING
    _timers[ INTERIOR_BLOCK_WORKSPACE ].tock();
#endif
    interiorBlockSchurMultiply( A, nodeList,
                                descendent_idx, base_block_interaction,
                                // Input matrix and its size
                                multInput,
                                transpose ? nRowsCompressed : nCols, rank,
                                // Map of matrix rows back to rows in the
                                // off-diagonal matrix the number of rows 
                                // in this matrix
                                _inverseRowMap, nRowsCompressed,
                                workspace.workspaceData( 0 ),
                                workspace.workspaceData( 0 )
                                  + totalColumns * rank,
                                outputData,
                                node.endColumn() + 1,
                                transpose );
  }
}

#if 0
//////////////////////////////////////////////////////////////////////
// Version of interaction multiply that uses the interior block
// formulation
//
// Assumes that findInteriorBlockDescendents has been called
//////////////////////////////////////////////////////////////////////
void FactorManager::interiorBlockInteractionMultiply(
                              const SPARSE_MATRIX::SparseColumnMatrix &A,
                              Supernode &node,
                              const Real *inputData,
                              Real *outputData,
                              int fixedRank )
{
  int                        nCols = node.numColumns();
  int                        numRowsUncompressed;
  int                        interaction_idx;
  int                        descendent_idx;
  int                        rank;

  int                        base_block_interaction;
  int                        block_interaction;

#ifdef DO_TIMING
  _timers[ INTERIOR_BLOCK_PREAMBLE ].tick();
#endif

  IntArray                   nodeList;

  nodeList.reserve( _factor.size() );

  if ( !inputData )
  {
    inputData = _decompMultTransWorkspace;
  }

  if ( !outputData )
  {
    outputData = _decompWorkspace;
  }

#ifdef DO_TIMING
  _timers[ INTERIOR_BLOCK_PREAMBLE ].tock();
#endif

#if 0
  printf( "Node %d has %d interior block descendents\n",
          (int)node.nodeID(), (int)_currentBlockDescendents.size() );
#endif

  for ( int block_idx = 0; block_idx < _lowRankBlocks.size(); block_idx++ )
  {
    if ( !_lowRankBlockActive[ block_idx ] )
    {
      continue;
    }

#ifdef DO_TIMING
    _timers[ INTERIOR_BLOCK_PREAMBLE ].tick();
#endif

    interaction_idx = _lowRankBlocks[ block_idx ];

    rank = ( fixedRank > 0 ) ? fixedRank : _blockRanks[ block_idx ];

    const SupernodeInteraction &interaction
                          = node.offDiagonal()[ interaction_idx ];

    node.baseSystemMultiply( A, _factor, interaction_idx,
                             inputData, rank, outputData );

    // Determine the number of rows in the full uncompressed block
    numRowsUncompressed = _factor[ interaction._nodeID ].numColumns();

#ifdef DO_TIMING
    _timers[ INTERIOR_BLOCK_PREAMBLE ].tock();
#endif

    // Interaction multiply with each interior block descendent
    for ( int i = 0; i < _currentBlockDescendents.size(); i++ )
    {
#ifdef DO_TIMING
      _timers[ INTERIOR_BLOCK_PREAMBLE ].tick();
#endif

      descendent_idx = _currentBlockDescendents[ i ];

      const InteriorBlock   &block = _interiorBlocks[ descendent_idx ];

      int                    totalColumns;

#if 0
      printf( "Interaction multiply with block %d which has %d interactions\n",
              descendent_idx, (int)block._interactions.size() );
#endif

      // Determine the interaction index in this block which refers to node
      // (it must exist) and the one which applies to the current low rank
      // interaction being processed (doesn't necessarily have to exist)
      base_block_interaction = block._interactionMap[ node.nodeID() ];
      block_interaction = block._interactionMap[ interaction._nodeID ];

      TRACE_ASSERT( base_block_interaction >= 0,
                    "No interaction found in interior block" );

#ifdef DO_TIMING
      _timers[ INTERIOR_BLOCK_PREAMBLE ].tock();
#endif

      if ( block_interaction < 0 )
      {
        continue;
      }

#ifdef DO_TIMING
      _timers[ INTERIOR_BLOCK_NODELIST ].tick();
#endif

      nodeList.clear();

      // Build the list of nodes in this block which contribute to this
      // low rank interaction
      buildInteriorBlockNodeList( descendent_idx, base_block_interaction,
                                  block_interaction, nodeList );

#ifdef DO_TIMING
      _timers[ INTERIOR_BLOCK_NODELIST ].tock();
#endif

      if ( nodeList.size() == 0 )
      {
        continue;
      }

#ifdef DO_TIMING
      _timers[ INTERIOR_BLOCK_NODELIST ].tick();
#endif
      totalColumns = nodeColumns( nodeList );
#ifdef DO_TIMING
      _timers[ INTERIOR_BLOCK_NODELIST ].tock();
#endif

#ifdef DO_TIMING
      _timers[ INTERIOR_BLOCK_WORKSPACE ].tick();
#endif

      // Generate workspaces for interior block multiplication
      PairArray              workspaceSizes( 2 );

      workspaceSizes[ 0 ] = IndexPair( totalColumns, rank );
      workspaceSizes[ 1 ] = IndexPair( numRowsUncompressed, rank );

      RealWorkspace          workspace( _realWorkspaceManager, workspaceSizes );

#ifdef DO_TIMING
      _timers[ INTERIOR_BLOCK_WORKSPACE ].tock();
#endif

      // Perform the multiplication using these workspaces
      interiorBlockSchurMultiply( A, nodeList, descendent_idx,
                                  base_block_interaction, block_interaction,
                                  // Input matrix and its size
                                  inputData, nCols, rank,
                                  workspace.workspaceData( 0 ),
                                  workspace.workspaceData( 1 ) );

#if 0
      // FIXME
      {
        MATRIX tmp1( nCols, rank, inputData );
        MATRIX tmp2( numRowsUncompressed, rank, workspace.workspaceData( 1 ) );

        tmp1.write( "test_schurMult_input.matrix" );
        tmp2.write( "test_schurMult_output.matrix" );
      }
#endif

      // Subtract the non-zero rows from the workspace from the final
      // output matrix
#if 0
      printf( "Copying %d rows\n", (int)interaction._rowList.size() );
#endif
#ifdef DO_TIMING
      _timers[ INTERIOR_BLOCK_COPY ].tick();
#endif
      MATRIX::copyAddRows( workspace.workspaceData( 1 ), outputData, 
                           interaction._rowList, rank,
                           // Subtract
                           -1.0 );
#ifdef DO_TIMING
      _timers[ INTERIOR_BLOCK_COPY ].tock();
#endif

#if 0
      // FIXME
      {
        MATRIX tmp1( interaction._rowList.size(), rank );

        MATRIX::copyAddRows( workspace.workspaceData( 1 ), tmp1.data(),
                             interaction._rowList, rank, -1.0 );

        tmp1.write( "test_schurMult_copy.matrix" );
      }
#endif
    }

    inputData += nCols * rank;
    outputData += interaction._rowList.size() * rank;
  }
}
#endif

//////////////////////////////////////////////////////////////////////
// Similar to the above function, but does the multiplication with
// a sub-matrix of the main diagonal (indexed by block_idx)
//////////////////////////////////////////////////////////////////////
void FactorManager::diagonalBlockMultiply(
                              const SPARSE_MATRIX::SparseColumnMatrix &A,
                              Supernode &node,
                              bool useFactor,
                              int block_idx,
                              int rank,
                              const Real *inputData,
                              Real *outputData,
                              bool write )
{
  const DenseBlock &block = node.diagonalLowRankBlocks()[ block_idx ];

  int                        descendent_idx;
  int                        ancestor_idx;

  const Real                *multInput;
  MATRIX                     solveWorkspace;

  IntArray                   nodeList;

  // Workspace sizes
  PairArray                  dataSizes( 3 );

  if ( !inputData ) {
    inputData = _decompMultTransWorkspace;
  }

  if ( !outputData ) {
    outputData = _decompWorkspace;
  }

  // Make space for storing node lists if we're doing interior block multiplies
  if ( _useInteriorBlocks ) {
    nodeList.reserve( _factor.size() );
  }

  // If we are decomposing the factored block, apply the inverse
  // transpose to the factor here
  TIMING_START( "diagonalBlockMultiply: in place solve" );
  if ( useFactor ) {
    solveWorkspace.resizeAndWipe( block.numColumns(), rank );

    MATRIX::copy( solveWorkspace.data(), inputData, block.numColumns(), rank );

    TRACE_ASSERT( node.inPlaceDiagonal(),
                  "Extended variable factor decomposition not supported" );

    node.inPlaceOffDiagonalSolve( block_idx, solveWorkspace.data(),
                                  rank, true /* transpose */ );

    multInput = solveWorkspace.data();
  } else {
    multInput = inputData;
  }
  TIMING_STOP( "diagonalBlockMultiply: in place solve" );

  // Multiply by the original system
  TIMING_START( "diagonalBlockMultiply: base system" );
  node.baseSystemMultiply( A, multInput, rank,
                           block._rowRange, block._columnRange, outputData );
  TIMING_STOP( "diagonalBlockMultiply: base system" );

  diagonalBlockMultiply_nodeDescendents( node, block, block_idx,
                                         dataSizes, rank,
                                         multInput, outputData );

  if ( _useInteriorBlocks ) {
    diagonalBlockMultiply_blockDescendents( A, node, block, block_idx,
                                            rank, multInput, outputData,
                                            nodeList, false );
  }

  // Apply diagonal contributions to the input matrix.  These are
  // contributions that result from the compression of other off-diagonal
  // blocks in this node (and previous nodes)
  TIMING_START( "diagonalBlockMultiply: diagonalContributionMult" );
  node.diagonalContributionMult( block_idx, multInput, rank,
                                 _slackData, _realWorkspaceManager,
                                 outputData );
  TIMING_STOP( "diagonalBlockMultiply: diagonalContributionMult" );

  // If we are doing in place decomposition, we also have to account
  // for any previous in place blocks
  TIMING_START( "diagonalBlockMultiply: multiplyPreviousInPlaceBlocks" );
  if ( node.inPlaceDiagonal() ) {
    node.multiplyPreviousInPlaceBlocks( block_idx, multInput,
                                        block.numColumns(), rank,
                                        outputData,
                                        // Left, non-transposed multiply
                                        true, false );
  }
  TIMING_STOP( "diagonalBlockMultiply: multiplyPreviousInPlaceBlocks" );
}

//////////////////////////////////////////////////////////////////////
// Processes normal descendent multiplications when forming a low-rank
// diagonal block
//////////////////////////////////////////////////////////////////////
void FactorManager::diagonalBlockMultiply_nodeDescendents(
                                                const Supernode &node,
                                                const DenseBlock &block,
                                                int block_idx,
                                                PairArray &dataSizes,
                                                int rank,
                                                const Real *multInput,
                                                Real *outputData )
{
  int                        descendent_idx;
  int                        ancestor_idx;

  // Look through all descendents and multiply by their contributions
  // to this block
  for ( int i = 0; i < _currentDescendents.size(); i++ )
  {
    descendent_idx = _currentDescendents[ i ];

    if ( isInInteriorBlock( descendent_idx ) ) {
      continue;
    }

    Supernode               &descendent = _factor[ descendent_idx ];

    // Get the index of the interaction the descendent has with "node",
    // (its ancestor)
    ancestor_idx = _nextInteractionCache[ descendent_idx ];

    const DenseBlock &interactionRanges
                                  = _diagonalBlockRanges[ i ][ block_idx ];

    if ( !interactionRanges.valid() )
    {
      // This descendent does not actually interact
      // with this particular sub-block of node's diagonal
      continue;
    }

    // Determine workspace sizes
    //
    // sub-matrix workspace
    dataSizes[ 0 ] = IndexPair( range_size( interactionRanges._columnRange ),
                                rank );
    // Initial interaction multiply workspace
    dataSizes[ 1 ] = IndexPair( descendent.numColumns(), rank );
    // Final interaction multiply workspace
    dataSizes[ 2 ] = IndexPair( range_size( interactionRanges._rowRange ),
                                rank );

    RealWorkspace            workspace( _realWorkspaceManager, dataSizes );

    Real                    *subMatrix = workspace.workspaceData( 0 );
    Real                    *multWorkspaceInitial = workspace.workspaceData( 1 );
    Real                    *multWorkspaceFinal = workspace.workspaceData( 2 );

    TIMING_START( "diagonalBlockMultiply: BuildInteractionSubMatrix" );
    Supernode::BuildInteractionSubMatrix( 
                                  descendent.offDiagonal()[ ancestor_idx ],
                                  multInput, subMatrix, rank,
                                  // Range of rows to copy
                                  interactionRanges._columnRange,
                                  // Use an offset, since multInput is only
                                  // as big as it needs to be (it doesn't
                                  // fill the full row space of the interaction)
                                  block._columnRange.first );
    TIMING_STOP( "diagonalBlockMultiply: BuildInteractionSubMatrix" );

    TIMING_START( "diagonalBlockMultiply: interactionMultiply" );
    descendent.interactionMultiply( ancestor_idx, node, block, subMatrix, rank,
                                    interactionRanges._rowRange,
                                    interactionRanges._columnRange,
                                    multWorkspaceInitial,
                                    multWorkspaceFinal );
    TIMING_STOP( "diagonalBlockMultiply: interactionMultiply" );

    // Scatter back to the appropriate workspace
    TIMING_START( "diagonalBlockMultiply: ScatterLowRankUpdate" );
    Supernode::ScatterLowRankUpdate( multWorkspaceFinal, outputData,
                                     descendent.offDiagonal()[ ancestor_idx ],
                                     interactionRanges._rowRange,
                                     rank, block._rowRange.first );
    TIMING_STOP( "diagonalBlockMultiply: ScatterLowRankUpdate" );
  }
}

//////////////////////////////////////////////////////////////////////
// Processes interior block descendent multiplications when forming a
// low-rank diagonal block
//////////////////////////////////////////////////////////////////////
void FactorManager::diagonalBlockMultiply_blockDescendents(
                              const SPARSE_MATRIX::SparseColumnMatrix &A,
                              const Supernode &node,
                              const DenseBlock &block,
                              int block_idx,
                              int rank,
                              const Real *multInput,
                              Real *outputData,
                              IntArray &nodeList,
                              bool transpose )
{
  int                        descendent_idx;
  int                        base_block_interaction;

  // Multiply with each interior block descendent
  for ( int i = 0; i < _currentBlockDescendents.size(); i++ ) {
    descendent_idx = _currentBlockDescendents[ i ];

    const InteriorBlock     &interiorBlock = _interiorBlocks[ descendent_idx ];

    int                      totalColumns;

    // Determine the interaction index in this block which refers to node
    // (it must exist) and the one which applies to the current low rank
    // interaction being processed (doesn't necessarily have to exist)
    base_block_interaction = interiorBlock._interactionMap[ node.nodeID() ];

    TRACE_ASSERT( base_block_interaction >= 0,
                  "No interaction found in interior block" );

    nodeList.clear();

    // Build the list of nodes in this block which contribute to
    // the low rank interaction
    buildInteriorBlockNodeList( descendent_idx, base_block_interaction,
                                base_block_interaction, nodeList );

    TRACE_ASSERT( nodeList.size() > 0, "Zero-length node list" );

    totalColumns = nodeColumns( nodeList );

    // Generate workspaces for interior block multiplication
#if 0
    RealWorkspace            workspace( _realWorkspaceManager,
                                        IndexPair( totalColumns, rank ) );
#endif
    RealWorkspace            workspace( _realWorkspaceManager,
                                        IndexPair( totalColumns, 2 * rank ) );

    // Perform the actual multiplication, subtracting from the
    // output workspace
#if 0
    interiorBlockSchurMultiply( A, nodeList, descendent_idx,
                                base_block_interaction, base_block_interaction,
                                // Input matrix and its size
                                multInput,
                                transpose ? block.numRows() : block.numColumns(),
                                rank,
                                // The inverse row map is just a placeholder here,
                                // since we are providing row/column ranges
                                _inverseRowMap, block.numRows(),
                                workspace.workspaceData( 0 ),
                                outputData,
                                // Row/column ranges
                                block._rowRange, block._columnRange,
                                transpose );
#endif
    interiorBlockSchurMultiply( A, nodeList, descendent_idx,
                                base_block_interaction, base_block_interaction,
                                // Input matrix and its size
                                multInput,
                                transpose ? block.numRows() : block.numColumns(),
                                rank,
                                // The inverse row map is just a placeholder here,
                                // since we are providing row/column ranges
                                _inverseRowMap, block.numRows(),
                                workspace.workspaceData( 0 ),
                                workspace.workspaceData( 0 )
                                  + totalColumns * rank,
                                outputData,
                                // Row/column ranges
                                block._rowRange, block._columnRange,
                                transpose );
  }
}

//////////////////////////////////////////////////////////////////////
// Multiplies data in _decompWorkspace on the left with the
// transpose of each low rank interaction matrix for the current node.
// Puts the result in _decompMultTransWorkspace.
//
// (half a step of power iteration)
//////////////////////////////////////////////////////////////////////
void FactorManager::interactionTransMultiply(
                              const SPARSE_MATRIX::SparseColumnMatrix &A,
                              Supernode &node,
                              bool useFactor,
                              const Real *inputData,
                              Real *outputData,
                              int fixedRank )
{
  int                        nCols = node.numColumns();
  int                        interaction_idx;
  int                        rank;

  IntArray                   nodeList;

  // Workspace sizes
  PairArray                  dataSizes( 3 );

  if ( !inputData ) {
    inputData = _decompWorkspace;
  }

  if ( !outputData ) {
    outputData = _decompMultTransWorkspace;
  }

  // Make space for storing node lists if we're doing interior block multiplies
  if ( _useInteriorBlocks ) {
    nodeList.reserve( _factor.size() );
  }

  for ( int block_idx = 0; block_idx < _lowRankBlocks.size(); block_idx++ )
  {
    if ( !_lowRankBlockActive[ block_idx ] )
    {
      continue;
    }

    interaction_idx = _lowRankBlocks[ block_idx ];

    rank = ( fixedRank > 0 ) ? fixedRank : _blockRanks[ block_idx ];

    const SupernodeInteraction &interaction
                          = node.offDiagonal()[ interaction_idx ];

    TIMING_START( "interactionTransMultiply: base system" );
    node.baseSystemTransMultiply_left( A, interaction_idx, _factor,
                                       inputData, rank, outputData,
                                       true /* clear output */ );
    TIMING_STOP( "interactionTransMultiply: base system" );

    PairArray &descendents = _lowRankDescendents[ block_idx ];

    interactionTransMultiply_nodeDescendents( node, interaction, descendents,
                                              dataSizes, rank,
                                              inputData, outputData );

    if ( _useInteriorBlocks ) {
      interactionMultiply_blockDescendents( A, node, interaction,
                                            rank, inputData, outputData,
                                            nodeList, true );
    }

    // If we are decomposing the factor block, apply the factor
    // inverse here.
    TIMING_START( "interactionTransMultiply: diagonalSolve" );
    if ( useFactor ) {
      node.diagonalSolve( outputData, rank, _slackData );
    }
    TIMING_STOP( "interactionTransMultiply: diagonalSolve" );

    inputData += interaction._rowList.size() * rank;
    outputData += nCols * rank;
  }
}

//////////////////////////////////////////////////////////////////////
// Same as the above, but multiplies with the entire off-diagonal
// all at once
//////////////////////////////////////////////////////////////////////
void FactorManager::interactionTransMultiply_allRows(
                                const SPARSE_MATRIX::SparseColumnMatrix &A,
                                Supernode &node,
                                bool useFactor,
                                const Real *inputData,
                                Real *outputData,
                                int fixedRank )
{
  int                        nCols = node.numColumns();
  int                        rank;

  IntArray                   nodeList;

  // Workspace sizes
  PairArray                  dataSizes( 3 );

  if ( !inputData ) {
    inputData = _decompWorkspace;
  }

  if ( !outputData ) {
    outputData = _decompMultTransWorkspace;
  }

  // Make space for storing node lists if we're doing interior block multiplies
  if ( _useInteriorBlocks ) {
    nodeList.reserve( _factor.size() );
  }

  TRACE_ASSERT( fixedRank > 0 || _blockRanks.size() > 0 );
  rank = ( fixedRank > 0 ) ? fixedRank : _blockRanks[ 0 ];

  TIMING_START( "interactionTransMultiply_allRows: base system" );
  node.baseSystemTransMultiply_left( A, _factor, inputData, rank, outputData );
  TIMING_STOP( "interactionTransMultiply_allRows: base system" );

  interactionMultiply_nodeDescendents( node, dataSizes, rank,
                                       inputData, outputData,
                                       true /* transpose */ );

  if ( _useInteriorBlocks ) {
    interactionMultiply_blockDescendents( A, node, rank,
                                          inputData, outputData,
                                          nodeList, true /* transpose */ );
  }

  // If we are decomposing the factor block, apply the factor
  // inverse here.
  TIMING_START( "interactionTransMultiply_allRows: diagonalSolve" );
  if ( useFactor ) {
    node.diagonalSolve( outputData, rank, _slackData );
  }
  TIMING_STOP( "interactionTransMultiply_allRows: diagonalSolve" );
}

//////////////////////////////////////////////////////////////////////
// Processes normal descendent multiplications when forming a low-rank
// interaction
//////////////////////////////////////////////////////////////////////
void FactorManager::interactionTransMultiply_nodeDescendents(
                                const Supernode &node,
                                const SupernodeInteraction &interaction,
                                const PairArray &descendents,
                                PairArray &dataSizes,
                                int rank,
                                const Real *multInput,
                                Real *outputData )
{
  int                        descendent_idx;
  int                        ancestor_idx;
  int                        desc_interaction_idx;

  // Interaction multiplies with each descendent
  for ( int i = 0; i < descendents.size(); i++ )
  {
    descendent_idx = descendents[ i ].first;

    if ( isInInteriorBlock( descendent_idx ) ) {
      continue;
    }

    desc_interaction_idx = descendents[ i ].second;
    ancestor_idx = _nextInteractionCache[ descendent_idx ];

    Supernode         &descendent = _factor[ descendent_idx ];

    // Determine workspace sizes
    //
    // Sub-matrix workspace
    dataSizes[ 0 ] = IndexPair( interaction._rowList.size(), rank );
    // Initial interaction multiply workspace
    dataSizes[ 1 ] = IndexPair( descendent.numColumns(), rank );
    // Final interaction multiply workspace
    dataSizes[ 2 ] = IndexPair( descendent.interactionRows( ancestor_idx ),
                                rank );

    RealWorkspace         workspace( _realWorkspaceManager, dataSizes );

    Real                  *subMatrix = workspace.workspaceData( 0 );
    Real                  *multWorkspaceInitial = workspace.workspaceData( 1 );
    Real                  *multWorkspaceFinal = workspace.workspaceData( 2 );

    TIMING_START( "interactionTransMultiply: interactionTransMultiply_left" );
    descendent.interactionTransMultiply_left( desc_interaction_idx,
                                              ancestor_idx,
                                              interaction._rowList,
                                              node,
                                              multInput, rank,
                                              subMatrix,
                                              multWorkspaceInitial,
                                              multWorkspaceFinal );
    TIMING_STOP( "interactionTransMultiply: interactionTransMultiply_left" );

    // Scatter back to the desired row space
    TIMING_START( "interactionTransMultiply: ScatterLowRankTransUpdate_row" );
    Supernode::ScatterLowRankTransUpdate_row( ancestor_idx,
                                              descendent, node,
                                              multWorkspaceFinal,
                                              outputData, rank );
    TIMING_STOP( "interactionTransMultiply: ScatterLowRankTransUpdate_row" );
  }
}

#if 0
//////////////////////////////////////////////////////////////////////
// Processes interior block descendent multiplications when forming
// a low-rank interaction
//////////////////////////////////////////////////////////////////////
void FactorManager::interactionTransMultiply_blockDescendents(
                                const SPARSE_MATRIX::SparseColumnMatrix &A,
                                const Supernode &node,
                                const SupernodeInteraction &interaction,
                                int rank,
                                const Real *multInput,
                                Real *outputData,
                                IntArray &nodeList )
{
}
#endif

//////////////////////////////////////////////////////////////////////
// Similar to the above function, but does the multiplication with
// a sub-matrix of the main diagonal (indexed by block_idx)
//////////////////////////////////////////////////////////////////////
void FactorManager::diagonalBlockTransMultiply(
                              const SPARSE_MATRIX::SparseColumnMatrix &A,
                              Supernode &node,
                              bool useFactor,
                              int block_idx,
                              int rank,
                              const Real *inputData,
                              Real *outputData,
                              bool write )
{
  const DenseBlock &block = node.diagonalLowRankBlocks()[ block_idx ];

  int                        descendent_idx;
  int                        ancestor_idx;

  IntArray                   nodeList;

  // Workspace sizes
  PairArray                  dataSizes( 3 );

  if ( !inputData ) {
    //inputData = _decompMultTransWorkspace;
    inputData = _decompWorkspace;
  }

  if ( !outputData ) {
    //outputData = _decompWorkspace;
    outputData = _decompMultTransWorkspace;
  }

  // Make space for storing node lists if we're doing interior block multiplies
  if ( _useInteriorBlocks ) {
    nodeList.reserve( _factor.size() );
  }

#if 0
  {
    MATRIX::write( inputData, block.numRows(), rank,
                   "super_numeric/diag_multInput.matrix" );
  }
#endif

  // Multiply by the original system
  TIMING_START( "diagonalBlockTransMultiply: base system" );
  node.baseSystemTransMultiply_left( A, inputData, rank,
                                     block._rowRange, block._columnRange,
                                     outputData );
  TIMING_STOP( "diagonalBlockTransMultiply: base system" );

  diagonalBlockTransMultiply_nodeDescendents( node, block, block_idx,
                                              dataSizes, rank,
                                              inputData, outputData );

  if ( _useInteriorBlocks ) {
    diagonalBlockMultiply_blockDescendents( A, node, block, block_idx,
                                            rank, inputData, outputData,
                                            nodeList, true );
  }

#if 0
  // Look through all descendents and multiply by their contributions
  // to this block
  for ( int i = 0; i < _currentDescendents.size(); i++ )
  {
    descendent_idx = _currentDescendents[ i ];

    Supernode               &descendent = _factor[ descendent_idx ];

    // Get the index of the interaction the descendent has with "node",
    // (its ancestor)
    ancestor_idx = _nextInteractionCache[ descendent_idx ];

    const DenseBlock &interactionRanges
                                  = _diagonalBlockRanges[ i ][ block_idx ];

    if ( !interactionRanges.valid() )
    {
      // This descendent does not actually interact
      // with this particular sub-block of node's diagonal
      continue;
    }

    // Determine workspace sizes
    //
    // Sub-matrix workspace
    dataSizes[ 0 ] = IndexRange( range_size( interactionRanges._rowRange ),
                                 rank );
    // Initial interaction multiply workspace
    dataSizes[ 1 ] = IndexRange( descendent.numColumns(), rank );
    // Final interaction multiply workspace
    dataSizes[ 2 ] = IndexRange( range_size( interactionRanges._columnRange ),
                                 rank );

    RealWorkspace         workspace( _realWorkspaceManager, dataSizes );

    Real                  *subMatrix = workspace.workspaceData( 0 );
    Real                  *multWorkspaceInitial = workspace.workspaceData( 1 );
    Real                  *multWorkspaceFinal = workspace.workspaceData( 2 );

    TIMING_START( "diagonalBlockTransMultiply: interactionTransMultiply" );
    descendent.interactionTransMultiply_left( ancestor_idx, node, block,
                                              inputData, rank,
                                              interactionRanges._rowRange,
                                              interactionRanges._columnRange,
                                              subMatrix,
                                              multWorkspaceInitial,
                                              multWorkspaceFinal );
    TIMING_STOP( "diagonalBlockTransMultiply: interactionTransMultiply" );

    // FIXME: I *think* this version of this function works here
    TIMING_START( "diagonalBlockTransMultiply: ScatterLowRankUpdate" );
    Supernode::ScatterLowRankUpdate( multWorkspaceFinal, outputData,
                                     descendent.offDiagonal()[ ancestor_idx ],
                                     // Row range instead of column range
                                     // due to transposition
                                     interactionRanges._columnRange, rank,
                                     block._columnRange.first );
    TIMING_STOP( "diagonalBlockTransMultiply: ScatterLowRankUpdate" );
  }
#endif

  // Apply diagonal contributions to the input matrix.  These are
  // contributions that result from the compression of other off-diagonal
  // blocks in this node (and previous nodes)
  TIMING_START( "diagonalBlockTransMultiply: diagonalContributionMult" );
  node.diagonalContributionMult( block_idx, inputData, rank,
                                 _slackData, _realWorkspaceManager,
                                 outputData, true /* transpose */ );
  TIMING_STOP( "diagonalBlockTransMultiply: diagonalContributionMult" );

  // If we are doing in place decomposition, we also have to account for
  // any previous in place blocks
  TIMING_START( "diagonalBlockTransMultiply: multiplyPreviousInPlaceBlocks" );
  if ( node.inPlaceDiagonal() ) {
    node.multiplyPreviousInPlaceBlocks( block_idx, inputData,
                                        block.numRows(), rank,
                                        outputData,
                                        // Left, transposed multiply
                                        true, true );
  }
  TIMING_STOP( "diagonalBlockTransMultiply: multiplyPreviousInPlaceBlocks" );

  // If we are decomposing the factor block, apply the factor
  // inverse here
  TIMING_START( "diagonalBlockTransMultiply: inPlaceOffDiagonalSolve" );
  if ( useFactor ) {
    TRACE_ASSERT( node.inPlaceDiagonal(),
                  "Factor decomposition with extended nodes not supported" );

    node.inPlaceOffDiagonalSolve( block_idx, outputData, rank,
                                  false /* no transpose */ );
  }
  TIMING_STOP( "diagonalBlockTransMultiply: inPlaceOffDiagonalSolve" );

#if 0
  {
    MATRIX::write( outputData, block.numColumns(), rank,
                   "super_numeric/diag_multOutput.matrix" );
    abort();
  }
#endif
}

//////////////////////////////////////////////////////////////////////
// Processes normal descendent multiplications when forming a low-rank
// diagonal block
//////////////////////////////////////////////////////////////////////
void FactorManager::diagonalBlockTransMultiply_nodeDescendents(
                                                const Supernode &node,
                                                const DenseBlock &block,
                                                int block_idx,
                                                PairArray &dataSizes,
                                                int rank,
                                                const Real *multInput,
                                                Real *outputData )
{
  int                        descendent_idx;
  int                        ancestor_idx;

  // Look through all descendents and multiply by their contributions
  // to this block
  for ( int i = 0; i < _currentDescendents.size(); i++ )
  {
    descendent_idx = _currentDescendents[ i ];

    if ( isInInteriorBlock( descendent_idx ) ) {
      continue;
    }

    Supernode               &descendent = _factor[ descendent_idx ];

    // Get the index of the interaction the descendent has with "node",
    // (its ancestor)
    ancestor_idx = _nextInteractionCache[ descendent_idx ];

    const DenseBlock &interactionRanges
                                  = _diagonalBlockRanges[ i ][ block_idx ];

    if ( !interactionRanges.valid() )
    {
      // This descendent does not actually interact
      // with this particular sub-block of node's diagonal
      continue;
    }

    // Determine workspace sizes
    //
    // Sub-matrix workspace
    dataSizes[ 0 ] = IndexRange( range_size( interactionRanges._rowRange ),
                                 rank );
    // Initial interaction multiply workspace
    dataSizes[ 1 ] = IndexRange( descendent.numColumns(), rank );
    // Final interaction multiply workspace
    dataSizes[ 2 ] = IndexRange( range_size( interactionRanges._columnRange ),
                                 rank );

    RealWorkspace         workspace( _realWorkspaceManager, dataSizes );

    Real                  *subMatrix = workspace.workspaceData( 0 );
    Real                  *multWorkspaceInitial = workspace.workspaceData( 1 );
    Real                  *multWorkspaceFinal = workspace.workspaceData( 2 );

    TIMING_START( "diagonalBlockTransMultiply: interactionTransMultiply" );
    descendent.interactionTransMultiply_left( ancestor_idx, node, block,
                                              multInput, rank,
                                              interactionRanges._rowRange,
                                              interactionRanges._columnRange,
                                              subMatrix,
                                              multWorkspaceInitial,
                                              multWorkspaceFinal );
    TIMING_STOP( "diagonalBlockTransMultiply: interactionTransMultiply" );

    // FIXME: I *think* this version of this function works here
    TIMING_START( "diagonalBlockTransMultiply: ScatterLowRankUpdate" );
    Supernode::ScatterLowRankUpdate( multWorkspaceFinal, outputData,
                                     descendent.offDiagonal()[ ancestor_idx ],
                                     // Row range instead of column range
                                     // due to transposition
                                     interactionRanges._columnRange, rank,
                                     block._columnRange.first );
    TIMING_STOP( "diagonalBlockTransMultiply: ScatterLowRankUpdate" );
  }
}

//////////////////////////////////////////////////////////////////////
// Forms the second part of a low rank interaction by multiplying
// the interaction on the left with its basis transposed
//////////////////////////////////////////////////////////////////////
void FactorManager::fullBasisProject(
                              const SPARSE_MATRIX::SparseColumnMatrix &A,
                              Supernode &node,
                              bool writeFactors,
                              int fixedRank,
                              bool transposeBasis,
                              bool transpose,
                              bool useFactor,
                              OffDiagonalCompressionType basisType )
{
  if ( transpose ) {
    switch ( basisType ) {
      case COMPRESS_FULL_OFF_DIAGONAL:
        fullBasisProject_transpose_allRows( A, node, writeFactors,
                                            fixedRank, transposeBasis,
                                            useFactor );
        break;
      case COMPRESS_INDIVIDUAL_INTERACTIONS:
      default:
        fullBasisProject_transpose( A, node, writeFactors,
                                    fixedRank, transposeBasis,
                                    useFactor );
        break;
    }
  }
  else {
    switch ( basisType ) {
      case COMPRESS_FULL_OFF_DIAGONAL:
        fullBasisProject_standard_allRows( A, node, writeFactors,
                                           fixedRank, transposeBasis,
                                           useFactor );
        break;
      case COMPRESS_INDIVIDUAL_INTERACTIONS:
      default:
        fullBasisProject_standard( A, node, writeFactors,
                                   fixedRank, transposeBasis,
                                   useFactor );
        break;
    }
  }
}

//////////////////////////////////////////////////////////////////////
// Version of fullBasisProject for decomposing lower triangular
// blocks directly
//////////////////////////////////////////////////////////////////////
void FactorManager::fullBasisProject_standard(
                                const SPARSE_MATRIX::SparseColumnMatrix &A,
                                Supernode &node,
                                bool writeFactors,
                                int fixedRank,
                                bool transposeBasis,
                                bool useFactor )
{
  int                        interaction_idx;
  int                        descendent_idx;
  int                        desc_interaction_idx;
  int                        ancestor_idx;
  int                        nCols = node.numColumns();
  int                        nRows1;
  int                        rank;
  const Real                *Q;
  const Real                *U;
#if 0
  Real                      *Q;
  Real                      *U;
#endif

  // Workspace sizes
  PairArray                  dataSizes( 3 );

  TRACE_ASSERT( transposeBasis, "Not implemented" );

  // FIXME: debugging another way of doing this
  // Make all blocks active so that we can use interactionTransMultiply
#ifdef DO_TIMING
  _timers[ BASIS_PROJECT ].tick();
#endif
  activateAllBlocks();
  interactionTransMultiply( A, node, useFactor,
                            // Standard workspaces, since bases are already
                            // stored in _decompWorkspace
                            NULL, NULL,
                            // Use ranks in _blockRanks
                            0 );
#ifdef DO_TIMING
  _timers[ BASIS_PROJECT ].tock();
#endif

  // Try to make things positive definite
  //reprojectDecomposition( A, node );

  Q = _decompWorkspace;
  U = _decompMultTransWorkspace;

  for ( int block_idx = 0; block_idx < _lowRankBlocks.size(); block_idx++ )
  {
    interaction_idx = _lowRankBlocks[ block_idx ];

    rank = ( fixedRank > 0 ) ? fixedRank : _blockRanks[ block_idx ];

    const SupernodeInteraction &interaction
                        = node.offDiagonal()[ interaction_idx ];

    // FIXME: need to do this for now
    MATRIX Utrans( rank, nCols );
    MATRIX::transposeBLAS( Utrans.data(), U, nCols, rank );

    // Get a workspace for expanding V
    RealWorkspace workspace(
              _realWorkspaceManager,
              IndexPair( _factor[ interaction._nodeID ].numColumns(), rank ) );

#ifdef DO_TIMING
    _timers[ ASSIGN_BASIS ].tick();
#endif
    node.assignExtraData( interaction_idx, _factor, Q /* V */,
                          Utrans.data(),
                          rank, _slackData,
                          _slackDataOffset, _availableSlackData,
                          workspace.workspaceData( 0 ) /* To expand V */ );
#ifdef DO_TIMING
    _timers[ ASSIGN_BASIS ].tock();
#endif

    TRACE_ASSERT( _availableSlackData >= 0, "Ran out of extra data" );

    Q += interaction._rowList.size() * rank;
    U += nCols * rank;
  }
}

//////////////////////////////////////////////////////////////////////
// Version of the above to use when decomposing all rows at once
//////////////////////////////////////////////////////////////////////
void FactorManager::fullBasisProject_standard_allRows(
                                const SPARSE_MATRIX::SparseColumnMatrix &A,
                                Supernode &node,
                                bool writeFactors,
                                int fixedRank,
                                bool transposeBasis,
                                bool useFactor )
{
  int                        nCols = node.numColumns();
  int                        rank;
  const Real                *Q;
  const Real                *U;

  TRACE_ASSERT( transposeBasis, "Not implemented" );

  // Make all blocks active so that we can use interactionTransMultiply
#ifdef DO_TIMING
  _timers[ BASIS_PROJECT ].tick();
#endif
  activateAllBlocks();
  interactionTransMultiply_allRows( A, node, useFactor,
                                    // Standard workspaces, since bases are
                                    // already stored in _decompWorkspace
                                    NULL, NULL,
                                    // Use ranks in _blockRanks
                                    0 );
#ifdef DO_TIMING
  _timers[ BASIS_PROJECT ].tock();
#endif

  Q = _decompWorkspace;
  U = _decompMultTransWorkspace;

  rank = ( fixedRank > 0 ) ? fixedRank : _blockRanks[ 0 ];

  MATRIX Utrans( rank, nCols );
  MATRIX::transposeBLAS( Utrans.data(), U, nCols, rank );

#ifdef DO_TIMING
  _timers[ ASSIGN_BASIS ].tick();
#endif
  node.assignExtraData( Q /* V */, Utrans.data(), rank, _slackData,
                        _slackDataOffset, _availableSlackData );
#ifdef DO_TIMING
  _timers[ ASSIGN_BASIS ].tock();
#endif

  TRACE_ASSERT( _availableSlackData >= 0, "Ran out of extra data" );
}

//////////////////////////////////////////////////////////////////////
// Version of fullBasisProject for decomposing transposes of lower
// triangular blocks (ie. decomposing upper triangular blocks)
//////////////////////////////////////////////////////////////////////
void FactorManager::fullBasisProject_transpose(
                                const SPARSE_MATRIX::SparseColumnMatrix &A,
                                Supernode &node,
                                bool writeFactors,
                                int fixedRank,
                                bool transposeBasis,
                                bool useFactor )
{
  int                        interaction_idx;
  int                        nCols = node.numColumns();
  int                        rank;
  const Real                *Q;
  const Real                *U;
#if 0
  Real                      *Q;
  Real                      *U;
#endif

  TRACE_ASSERT( transposeBasis, "Not implemented" );

#if 0
  Q = _decompWorkspace;
  U = _decompMultTransWorkspace;
  for ( int block_idx = 0; block_idx < _lowRankBlocks.size(); block_idx++ )
  {
    interaction_idx = _lowRankBlocks[ block_idx ];

    rank = ( fixedRank > 0 ) ? fixedRank : _blockRanks[ block_idx ];

    const SupernodeInteraction &interaction
                        = node.offDiagonal()[ interaction_idx ];

    char buf[ 1024 ];
    if ( _useInteriorBlocks ) {
      sprintf( buf, "node_%d_interaction_%d_prefinal_int.matrix",
               node.nodeID(), _lowRankBlocks[ block_idx ] );
    } else {
      sprintf( buf, "node_%d_interaction_%d_prefinal.matrix",
               node.nodeID(), _lowRankBlocks[ block_idx ] );
    }

    MATRIX::write( U, nCols, rank, buf );

    Q += interaction._rowList.size() * rank;
    U += nCols * rank;
  }
#endif

  // FIXME: debugging another way of doing this
  // Make all blocks active so that we can use interactionTransMultiply
#ifdef DO_TIMING
  _timers[ BASIS_PROJECT ].tick();
#endif
  activateAllBlocks();

  interactionMultiply( A, node, useFactor,
                       // Standard workspaces, since bases are already
                       // stored in _decompMultTransWorkspace
                       NULL, NULL,
                       // Use ranks in _blockRanks
                       0 );
#ifdef DO_TIMING
  _timers[ BASIS_PROJECT ].tock();
#endif

#if 0
  Q = _decompWorkspace;
  U = _decompMultTransWorkspace;
  for ( int block_idx = 0; block_idx < _lowRankBlocks.size(); block_idx++ )
  {
    interaction_idx = _lowRankBlocks[ block_idx ];

    rank = ( fixedRank > 0 ) ? fixedRank : _blockRanks[ block_idx ];

    const SupernodeInteraction &interaction
                        = node.offDiagonal()[ interaction_idx ];

    char buf[ 1024 ];
    if ( _useInteriorBlocks ) {
      sprintf( buf, "node_%d_interaction_%d_prereproj_int.matrix",
               node.nodeID(), _lowRankBlocks[ block_idx ] );
    } else {
      sprintf( buf, "node_%d_interaction_%d_prereproj.matrix",
               node.nodeID(), _lowRankBlocks[ block_idx ] );
    }

    MATRIX::write( Q, interaction._rowList.size(), rank, buf );

    if ( _useInteriorBlocks ) {
      sprintf( buf, "node_%d_interaction_%d_prereproj_U_int.matrix",
               node.nodeID(), _lowRankBlocks[ block_idx ] );
    } else {
      sprintf( buf, "node_%d_interaction_%d_prereproj_U.matrix",
               node.nodeID(), _lowRankBlocks[ block_idx ] );
    }

    MATRIX::write( U, nCols, rank, buf );

    Q += interaction._rowList.size() * rank;
    U += nCols * rank;
  }
#endif

  // Try to make things positive definite
  reprojectDecomposition( A, node );

#if 0
  Q = _decompWorkspace;
  U = _decompMultTransWorkspace;
  for ( int block_idx = 0; block_idx < _lowRankBlocks.size(); block_idx++ )
  {
    interaction_idx = _lowRankBlocks[ block_idx ];

    rank = ( fixedRank > 0 ) ? fixedRank : _blockRanks[ block_idx ];

    const SupernodeInteraction &interaction
                        = node.offDiagonal()[ interaction_idx ];

    char buf[ 1024 ];
    if ( _useInteriorBlocks ) {
      sprintf( buf, "node_%d_interaction_%d_postreproj_int.matrix",
               node.nodeID(), _lowRankBlocks[ block_idx ] );
    } else {
      sprintf( buf, "node_%d_interaction_%d_postreproj.matrix",
               node.nodeID(), _lowRankBlocks[ block_idx ] );
    }

    MATRIX::write( Q, interaction._rowList.size(), rank, buf );

    if ( _useInteriorBlocks ) {
      sprintf( buf, "node_%d_interaction_%d_postreproj_U_int.matrix",
               node.nodeID(), _lowRankBlocks[ block_idx ] );
    } else {
      sprintf( buf, "node_%d_interaction_%d_postreproj_U.matrix",
               node.nodeID(), _lowRankBlocks[ block_idx ] );
    }

    MATRIX::write( U, nCols, rank, buf );

    Q += interaction._rowList.size() * rank;
    U += nCols * rank;
  }
#endif

  Q = _decompWorkspace;
  U = _decompMultTransWorkspace;

  for ( int block_idx = 0; block_idx < _lowRankBlocks.size(); block_idx++ )
  {
    interaction_idx = _lowRankBlocks[ block_idx ];

    rank = ( fixedRank > 0 ) ? fixedRank : _blockRanks[ block_idx ];

    const SupernodeInteraction &interaction
                        = node.offDiagonal()[ interaction_idx ];

    // FIXME: need to do this for now
    MATRIX Utrans( rank, nCols );
    MATRIX::transposeBLAS( Utrans.data(), U, nCols, rank );

    // Get a workspace for expanding V
    RealWorkspace workspace(
              _realWorkspaceManager,
              IndexPair( _factor[ interaction._nodeID ].numColumns(), rank ) );

#ifdef DO_TIMING
    _timers[ ASSIGN_BASIS ].tick();
#endif
    node.assignExtraData( interaction_idx, _factor, Q /* V */,
                          Utrans.data(),
                          rank, _slackData,
                          _slackDataOffset, _availableSlackData,
                          workspace.workspaceData( 0 ) /* To expand V */ );
#ifdef DO_TIMING
    _timers[ ASSIGN_BASIS ].tock();
#endif

    TRACE_ASSERT( _availableSlackData >= 0, "Ran out of extra data" );

#if 0
    {
      char buf[ 1024 ];

      if ( _useInteriorBlocks ) {
        sprintf( buf, "V_node_%d_interaction_%d_int.matrix", node.nodeID(),
                 interaction_idx );
      } else {
        sprintf( buf, "V_node_%d_interaction_%d.matrix", node.nodeID(),
                 interaction_idx );
      }

      MATRIX::write( Q, interaction._rowList.size(), rank, buf );

      if ( _useInteriorBlocks ) {
        sprintf( buf, "U_node_%d_interaction_%d_int.matrix", node.nodeID(),
                 interaction_idx );
      } else {
        sprintf( buf, "U_node_%d_interaction_%d.matrix", node.nodeID(),
                 interaction_idx );
      }

      MATRIX::write( U, nCols, rank, buf );
    }
#endif

    Q += interaction._rowList.size() * rank;
    U += nCols * rank;
  }
}

//////////////////////////////////////////////////////////////////////
// Version to use when decomposing all rows at once
//////////////////////////////////////////////////////////////////////
void FactorManager::fullBasisProject_transpose_allRows(
                                const SPARSE_MATRIX::SparseColumnMatrix &A,
                                Supernode &node,
                                bool writeFactors,
                                int fixedRank,
                                bool transposeBasis,
                                bool useFactor )
{
  int                        nCols = node.numColumns();
  int                        rank;
  const Real                *Q;
  const Real                *U;

  TRACE_ASSERT( transposeBasis, "Not implemented" );

  // Make all blocks active so that we can use interactionMultiply
#ifdef DO_TIMING
  _timers[ BASIS_PROJECT ].tick();
#endif
  activateAllBlocks();

#if 0
  // FIXME: debugging
  {
    cout << SDUMP( node.nodeID() ) << endl;
    MATRIX::write( _decompMultTransWorkspace, nCols,
                   ( fixedRank > 0 ) ? fixedRank : _blockRanks[ 0 ],
                   "U.matrix" );
    abort();
  }
#endif

  interactionMultiply_allRows( A, node, useFactor,
                               // Standard workspaces, since bases are already
                               // stored in _decompMultTransWorkspace
                               NULL, NULL,
                               // Use ranks in _blockRanks
                               0 );
#ifdef DO_TIMING
  _timers[ BASIS_PROJECT ].tock();
#endif

  Q = _decompWorkspace;
  U = _decompMultTransWorkspace;

  rank = ( fixedRank > 0 ) ? fixedRank : _blockRanks[ 0 ];

  MATRIX Utrans( rank, nCols );
  MATRIX::transposeBLAS( Utrans.data(), U, nCols, rank );

#ifdef DO_TIMING
  _timers[ ASSIGN_BASIS ].tick();
#endif
  node.assignExtraData( Q /* V */, Utrans.data(), rank, _slackData,
                        _slackDataOffset, _availableSlackData );
#ifdef DO_TIMING
  _timers[ ASSIGN_BASIS ].tock();
#endif

  TRACE_ASSERT( _availableSlackData >= 0, "Ran out of extra data" );
}

//////////////////////////////////////////////////////////////////////
// Reprojection procedure to guarantee that the modified system
// remains positive definite.
//////////////////////////////////////////////////////////////////////
void FactorManager::reprojectDecomposition(
                          const SPARSE_MATRIX::SparseColumnMatrix &A,
                          Supernode &node )
{
  // Currently, the (orthonormal) "V" parts of the decompositions
  // are stored in _decompWorkspace, and the "U" parts are stored
  // in _decompMultTransWorkspace
  
  // Workspace for building a global "U" basis
  MATRIX                     Ubasis;
  MATRIX                     Ucopy;
  MATRIX                     Uoutput;
  MATRIX                     innerProduct;
  MATRIX                     qrData;

  int                        nCols = node.numColumns();
  int                        totalRank = 0;
  int                        maxRank = 0;
  int                        totalRows = 0;
  int                        totalCols;
  int                        rank;

  const Real                *baseInputData = _decompMultTransWorkspace;
  Real                      *baseData;

  for ( int block_idx = 0; block_idx < _lowRankBlocks.size(); block_idx++ ) {
    totalRank += _blockRanks[ block_idx ];
  }

  Ubasis.resizeAndWipe( nCols, totalRank );
  Uoutput.resizeAndWipe( nCols, totalRank );
  innerProduct.resizeAndWipe( totalRank, totalRank );
  qrData.resizeAndWipe( max( nCols, totalRank ), 1 );

  // Copy all "U" components to the same matrix
  baseData = Ubasis.data();
  for ( int block_idx = 0; block_idx < _lowRankBlocks.size(); block_idx++ ) {
    const SupernodeInteraction &interaction
                          = node.offDiagonal()[ _lowRankBlocks[ block_idx ] ];

    totalRows += interaction._rowList.size();

    rank = _blockRanks[ block_idx ];

    maxRank = max( rank, maxRank );

#if 0
    printf( "Block %d is %d x %d\n",
            block_idx, (int)interaction._rowList.size(), rank );
#endif

    MATRIX::copy( baseData, baseInputData, nCols, rank,
                  // Leading dimension of output
                  totalRank,
                  // Leading dimension of input
                  rank );

    baseInputData += nCols * rank;
    baseData += rank;
  }
  Ucopy = Ubasis;

  totalCols = nCols * _lowRankBlocks.size();

#if 0
  {
    char buf[ 1024 ];
    if ( _useInteriorBlocks ) {
      sprintf( buf, "fullBasisSet_node_%d_int.matrix", node.nodeID() );
    } else {
      sprintf( buf, "fullBasisSet_node_%d.matrix", node.nodeID() );
    }

    Ucopy.write( buf );
  }
#endif

#if 0
  Ucopy.write( "super_numeric/Ucopy.matrix" );
#endif

#if 0
  printf( "Generating orthonormal basis of size %d x %d\n",
          Ubasis.rows(), Ubasis.cols() );
#endif

  // Form an orthonormal basis from the "U" components of the decomposition
  MATRIX::qr( Ubasis.data(), Ubasis.rows(), Ubasis.cols(), qrData.data(),
              NULL, -1 );
  MATRIX::extractQRfactor( Ubasis.data(), Ubasis.rows(), Ubasis.cols(),
                           qrData.data(), NULL, -1 );

#if 0
  Ubasis.write( "super_numeric/Ubasis.matrix" );
#endif

  // Make this the basis for all "U" components
  if ( totalRows * totalRank > _decompWorkspaceSz ) {
    cout << SDUMP( totalRows ) << "    " << SDUMP( totalRank ) << endl;
    cout << SDUMP( totalRows ) << "    " << SDUMP( totalRank ) << endl;
    cout << SDUMP( totalRows ) << "    " << SDUMP( totalRank ) << endl;
    cout << SDUMP( totalRows ) << "    " << SDUMP( totalRank ) << endl;
    cout << SDUMP( totalRows ) << "    " << SDUMP( totalRank ) << endl;
    cout << SDUMP( totalRows ) << "    " << SDUMP( totalRank ) << endl;
    cout << SDUMP( totalRows ) << "    " << SDUMP( totalRank ) << endl;
    cout << SDUMP( totalRows * totalRank ) << endl;
    cout << SDUMP( _decompWorkspaceSz ) << endl;
  }

  if ( totalCols * totalRank >= _decompMultTransWorkspaceSz ) {
    cout << SDUMP( totalCols ) << "    " << SDUMP( totalRank ) << endl;
    cout << SDUMP( totalCols * totalRank ) << endl;
    cout << SDUMP( _decompMultTransWorkspaceSz ) << endl;
  }

  TRACE_ASSERT( totalCols * totalRank <= _decompMultTransWorkspaceSz );
  TRACE_ASSERT( totalRows * totalRank <= _decompWorkspaceSz );
  baseData = _decompMultTransWorkspace;
  for ( int block_idx = 0; block_idx < _lowRankBlocks.size(); block_idx++ ) {
    const SupernodeInteraction  &interaction
                        = node.offDiagonal()[ _lowRankBlocks[ block_idx ] ];

    // Change the rank of the block - everything now has the same rank
    _blockRanks[ block_idx ] = totalRank;

    MATRIX::copy( baseData, Ubasis.data(), nCols, totalRank );

    baseData += nCols * totalRank;
  }

  // Project to get the new "V" components
  interactionMultiply( A, node, true /* use factor */,
                       NULL, NULL, -1 );

#if 0
  // Try multiplying every block by the new basis
  MATRIX multiplyInputWorkspace( _lowRankBlocks.size() * nCols, totalRank );
  MATRIX multiplyOutputWorkspace( totalRows, totalRank );
  baseData = multiplyInputWorkspace.data();
  for ( int block_idx = 0; block_idx < _lowRankBlocks.size(); block_idx++ ) {
    MATRIX::copy( baseData, Ubasis.data(), nCols, totalRank );

    baseData += nCols * totalRank;
  }

  interactionMultiply( A, node, true /* use factor */,
                       // Override inputs/outputs
                       multiplyInputWorkspace.data(),
                       multiplyOutputWorkspace.data(),
                       // Use a fixed rank
                       totalRank );

  multiplyOutputWorkspace.write( "super_numeric/reprojection.matrix" );
#endif

#if 0
  //////////////////////////////////////////////////////////////////////
  // Let's try writing V as one tall, thin, orthonormal basis
  //////////////////////////////////////////////////////////////////////
  MATRIX                     Vbasis( totalRows, maxRank );

  printf( "Generating orthonormal V basis; max rank = %d\n", maxRank );

  baseInputData = _decompWorkspace;
  baseData = Vbasis.data();

  for ( int block_idx = 0; block_idx < _lowRankBlocks.size(); block_idx++ ) {
    int                          interaction_idx = _lowRankBlocks[ block_idx ];
    const SupernodeInteraction  &interaction
                                    = node.offDiagonal()[ interaction_idx ];

    rank = _blockRanks[ block_idx ];

    MATRIX::copy( baseData, baseInputData,
                  interaction._rowList.size(), rank,
                  // Input leading dimension
                  maxRank,
                  // Output leading dimension
                  rank );

    baseData += interaction._rowList.size() * maxRank;
    baseInputData += interaction._rowList.size() * rank;
  }

  // Turn this in to a basis
  qrData.resizeAndWipe( max( totalRank, maxRank ), 1 );
  MATRIX::qr( Vbasis.data(), totalRows, maxRank, qrData.data(),
              NULL, -1 );
  MATRIX::extractQRfactor( Vbasis.data(), totalRows, maxRank, qrData.data(),
                           NULL, -1 );

  // Copy this to the decomposition workspace
  TRACE_ASSERT( totalRows * maxRank <= _decompWorkspaceSz );
  MATRIX::copy( _decompWorkspace, Vbasis.data(), totalRows, maxRank );

  MATRIX                     projection( _lowRankBlocks.size() * nCols,
                                         maxRank );

  interactionTransMultiply( A, node, true /* use factor */,
                            _decompWorkspace, projection.data(),
                            // Fixed rank for everything
                            maxRank );

  // Now, sum these up to make a single contribution
  baseInputData = projection.data() + nCols * maxRank;
  for ( int block_idx = 1; block_idx < _lowRankBlocks.size(); block_idx++ ) {
    MATRIX::axpy( projection.data(), baseInputData, nCols, maxRank );

    baseInputData += nCols * maxRank;
  }

  // Turn this in to a basis
  qrData.resizeAndWipe( max( nCols, maxRank ), 1 );
  MATRIX::qr( projection.data(), nCols, maxRank, qrData.data(),
              NULL, -1 );
  MATRIX::extractQRfactor( projection.data(), nCols, maxRank, qrData.data(),
                           NULL, -1 );

  // The "U" we want is now stored in the leading block of projection.data()
  // so overwrite _decompMultTransWorkspace
  TRACE_ASSERT(
    projection.rows() * projection.cols() <= _decompMultTransWorkspaceSz );

  baseData = _decompMultTransWorkspace;
  for ( int block_idx = 0; block_idx < _lowRankBlocks.size(); block_idx++ ) {
    MATRIX::copy( baseData, projection.data(), nCols, maxRank );

    // Set the rank as well
    _blockRanks[ block_idx ] = maxRank;

    baseData += nCols * maxRank;
  }

  // Do one more multiply to form the proper projection in to this
  // new "U" basis
  interactionMultiply( A, node, true /* use factor */,
                       NULL, NULL, -1 );
#endif
}

//////////////////////////////////////////////////////////////////////
// Same as the above function, b
//////////////////////////////////////////////////////////////////////
void FactorManager::diagonalBasisProject(
                          const SPARSE_MATRIX::SparseColumnMatrix &A,
                          Supernode &node, int block_idx,
                          bool writeFactors,
                          bool transpose,
                          bool useFactor,
                          int fixedRank )
{
  if ( transpose ) {
    diagonalBasisProject_transpose( A, node, block_idx, writeFactors,
                                    useFactor, fixedRank );
  } else {
    diagonalBasisProject_standard( A, node, block_idx, writeFactors,
                                   useFactor, fixedRank );
  }
}

//////////////////////////////////////////////////////////////////////
// Version of diagonalBasisProject for decomposing lower triangular
// blocks directly
//////////////////////////////////////////////////////////////////////
void FactorManager::diagonalBasisProject_standard(
                          const SPARSE_MATRIX::SparseColumnMatrix &A,
                          Supernode &node, int block_idx,
                          bool writeFactors,
                          bool useFactor,
                          int fixedRank )
{
  int                        descendent_idx;
  int                        ancestor_idx;
  int                        rank;

  // Workspace sizes
  PairArray                  dataSizes( 3 );

  Real                      *Q = _decompWorkspace;
  Real                      *U = _decompMultTransWorkspace;

  const DenseBlock &block = node.diagonalLowRankBlocks()[ block_idx ];

  rank = ( fixedRank > 0 ) ? fixedRank : _blockRanks[ 0 ];

  diagonalBlockTransMultiply( A, node, useFactor, block_idx,
                              rank, NULL, NULL );

  MATRIX                     Utrans( rank, block.numColumns() );

  MATRIX::transposeBLAS( Utrans.data(), U, block.numColumns(), rank );

  // Copy V and Utrans data in to the appropriate interaction, allocate
  // data for new node, etc.
  node.assignExtraDiagonalData( block_idx, _factor, Q /* V */,
                                Utrans.data() /* Utrans */,
                                rank, _slackData,
                                _slackDataOffset, _availableSlackData );

  TRACE_ASSERT( _availableSlackData >= 0, "Ran out of extra data" );

#if 0
  // Multiply the transpose of Q by the original system submatrix
  node.baseSystemTransMultiply( A, Q, rank,
                                block._rowRange, block._columnRange,
                                _decompMultTransWorkspace );

  // Look through all descendents and multiply by their contributions
  // to this block
  for ( int i = 0; i < _currentDescendents.size(); i++ )
  {
    descendent_idx = _currentDescendents[ i ];

    Supernode               &descendent = _factor[ descendent_idx ];

    // Get the index of the interaction the descendent has with "node",
    // (its ancestor)
    ancestor_idx = _nextInteractionCache[ descendent_idx ];

    const DenseBlock &interactionRanges
                                  = _diagonalBlockRanges[ i ][ block_idx ];

    const SupernodeInteraction &interaction
                                  = descendent.offDiagonal()[ ancestor_idx ];

    if ( !interactionRanges.valid() )
    {
      // This descendent does not actually interact
      // with this particular sub-block of node's diagonal
      continue;
    }

    // Determine workspace sizes
    //
    // sub-matrix workspace
    dataSizes[ 0 ] = IndexRange( range_size( interactionRanges._rowRange ),
                                 rank );
    // Initial interaction multiply workspace
    dataSizes[ 1 ] = IndexRange( rank, descendent.numColumns() );
    // Final interaction multiply workspace
    dataSizes[ 2 ] = IndexRange( rank,
                                 range_size( interactionRanges._columnRange ) );

    RealWorkspace          workspace( _realWorkspaceManager, dataSizes );

    Real                  *subMatrix = workspace.workspaceData( 0 );
    Real                  *multWorkspaceInitial = workspace.workspaceData( 1 );
    Real                  *multWorkspaceFinal = workspace.workspaceData( 2 );

#if 0
    descendent.interactionTransMultiply( ancestor_idx,
                                         node, block, Q, rank,
                                         interactionRanges._rowRange,
                                         interactionRanges._columnRange,
                                         _subMatrix,
                                         _copyWorkspace,
                                         _decompMultWorkspace );

    // Scatter to the desired column space
    Supernode::ScatterLowRankTransUpdate( _decompMultWorkspace,
                                          _decompMultTransWorkspace,
                                          interaction,
                                          interactionRanges._columnRange,
                                          rank,
                                          interactionRanges.numColumns(),
                                          block.numColumns(),
                                          // Offset to apply to column indices
                                          block._columnRange.first );
#endif
    descendent.interactionTransMultiply( ancestor_idx,
                                         node, block, Q, rank,
                                         interactionRanges._rowRange,
                                         interactionRanges._columnRange,
                                         subMatrix,
                                         multWorkspaceInitial,
                                         multWorkspaceFinal );

    // Scatter to the desired column space
    Supernode::ScatterLowRankTransUpdate( multWorkspaceFinal,
                                          _decompMultTransWorkspace,
                                          interaction,
                                          interactionRanges._columnRange,
                                          rank,
                                          interactionRanges.numColumns(),
                                          block.numColumns(),
                                          // Offset to apply to column indices
                                          block._columnRange.first );
  }

  // Apply the basis to diagonal contributions resulting from variable
  // compression and add to the result
#if 0
  node.diagonalContributionMult( block_idx, Q, rank,
                                 _slackData, _expansionWorkspace,
                                 _copyWorkspace, _decompMultTransWorkspace,
                                 false, /* Don't transpose the update */
                                 false /* Right side multiplication */ );
#endif
  node.diagonalContributionMult( block_idx, Q, rank,
                                 _slackData, _realWorkspaceManager,
                                 _decompMultTransWorkspace,
                                 false, /* Don't transpose the update */
                                 false /* Right side multiplication */ );

  // Copy V and Utrans data in to the appropriate interaction, allocate
  // data for new node, etc.
  node.assignExtraDiagonalData( block_idx, _factor, Q /* V */,
                                _decompMultTransWorkspace /* Utrans */,
                                rank, _slackData,
                                _slackDataOffset, _availableSlackData );

  TRACE_ASSERT( _availableSlackData >= 0, "Ran out of extra data" );
#endif
}

//////////////////////////////////////////////////////////////////////
// Version of diagonalBasisProject for decomposing transposes of
// lower triangular blocks (ie. decomposing upper triangular blocks)
//////////////////////////////////////////////////////////////////////
void FactorManager::diagonalBasisProject_transpose(
                          const SPARSE_MATRIX::SparseColumnMatrix &A,
                          Supernode &node, int block_idx,
                          bool writeFactors,
                          bool useFactor,
                          int fixedRank )
{
  int                        descendent_idx;
  int                        ancestor_idx;
  int                        rank;

#if 0
  // Workspace sizes
  PairArray                  dataSizes( 3 );
#endif

  Real                      *Q = _decompWorkspace;
  Real                      *U = _decompMultTransWorkspace;

  const DenseBlock &block = node.diagonalLowRankBlocks()[ block_idx ];

  rank = ( fixedRank > 0 ) ? fixedRank : _blockRanks[ 0 ];

  diagonalBlockMultiply( A, node, useFactor, block_idx,
                         rank, NULL, NULL );

  MATRIX                     Utrans( rank, block.numColumns() );

  MATRIX::transposeBLAS( Utrans.data(), U, block.numColumns(), rank );

  // Copy V and Utrans data in to the appropriate interaction, allocate
  // data for new node, etc.
  node.assignExtraDiagonalData( block_idx, _factor, Q /* V */,
                                Utrans.data() /* Utrans */,
                                rank, _slackData,
                                _slackDataOffset, _availableSlackData );

  TRACE_ASSERT( _availableSlackData >= 0, "Ran out of extra data" );
}

//////////////////////////////////////////////////////////////////////
// Performs a forward solve (L x = b) using the computed
// numerical factor
//////////////////////////////////////////////////////////////////////
void FactorManager::forwardSolve( VECTOR &x, int startNode )
{
  TIMING_START( "Forward solve" );

  PairArray                  dataSizes( 1 );

  dataSizes[ 0 ] = IndexPair( x.size(), 1 );

  RealWorkspace              workspace( _realWorkspaceManager, dataSizes );

  for ( int node_idx = startNode; node_idx < _factor.size(); node_idx++ )
  {
    _factor[ node_idx ].forwardSolve( x, _solveWorkspace,
                                      _extendedSolveWorkspace,
                                      _factor, _slackData,
                                      workspace.workspaceData( 0 ) );
                                      //_expansionWorkspace );
  }

  TIMING_STOP( "Forward solve" );
}

//////////////////////////////////////////////////////////////////////
// Solve functions when we are using interior blocks
//////////////////////////////////////////////////////////////////////
void FactorManager::forwardSolve( const SPARSE_MATRIX::SparseColumnMatrix &A,
                                  VECTOR &x )
{
  TIMING_START( "Forward solve" );

  PairArray                  dataSizes( 1 );

  dataSizes[ 0 ] = IndexPair( x.size(), 1 );

  RealWorkspace              workspace( _realWorkspaceManager, dataSizes );

  int                        block_idx;

  for ( int node_idx = 0; node_idx < _factor.size(); node_idx++ )
  {
    block_idx = _interiorBlockMap[ node_idx ];

    if ( block_idx != EMPTY ) {
      // Do a partial solve with this full block
      interiorBlockForwardSolve( A, x, block_idx,
                                 workspace.workspaceData( 0 ) );

      // Skip to the end of the block
      node_idx = _interiorBlocks[ block_idx ]._nodeRange.second;
    } else {
      _factor[ node_idx ].forwardSolve( x, _solveWorkspace,
                                        _extendedSolveWorkspace,
                                        _factor, _slackData,
                                        workspace.workspaceData( 0 ) );
                                        //_expansionWorkspace );
    }
  }

  TIMING_STOP( "Forward solve" );
}

//////////////////////////////////////////////////////////////////////
// Same as the above, but a backward solve (L' x = b)
//////////////////////////////////////////////////////////////////////
void FactorManager::backwardSolve( VECTOR &x, int startNode )
{
  TIMING_START( "Backward solve" );

  PairArray                  dataSizes( 1 );

  dataSizes[ 0 ] = IndexPair( x.size(), 1 );

  RealWorkspace              workspace( _realWorkspaceManager, dataSizes );

  for ( int node_idx = _factor.size() - 1; node_idx >= startNode; node_idx-- )
  {
    _factor[ node_idx ].backwardSolve( x, _solveWorkspace,
                                       _extendedSolveWorkspace,
                                       _factor, _slackData,
                                       workspace.workspaceData( 0 ) );
                                       //_expansionWorkspace );
  }

  TIMING_STOP( "Backward solve" );
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void FactorManager::backwardSolve( const SPARSE_MATRIX::SparseColumnMatrix &A,
                                   VECTOR &x )
{
  TIMING_START( "Backward solve" );

  PairArray                  dataSizes( 1 );

  dataSizes[ 0 ] = IndexPair( x.size(), 1 );

  RealWorkspace              workspace( _realWorkspaceManager, dataSizes );

  int                        block_idx;

  for ( int node_idx = _factor.size() - 1; node_idx >= 0; node_idx-- )
  {
    block_idx = _interiorBlockMap[ node_idx ];

    if ( block_idx != EMPTY ) {
      interiorBlockBackwardSolve( A, x, block_idx,
                                  workspace.workspaceData( 0 ) );

      // Skip to the start of the block
      node_idx = _interiorBlocks[ block_idx ]._nodeRange.first;
    } else {
      _factor[ node_idx ].backwardSolve( x, _solveWorkspace,
                                         _extendedSolveWorkspace,
                                         _factor, _slackData,
                                         workspace.workspaceData( 0 ) );
                                         //_expansionWorkspace );
    }
  }

  TIMING_STOP( "Backward solve" );
}

//////////////////////////////////////////////////////////////////////
// Helper functions for running solves with interior blocks
//////////////////////////////////////////////////////////////////////
void FactorManager::interiorBlockForwardSolve(
                    const SPARSE_MATRIX::SparseColumnMatrix &A,
                    VECTOR &x, int block_idx,
                    Real *workspace )
{
  TIMING_START( "interiorBlockForwardSolve preamble" );

  const InteriorBlock       &block = _interiorBlocks[ block_idx ];
  const IndexRange          &nodeRange = block._nodeRange;
  int                        startCol;
  int                        endCol;
  int                        nCols;

  VECTOR                     extraData;

  Real                      *diagonalData;
  Real                      *diagonalWorkspace;

  startCol = _factor[ nodeRange.first ].startColumn();
  endCol = _factor[ nodeRange.second ].endColumn();
  nCols = endCol - startCol + 1;

  extraData.resizeAndWipe( nCols );

  // Start by performing a forward solve on the desired block
  diagonalData = x.data() + startCol;
  diagonalWorkspace = workspace + startCol;

  TIMING_STOP( "interiorBlockForwardSolve preamble" );

  FLOP_COUNT_START( "interiorBlockForwardSolve forwardSubSolve" );
  TIMING_START( "interiorBlockForwardSolve forwardSubSolve" );
  forwardSubSolve( nodeRange, 1, diagonalData, diagonalWorkspace );
  TIMING_STOP( "interiorBlockForwardSolve forwardSubSolve" );
  FLOP_COUNT_END( "interiorBlockForwardSolve forwardSubSolve" );

  // Next, we have to update remaining indices in the right hand
  // side.  Since we do not solve the factored version of off-diagonal
  // data, an extra solve is needed here
  FLOP_COUNT_START( "interiorBlockForwardSolve backwardSubSolve" );
  TIMING_START( "interiorBlockForwardSolve backwardSubSolve" );
  MATRIX::copy( extraData.data(), diagonalData, nCols, 1 );
  backwardSubSolve( nodeRange, 1, extraData.data(), diagonalWorkspace );
  TIMING_STOP( "interiorBlockForwardSolve backwardSubSolve" );
  FLOP_COUNT_END( "interiorBlockForwardSolve backwardSubSolve" );

#if 0
  // Copy all interaction data to the workspace.  For now, we will just
  // ignore row lists and copy full interactions.
  //
  // TODO: figure out how much this impacts performance
  for ( int interaction_idx = 0; interaction_idx < block._interactions.size();
        interaction_idx++ ) {
    const InteriorBlock::Interaction &interaction
                                  = block._interactions[ interaction_idx ];

    TRACE_ASSERT( interaction.size() > 0 );

    int                      node_idx = interaction[ 0 ].first;
    int                      node_int_idx = interaction[ 0 ].second;
    int                      block_node_idx;
    int                      block_start_row;
    int                      block_columns;

    TRACE_ASSERT( in_range( nodeRange, node_idx ) );

    block_node_idx = _factor[ node_idx ]
                            .offDiagonal()[ node_int_idx ]._nodeID;
    block_start_row = _factor[ block_node_idx ].startColumn();
    block_columns = _factor[ block_node_idx ].numColumns();

    MATRIX::copy( workspace + block_start_row, x.data() + block_start_row,
                  block_columns, 1 );
  }
#endif

  // Multiply by the sparse matrix
  TIMING_START( "interiorBlockForwardSolve subMatrixMultiply" );
#if 0
  SPARSE_MATRIX::subMatrixMultiply( A, extraData.data(),
                                    // Start of the relevant part of the
                                    // rhs vector
                                    x.data() + endCol + 1,
                                    // Starting row in sparse matrix
                                    endCol + 1,
                                    // Starting column in sparse matrix
                                    startCol,
                                    // Number of rows
                                    factorSystemSize() - endCol - 1,
                                    // Number of columns
                                    nCols,
                                    // subtract from the workspace
                                    -1.0,
                                    // Don't clear or transpose
                                    false, false );
#endif
  SPARSE_MATRIX::subMatrixMultiply( A, extraData.data(),
                                    // Start of the relevant part of the
                                    // rhs vector
                                    x.data() + endCol + 1,
                                    // Starting row in sparse matrix
                                    endCol + 1,
                                    // Starting column in the sparse matrix
                                    startCol,
                                    // Row offsets for each column
                                    block._firstOffDiagonalRows,
                                    // Subtract from the workspace
                                    -1.0,
#if 0
                                    // Workspace for vector operations
                                    workspace,
#endif
                                    // Don't clear or transpose
                                    false, false );
  TIMING_STOP( "interiorBlockForwardSolve subMatrixMultiply" );

#if 0
  // Copy solved data back to the right hand side.  For now, we will just
  // ignore row lists and copy full interactions.
  //
  // TODO: figure out how much this impacts performance
  for ( int interaction_idx = 0; interaction_idx < block._interactions.size();
        interaction_idx++ ) {
    const InteriorBlock::Interaction &interaction
                                  = block._interactions[ interaction_idx ];

    TRACE_ASSERT( interaction.size() > 0 );

    int                      node_idx = interaction[ 0 ].first;
    int                      node_int_idx = interaction[ 0 ].second;
    int                      block_node_idx;
    int                      block_start_row;
    int                      block_columns;

    TRACE_ASSERT( in_range( nodeRange, node_idx ) );

    block_node_idx = _factor[ node_idx ]
                            .offDiagonal()[ node_int_idx ]._nodeID;
    block_start_row = _factor[ block_node_idx ].startColumn();
    block_columns = _factor[ block_node_idx ].numColumns();

    MATRIX::copy( x.data() + block_start_row, workspace + block_start_row,
                  block_columns, 1 );
  }
#endif
}

//////////////////////////////////////////////////////////////////////
// Helper function for backward solves with an interior block
//////////////////////////////////////////////////////////////////////
void FactorManager::interiorBlockBackwardSolve(
                    const SPARSE_MATRIX::SparseColumnMatrix &A,
                    VECTOR &x, int block_idx,
                    Real *workspace )
{
  TIMING_START( "interiorBlockBackwardSolve preamble" );

  const InteriorBlock       &block = _interiorBlocks[ block_idx ];
  const IndexRange          &nodeRange = block._nodeRange;
  int                        startCol;
  int                        endCol;
  int                        nCols;

  VECTOR                     extraData;

  Real                      *diagonalData;
  Real                      *diagonalWorkspace;

  startCol = _factor[ nodeRange.first ].startColumn();
  endCol = _factor[ nodeRange.second ].endColumn();
  nCols = endCol - startCol + 1;

  extraData.resizeAndWipe( nCols );

  // Start by performing a forward solve on the desired block
  diagonalData = x.data() + startCol;
  diagonalWorkspace = workspace + startCol;

  TIMING_STOP( "interiorBlockBackwardSolve preamble" );

  // Multiply by the off-diagonal factor contribution, which means
  // we first have to hit this with the sparse matrix for the
  // off-diagonal, then do a solve
  TIMING_START( "interiorBlockBackwardSolve subMatrixMultiply" );
#if 0
  SPARSE_MATRIX::subMatrixMultiply( A,
                                    // Start of the relevant part
                                    // of the rhs vector
                                    x.data() + endCol + 1,
                                    extraData.data(),
                                    // Starting row in the sparse matrix
                                    endCol + 1,
                                    // Starting column in the sparse matrix
                                    startCol,
                                    // Number of rows
                                    factorSystemSize() - endCol - 1,
                                    // Number of columns
                                    nCols,
                                    // Subtract
                                    -1.0,
                                    // Don't clear, and do transpose
                                    false, true );
#endif
  SPARSE_MATRIX::subMatrixMultiply( A,
                                    // Start of the relevant part
                                    // of the rhs vector
                                    x.data() + endCol + 1,
                                    extraData.data(),
                                    // Starting row in the sparse matrix
                                    endCol + 1,
                                    // Starting column in the sparse matrix
                                    startCol,
                                    // Row offsets for each column
                                    block._firstOffDiagonalRows,
                                    // Subtract
                                    -1.0,
#if 0
                                    // Workspace for vector operations
                                    workspace,
#endif
                                    // Don't clear, and do transpose
                                    false, true );
  TIMING_STOP( "interiorBlockBackwardSolve subMatrixMultiply" );

  // Finalize multiplication by doing a forward solve and subtract
  // from the diagonal block
  TIMING_START( "interiorBlockBackwardSolve forwardSubSolve" );
  forwardSubSolve( nodeRange, 1, extraData.data(), diagonalWorkspace );
  MATRIX::axpy( diagonalData, extraData.data(), nCols, 1 );
  TIMING_STOP( "interiorBlockBackwardSolve forwardSubSolve" );

  // Do the final solve
  TIMING_START( "interiorBlockBackwardSolve backwardSubSolve" );
  backwardSubSolve( nodeRange, 1, diagonalData, diagonalWorkspace );
  TIMING_STOP( "interiorBlockBackwardSolve backwardSubSolve" );
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void FactorManager::LDLdiagonalSolve( VECTOR &x, int startNode )
{
  for ( int node_idx = startNode; node_idx < _factor.size(); node_idx++ ) {
    _factor[ node_idx ].LDLdiagonalSolve( x );
  }
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void FactorManager::applyLDLforwardPermutation( VECTOR &x, int startNode )
{
  for ( int node_idx = startNode; node_idx < _factor.size(); node_idx++ ) {
    _factor[ node_idx ].applyLDLforwardPermutation( x );
  }
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void FactorManager::applyLDLbackwardPermutation( VECTOR &x, int startNode )
{
  for ( int node_idx = startNode; node_idx < _factor.size(); node_idx++ ) {
    _factor[ node_idx ].applyLDLbackwardPermutation( x );
  }
}

//////////////////////////////////////////////////////////////////////
// Solves a sub-system specified by the given node list.
// Multiple right hand sides are supported.
//
// This is here to support the "interior block" machinery
//////////////////////////////////////////////////////////////////////
void FactorManager::forwardSubSolve( const IntArray &nodeList, int nRHS,
                                     Real *rhs, Real *workspace )
{
  if ( workspace == NULL ) {
    workspace = _solveWorkspace;
  }

  for ( int node_idx = 0; node_idx < nodeList.size(); node_idx++ ) {
    _factor.at( nodeList[ node_idx ] ).forwardSubSolve(
                                    rhs, nRHS,
                                    _solveWorkspace,
                                    nodeList, _factor,
#ifdef DO_TIMING
                                    &_timers[ INTERIOR_BLOCK_MULT_TRISOLVE ],
                                    &_timers[ INTERIOR_BLOCK_MULT_PROPAGATE ]
#else
                                    NULL, NULL
#endif
                                    );
  }
}

//////////////////////////////////////////////////////////////////////
// With a node range, rather than a list
//////////////////////////////////////////////////////////////////////
void FactorManager::forwardSubSolve( IndexRange nodeRange, int nRHS,
                                     Real *rhs, Real *workspace )
{
  if ( workspace == NULL ) {
    workspace = _solveWorkspace;
  }

  for ( int node_idx = nodeRange.first; node_idx <= nodeRange.second;
        node_idx++ )
  {
    _factor.at( node_idx ).forwardSubSolve(
                                    rhs, nRHS,
                                    _solveWorkspace,
                                    nodeRange, _factor );
  }
}

//////////////////////////////////////////////////////////////////////
// Same as the above, but for backward solves
//////////////////////////////////////////////////////////////////////
void FactorManager::backwardSubSolve( const IntArray &nodeList, int nRHS,
                                      Real *rhs, Real *workspace )
{
  if ( workspace == NULL )
  {
    workspace = _solveWorkspace;
  }

  for ( int node_idx = nodeList.size() - 1; node_idx >= 0; node_idx-- )
  {
    _factor.at( nodeList[ node_idx ] ).backwardSubSolve(
                                    rhs, nRHS,
                                    _solveWorkspace,
                                    nodeList, _factor,
#ifdef DO_TIMING
                                    &_timers[ INTERIOR_BLOCK_MULT_TRISOLVE ],
                                    &_timers[ INTERIOR_BLOCK_MULT_PROPAGATE ]
#else
                                    NULL, NULL
#endif
                                    );
  }
}

//////////////////////////////////////////////////////////////////////
// With a node range rather than a list
//////////////////////////////////////////////////////////////////////
void FactorManager::backwardSubSolve( IndexRange nodeRange, int nRHS,
                                      Real *rhs, Real *workspace )
{
  if ( workspace == NULL )
  {
    workspace = _solveWorkspace;
  }

  for ( int node_idx = nodeRange.second; node_idx >= nodeRange.first;
        node_idx-- )
  {
    _factor.at( node_idx ).backwardSubSolve(
                                    rhs, nRHS,
                                    _solveWorkspace,
                                    nodeRange, _factor );
  }
}

#if 0
//////////////////////////////////////////////////////////////////////
// Returns ranks for a set of low rank blocks
//////////////////////////////////////////////////////////////////////
void FactorManager::blockRanks( const Supernode &node,
                                const IntArray &lowRankBlocks,
                                IntArray &ranks )
{
  int                        nCols = node.numColumns();
  int                        interaction_idx;
  int                        rank;

  ranks.clear();

  for ( int block_idx = 0; block_idx < lowRankBlocks.size(); block_idx++ )
  {
    interaction_idx = lowRankBlocks[ block_idx ];

    const SupernodeInteraction &interaction
      = node.offDiagonal()[ interaction_idx ];

    rank = blockRank( interaction._rowList.size(), nCols );
    rank += _overSampling;
    rank = min( rank, nCols );

    ranks.push_back( blockRank( interaction._rowList.size(), nCols ) );
  }
}
#endif

//////////////////////////////////////////////////////////////////////
// Calculates slack variable memory usage in both the
// off-diagonal part, and the Schur complement for the
// slack variables (prior to elimination).
//////////////////////////////////////////////////////////////////////
void FactorManager::calculateSlackUsage( int maxBlockSize ) const
{
  // Initialize usage arrays for each extended node
  vector<IntArray>           blockUsage( numExtendedNodes() );
  int                        start_node = _factor.size() - blockUsage.size();
  int                        node_idx;
  int                        new_idx;
  int                        rank;
  IntArray                   sizeEstimates( numExtendedNodes(), 0 );
  Real                       offDiagonalMB, schurUsageMB;

  long int                   offDiagonalUsage = 0;
  long int                   schurUsage = 0;

  for ( int ext_idx = 0; ext_idx < blockUsage.size(); ext_idx++ )
  {
    node_idx = start_node + ext_idx;

    const Supernode         &node = _factor[ node_idx ];

    blockUsage[ ext_idx ].resize( node.offDiagonal().size(), 0 );
  }

  // Get an estimate of the size for each node
  for ( int node_idx = 0; node_idx < start_node; node_idx++ )
  {
    const Supernode         &node = _factor[ node_idx ];

    for ( int interaction_idx = 0; interaction_idx < node.offDiagonal().size();
          interaction_idx++ )
    {
      const SupernodeInteraction &interaction
                      = node.offDiagonal()[ interaction_idx ];

      if ( interaction._type == STANDARD_NODE && interaction._active
        || interaction._type == EXTENDED_NODE )
      {
        continue;
      }

      rank = blockRank( interaction._rowList.size(), node.numColumns() );

      if ( interaction._type == STANDARD_NODE ) {
        new_idx = node_idx + interaction._extendedOffset;

        TRACE_ASSERT(
          new_idx
          == node.offDiagonal()[ interaction._extendedInteraction ]._nodeID,
          "Node ID mismatch" );

        new_idx -= start_node;

        // Try a different estimate here
        offDiagonalUsage += rank * ( interaction._rowList.size()
                                     + node.numColumns() );

        sizeEstimates[ new_idx ] = rank;
      }
      else if ( interaction._type == COMPRESSED_NODE ) {
        offDiagonalUsage += rank * ( interaction._rowList.size()
                                     + node.numColumns() );
      }
    }

    if ( !node.inPlaceDiagonal() ) {
      for ( int block_idx = 0; block_idx < node.diagonalLowRankBlocks().size();
            block_idx++ )
      {
        int interaction_idx = node.lowRankDiagonalInteraction( block_idx );
        const DenseBlock &block = node.diagonalLowRankBlocks()[ block_idx ];

        const SupernodeInteraction &interaction
                                      = node.offDiagonal()[ interaction_idx ];

        rank = diagonalBlockRank( block.numRows(), block.numColumns() );

        TRACE_ASSERT( interaction._compressed, "Interaction not compressed" );

        offDiagonalUsage += rank * interaction._compressedColumnList.size();

        new_idx = interaction._nodeID - start_node;

        sizeEstimates[ new_idx ] = rank;
      }
    }
  }

  offDiagonalMB = (Real)offDiagonalUsage * 8.0 / 1024.0 / 1024.0;

  printf( "Compressed off-diagonal usage: %lld non-zeros, %f MB\n",
          (long long int)offDiagonalUsage, offDiagonalMB );

  offDiagonalUsage = 0;

  for ( int node_idx = 0; node_idx < start_node; node_idx++ )
  {
    const Supernode         &node = _factor[ node_idx ];

    for ( int interaction_idx = 0; interaction_idx < node.offDiagonal().size();
          interaction_idx++ )
    {
      const SupernodeInteraction &interaction
                      = node.offDiagonal()[ interaction_idx ];

      if ( interaction._type != EXTENDED_NODE ) {
        continue;
      }

      const Supernode       &ancestor = _factor[ interaction._nodeID ];

      rank = sizeEstimates[ interaction._nodeID - start_node ];

      offDiagonalUsage += rank * node.numColumns();

      // Figure out how this fills in the schur complement
      for ( int next_idx = interaction_idx + 1;
            next_idx < node.offDiagonal().size(); next_idx++ )
      {
        const SupernodeInteraction &nextInteraction
                      = node.offDiagonal()[ next_idx ];

        // Find the corresponding node ID in the ancestor
        for ( int anc_idx = 0; anc_idx < ancestor.offDiagonal().size();
              anc_idx++ )
        {
          if ( ancestor.offDiagonal()[ anc_idx ]._nodeID
            == nextInteraction._nodeID )
          {
            blockUsage[ ancestor.nodeID() - start_node ][ anc_idx ]
              = max( blockUsage[ ancestor.nodeID() - start_node ][ anc_idx ],
                       sizeEstimates[ nextInteraction._nodeID - start_node ]
                     * rank );
          }
        }
      }
    }
  }

  // Total up the usage
  for ( int ext_idx = 0; ext_idx < blockUsage.size(); ext_idx++ )
  {
    schurUsage += sizeEstimates[ ext_idx ] * sizeEstimates[ ext_idx ];

    IntArray &blockSizes = blockUsage[ ext_idx ];

    for ( int interaction_idx = 0; interaction_idx < blockSizes.size();
          interaction_idx++ )
    {
      schurUsage += blockSizes[ interaction_idx ];
    }
  }

  offDiagonalMB = (Real)( offDiagonalUsage * sizeof( Real ) ) / 1024.0 / 1024.0;
  schurUsageMB = (Real)( schurUsage * sizeof( Real ) ) / 1024.0 / 1024.0;

  printf( "Off-diagonal usage: %lld non-zeros, %f MB\n",
          (long long int)offDiagonalUsage, offDiagonalMB );

  printf( "Schur complement usage: %lld non-zeros, %f MB\n",
          (long long int)schurUsage, schurUsageMB );
}

//////////////////////////////////////////////////////////////////////
// Given a set of blocks to be decomposed for this node,
// figure out what the blocks should actually be equal to
// and write the results to disk.
//////////////////////////////////////////////////////////////////////
void FactorManager::writeLowRankBlocks(
                                const SPARSE_MATRIX::SparseColumnMatrix &A,
                                Supernode &node )
{
  int                        nCols = node.numColumns();
  int                        nRows = 0;
  MATRIX                     ident( nCols, nCols );

  for ( int i = 0; i < nCols; i++ )
  {
    ident( i, i ) = 1.0;
  }

  for ( int block_idx = 0; block_idx < _lowRankBlocks.size(); block_idx++ )
  {
    const SupernodeInteraction
            &interaction = node.offDiagonal()[ _lowRankBlocks[ block_idx ] ];
    
    nRows += interaction._rowList.size();
  }

  MATRIX                     input( nCols * _lowRankBlocks.size(), nCols );
  MATRIX                     output( nRows, nCols );
  Real                      *inputData = input.data();

  // Stack a bunch of identity matrices on top of each other
  for ( int block_idx = 0; block_idx < _lowRankBlocks.size(); block_idx++ )
  {
    MATRIX::copy( inputData, ident.data(), nCols, nCols );

    inputData += nCols * nCols;
  }

  // Hijack some existing arrays and just use interactionMultiply
  IntArray                   rankTemp( _blockRanks );
  BoolArray                  activeTemp( _lowRankBlockActive );

  for ( int block_idx = 0; block_idx < _blockRanks.size(); block_idx++ )
  {
    _blockRanks[ block_idx ] = nCols;
    _lowRankBlockActive[ block_idx ] = true;
  }

  Real                      *temp1 = _decompWorkspace;
#if 0
  Real                      *temp2 = _decompMultWorkspace;
#endif
  Real                      *temp3 = _decompMultTransWorkspace;
#if 0
  Real                      *temp4 = _copyWorkspace;
  Real                      *temp5 = _subMatrix;
#endif

#if 0
  MATRIX                     decompMultWorkspace( nRows, nCols );
  MATRIX                     copyWorkspace( nRows, nCols );
  MATRIX                     subMatrix( nRows, nCols );
#endif

  _decompWorkspace = output.data();
  _decompMultTransWorkspace = input.data();

#if 0
  _decompMultWorkspace = decompMultWorkspace.data();
  _copyWorkspace = copyWorkspace.data();
  _subMatrix = subMatrix.data();
#endif

  interactionMultiply( A, node, false /* Only write schur complement*/,
                       NULL, NULL, 0 );

  Real                      *outputData = _decompWorkspace;

  for ( int block_idx = 0; block_idx < _lowRankBlocks.size(); block_idx++ )
  {
    const SupernodeInteraction
            &interaction = node.offDiagonal()[ _lowRankBlocks[ block_idx ] ];
    
    nRows = interaction._rowList.size();

    char    buf[ 1024 ];
    MATRIX  tmp( nRows, nCols, outputData );

    cout << "Writing a block" << endl;

    sprintf( buf, "super_numeric/node_%09d_real_interaction_%02d.matrix",
             node.nodeID(), _lowRankBlocks[ block_idx ] );
    tmp.write( buf );

    outputData += nRows * nCols;
  }

  // Set everything back to normal
  _blockRanks = rankTemp;
  _lowRankBlockActive = activeTemp;
  _decompWorkspace = temp1;
  _decompMultTransWorkspace = temp3;
#if 0
  _decompMultWorkspace = temp2;
  _copyWorkspace = temp4;
  _subMatrix = temp5;
#endif
}


//////////////////////////////////////////////////////////////////////
// Writes the given low rank diagonal block to disk
//////////////////////////////////////////////////////////////////////
void FactorManager::writeLowRankDiagonalBlock(
                                  const SPARSE_MATRIX::SparseColumnMatrix &A,
                                  Supernode &node, int block_idx )
{
  const DenseBlock          &block = node.diagonalLowRankBlocks()[ block_idx ];

  MATRIX                     ident = MATRIX::ident( block.numColumns() );
  MATRIX                     output( block.numRows(), block.numColumns() );

  char                       buf[ 1024 ];

  cout << "Diagonal block multiply" << endl;
  diagonalBlockMultiply( A, node, false /* Only write schur complement */,
                         block_idx, node.numColumns(),
                         ident.data(), output.data() );
  cout << "Done" << endl;

  cout << "Writing the block" << endl;
  sprintf( buf, "node_%09d_diagonal_lr_block_%02d.matrix",
           node.nodeID(), block_idx );
  output.write( buf );
  cout << "Done" << endl;
}

//////////////////////////////////////////////////////////////////////
// Breaks down standard node memory usage in terms of which
// nodes interactions actually act on
//////////////////////////////////////////////////////////////////////
void FactorManager::writeStandardMemoryUsage( int maxBlockSize )
{
  vector<long int>   interactionUsage( _factor.size() - _numExtendedNodes, 0 );
  long int           diagonalUsage;

  VECTOR             interactionUsageMB;

  for ( int node_idx = 0; node_idx < interactionUsage.size(); node_idx++ )
  {
    const Supernode &node = _factor[ node_idx ];

    for ( int interaction_idx = 0; interaction_idx < node.offDiagonal().size();
          interaction_idx++ )
    {
      const SupernodeInteraction
                    &interaction = node.offDiagonal()[ interaction_idx ];

      if ( !interaction._active || interaction._type == EXTENDED_NODE )
      {
        continue;
      }

      const IntArray &rowList = interaction._rowList;

      const Supernode &ancestor = _factor[ interaction._nodeID ];

      if ( ancestor.numColumns() >= maxBlockSize )
      {
        interactionUsage[ ancestor.nodeID() ]
                              += rowList.size() * node.numColumns();
      }
    }
  }

  interactionUsageMB.resizeAndWipe( interactionUsage.size() );

  for ( int node_idx = 0; node_idx < interactionUsage.size(); node_idx++ )
  {
    Real             usageMB;

    usageMB = ( (Real)interactionUsage[ node_idx ] ) * 8.0 / 1024.0 / 1024.0;

    interactionUsageMB( node_idx ) = usageMB;
  }

  interactionUsageMB.write( "usages.vector" );
}

//////////////////////////////////////////////////////////////////////
// Sets the block size for the current decomposition iteration
//////////////////////////////////////////////////////////////////////
void FactorManager::setDecompositionBlockSizes(
                                        const Supernode &node,
                                        int iteration,
                                        bool adaptive,
                                        OffDiagonalCompressionType basisType )
{
  int                        nCols = node.numColumns();

  switch ( basisType ) {
    case COMPRESS_FULL_OFF_DIAGONAL:
    {
      int                    nRows = node.countLowRankRows();

      bool                   hasDescendents = _currentDescendents.size() > 0;

      int                    baseSize = min( nCols, nRows );

      // We should have made one place here to store the basis rank
      TRACE_ASSERT( _blockRanks.size() > 0 );

      if ( adaptive ) {
        _blockRanks[ 0 ] = decompositionBlockSize(
                                      (int)sqrt( (Real)baseSize ), iteration,
                                      hasDescendents );
      } else {
        _blockRanks[ 0 ] = fixedDecompositionBlockSize(
                                      (int)sqrt( (Real)baseSize ),
                                      nRows, nCols, hasDescendents );
      }

      break;
    }
    case COMPRESS_INDIVIDUAL_INTERACTIONS:
    default:
    {
      for ( int block_idx = 0; block_idx < _lowRankBlocks.size(); block_idx++ )
      {
        const SupernodeInteraction &interaction
                            = node.offDiagonal()[ _lowRankBlocks[ block_idx ] ];

        bool hasDescendents = _lowRankDescendents[ block_idx ].size() > 0;

        int baseSize = min( nCols, (int)interaction._rowList.size() );

        // FIXME: disable the non-zero column stuff for now, and just
        // use the number of columns in the matrix
        if ( adaptive )
        {
          _blockRanks[ block_idx ] = decompositionBlockSize(
                                      //interaction._nonZeroColumns, iteration,
                                      (int)sqrt( (Real)baseSize ), iteration,
                                      hasDescendents );
        }
        else
        {
          _blockRanks[ block_idx ] = fixedDecompositionBlockSize(
                                      //interaction._nonZeroColumns, iteration,
                                      (int)sqrt( (Real)baseSize ),
                                      interaction._rowList.size(), nCols,
                                      hasDescendents );
        }

        // Clamp this to be no bigger than the size of the block
        _blockRanks[ block_idx ] = min( _blockRanks[ block_idx ], baseSize );

        // FIXME
        //_blockRanks[ block_idx ] = 16;
      }

      break;
    }
  }
}

//////////////////////////////////////////////////////////////////////
// Size of a block to use during low rank decomposition, given the
// number of non-zero columns in the block, and the current decomposition
// iteration
//////////////////////////////////////////////////////////////////////
int FactorManager::decompositionBlockSize( int nonZeroColumns, int iteration,
                                           bool hasDescendents,
                                           bool isDiagonalBlock )
{
  int                        blockSize;
  Real                       sizeDivisor;

#if 0
  iteration /= 3;

  blockSize = 2 << iteration;
  blockSize = nonZeroColumns / blockSize;
#endif
  if ( isDiagonalBlock )
  {
    if ( hasDescendents )
    {
      sizeDivisor = ceil( pow( 1.5, (Real)iteration ) );
      blockSize = 2 * nonZeroColumns / (int)sizeDivisor;

      blockSize += OVER_SAMPLING;

      blockSize = max( blockSize, 16 );
    }
    else
    {
      blockSize = nonZeroColumns + OVER_SAMPLING;
    }
  }
  else
  {
    if ( hasDescendents )
    {
      //blockSize = 2 << iteration;
      sizeDivisor = ceil( pow( 1.5, (Real)iteration ) );
      //blockSize = 4 * nonZeroColumns / blockSize;
      //blockSize = 4 * nonZeroColumns / (int)sizeDivisor;
      blockSize = 1 * nonZeroColumns / (int)sizeDivisor;

      blockSize += OVER_SAMPLING;

      blockSize = max( blockSize, 16 );
    }
    else
    {
      blockSize = nonZeroColumns + OVER_SAMPLING;
    }
  }

  return blockSize;
}

//////////////////////////////////////////////////////////////////////
// For fixed (non-adaptive) rank decomposition
//////////////////////////////////////////////////////////////////////
int FactorManager::fixedDecompositionBlockSize( int nonZeroColumns,
                                                int nRows, int nCols,
                                                bool hasDescendents,
                                                bool isDiagonalBlock )
{
  int                        blockSize;

#if 0
  if ( isDiagonalBlock )
  {
    if ( hasDescendents )
    {
      blockSize = nonZeroColumns * (int)log2( (Real)nonZeroColumns );
      blockSize /= 6;

      blockSize += OVER_SAMPLING;

      blockSize = max( blockSize, 16 );
    }
    else
    {
      blockSize = nonZeroColumns + OVER_SAMPLING;
    }
  }
  else
  {
    if ( hasDescendents )
    {
      blockSize = nonZeroColumns * (int)log2( (Real)nonZeroColumns );
      blockSize /= 2;

      blockSize += OVER_SAMPLING;

      blockSize = max( blockSize, 16 );
    }
    else
    {
      blockSize = nonZeroColumns + OVER_SAMPLING;
    }
  }
#endif

  if ( isDiagonalBlock ) {
    blockSize = diagonalBlockRank( nRows, nCols );
  } else {
    blockSize = blockRank( nRows, nCols );
  }

  return blockSize;
}

//////////////////////////////////////////////////////////////////////
// Heuristic block rank function
//////////////////////////////////////////////////////////////////////
int FactorManager::blockRank( int nRows, int nCols ) const
{
  // For now, just use the square root of the number of
  // columns as the decomposition rank
#if 0
  return max( 1, (int)ceil( _rankConstant * sqrt( (Real)nCols ) ) );
#endif
  Real                   nonZeroColumns;
  int                    blockSize;

  //nonZeroColumns = sqrt( (Real)min( nCols, nRows ) * _rankConstant );
  nonZeroColumns = sqrt( (Real)min( nCols, nRows ) );

  //blockSize = (int)( nonZeroColumns * log2( nonZeroColumns ) );
  blockSize
    = (int)( nonZeroColumns * log2( nonZeroColumns ) * _rankConstant );
  blockSize += OVER_SAMPLING;

  blockSize = max( blockSize, 16 );

  // Shouldn't be larger than the number of rows or columns
  blockSize = min( blockSize, nRows );
  blockSize = min( blockSize, nCols );

  return blockSize;
}

//////////////////////////////////////////////////////////////////////
// Heuristic block rank function for diagonal blocks
//////////////////////////////////////////////////////////////////////
int FactorManager::diagonalBlockRank( int nRows, int nCols ) const
{
  // For now, just use the square root of the number of
  // columns as the decomposition rank
#if 0
  return max( 1, (int)ceil( _rankConstant * sqrt( (Real)nCols ) ) );
#endif
  Real                   nonZeroColumns;
  int                    blockSize;

  //nonZeroColumns = sqrt( (Real)min( nCols, nRows ) * _rankConstant );
  nonZeroColumns = sqrt( (Real)min( nCols, nRows ) );

  //blockSize = (int)( nonZeroColumns * log2( nonZeroColumns ) );
  blockSize
    = (int)( nonZeroColumns * log2( nonZeroColumns ) * _diagonalRankConstant );
  blockSize += OVER_SAMPLING;

  blockSize = max( blockSize, 16 );

  // Shouldn't be larger than the number of rows or columns
  blockSize = min( blockSize, nRows );
  blockSize = min( blockSize, nCols );

  return blockSize;
}

//////////////////////////////////////////////////////////////////////
// Heuristic maximum rank function
//////////////////////////////////////////////////////////////////////
int FactorManager::maxBlockRank( int nRows, int nCols )
{
#if 0
  int                    blockSize;

  blockSize = max( 1, (int)ceil( 1500.0 * sqrt( (Real)nCols ) ) );
  blockSize = min( blockSize, nRows );
  blockSize = min( blockSize, nCols );

  return blockSize;
  //return nCols;
#endif
  return 5 * blockRank( nRows, nCols );
  //return 1500 * blockRank( nRows, nCols );
}

//////////////////////////////////////////////////////////////////////
// Heuristic maximum rank function for diagonal blocks
//////////////////////////////////////////////////////////////////////
int FactorManager::diagonalMaxBlockRank( int nRows, int nCols )
{
#if 0
  int                    blockSize;

  blockSize = max( 1, (int)ceil( 1500.0 * sqrt( (Real)nCols ) ) );
  blockSize = min( blockSize, nRows );
  blockSize = min( blockSize, nCols );

  return blockSize;
  //return nCols;
#endif
  return 5 * diagonalBlockRank( nRows, nCols );
  //return 1500 * blockRank( nRows, nCols );
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
int FactorManager::fixedBlockRank( int nRows, int nCols )
{
  return max( 1, (int)ceil( 10.0 * sqrt( (Real)nCols ) ) );
}

//////////////////////////////////////////////////////////////////////
// Adaptive low rank decomposition code
//////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////
// Error estimator based on power iteration.
//
// From [Liberty et al. 2007] "Randomized algorithms for..."
// boundMultiplier = 10.0 in this paper.
//////////////////////////////////////////////////////////////////////
void FactorManager::blockErrorEstimate(
                                const SPARSE_MATRIX::SparseColumnMatrix &A,
                                Supernode &node, int powerIterations,
                                Real boundMultiplier,
                                bool transpose, bool useFactor )
{
  int                        vectorSize;
  int                        nCols;
  int                        nBlockRows;
  int                        nColsTotal = 0;
  int                        interaction_idx;
  Real                      *baseData;
  Real                      *baseOutputData;
  Real                      *Q;
  Real                       vecNorm;

  nCols = node.numColumns();

  // Count the number of
  for ( int block_idx = 0; block_idx < _lowRankBlocks.size(); block_idx++ )
  {
    if ( !_lowRankBlockActive[ block_idx ] )
    {
      continue;
    }

    interaction_idx = _lowRankBlocks[ block_idx ];

    nColsTotal += transpose ?
                        node.offDiagonal()[ interaction_idx ]._rowList.size()
                      : nCols;
  }

  //////////////////////////////////////////////////////////////////////
  // Step 1: Estimate the norm of the matrix itself
  //////////////////////////////////////////////////////////////////////

  // Initialize a random vector
  buildRandomTestVectors( node, nColsTotal );

  for ( int iteration = 0; iteration < powerIterations; iteration++ )
  {
    if ( iteration == powerIterations - 1 )
    {
      baseData = _vectorWork1;

      for ( int block_idx = 0; block_idx < _lowRankBlocks.size(); block_idx++ )
      {
        if ( !_lowRankBlockActive[ block_idx ] )
        {
          continue;
        }

        vectorSize = transpose ?
            node.offDiagonal()[ _lowRankBlocks[ block_idx ] ]._rowList.size()
          : nCols;

        vecNorm = VECTOR::norm2( baseData, vectorSize );

        _blockErrors[ block_idx ][ 0 ] = 1.0 / vecNorm;

        baseData += vectorSize;
      }
    }

#if 0
    if ( _useInteriorBlocks )
    {
      if ( transpose ) {
        TRACE_ASSERT( NULL, "Not implemented" );
      }
      else {
        interiorBlockInteractionMultiply( A, node,
                                          _vectorWork1, _vectorWork2, 1 );
      }
    }
    else
    {
#endif
      if ( transpose ) {
        interactionTransMultiply( A, node, useFactor,
                                  _vectorWork1, _vectorWork2, 1 );
      }
      else {
        interactionMultiply( A, node, useFactor,
                             _vectorWork1, _vectorWork2, 1 );
      }
#if 0
    }
#endif

    if ( transpose ) {
      interactionMultiply( A, node, useFactor,
                           _vectorWork2, _vectorWork1, 1 );
    }
    else {
      interactionTransMultiply( A, node, useFactor,
                                _vectorWork2, _vectorWork1, 1 );
    }

    if ( iteration == powerIterations - 1 )
    {
      baseData = _vectorWork1;

      for ( int block_idx = 0; block_idx < _lowRankBlocks.size(); block_idx++ )
      {
        if ( !_lowRankBlockActive[ block_idx ] )
        {
          continue;
        }

        vectorSize = transpose ?
            node.offDiagonal()[ _lowRankBlocks[ block_idx ] ]._rowList.size()
          : nCols;

        vecNorm = VECTOR::norm2( baseData, nCols );

        vecNorm *= _blockErrors[ block_idx ][ 0 ];
        vecNorm = sqrt( vecNorm );

        _blockErrors[ block_idx ][ 0 ] = vecNorm;

        baseData += vectorSize;
      }
    }
  }

  //////////////////////////////////////////////////////////////////////
  // Step 2: Estimate the norm of the matrix less its low rank part
  //         A - Q * Q' * A
  //////////////////////////////////////////////////////////////////////
  buildRandomTestVectors( node, nColsTotal );

  // For each iteration, multiply by
  //    A' * ( I - Q * Q' ) * ( I - Q * Q' ) * A
  for ( int iteration = 0; iteration < powerIterations; iteration++ )
  {
    if ( iteration == powerIterations - 1 )
    {
      baseData = _vectorWork1;

      for ( int block_idx = 0; block_idx < _lowRankBlocks.size(); block_idx++ )
      {
        if ( !_lowRankBlockActive[ block_idx ] )
        {
          continue;
        }

        vectorSize = transpose ?
            node.offDiagonal()[ _lowRankBlocks[ block_idx ] ]._rowList.size()
          : nCols;

        vecNorm = VECTOR::norm2( baseData, vectorSize );

        _blockErrors[ block_idx ][ 1 ] = 1.0 / vecNorm;

        baseData += vectorSize;
      }
    }

    // Multiply by A
#if 0
    if ( _useInteriorBlocks )
    {
      if ( transpose ) {
        TRACE_ASSERT( NULL, "Not implemented" );
      }
      else {
        interiorBlockInteractionMultiply( A, node,
                                          _vectorWork1, _vectorWork2, 1 );
      }
    }
    else
    {
#endif
      if ( transpose ) {
        interactionTransMultiply( A, node, useFactor,
                                  _vectorWork1, _vectorWork2, 1 );
      }
      else {
        interactionMultiply( A, node, useFactor,
                             _vectorWork1, _vectorWork2, 1 );
      }
#if 0
    }
#endif

    baseData = _vectorWork2;
    baseOutputData = _vectorWork1;

    // Multiply each block by ( I - Q * Q' )^2
    for ( int block_idx = 0; block_idx < _lowRankBlocks.size(); block_idx++ )
    {
      if ( !_lowRankBlockActive[ block_idx ] )
      {
        continue;
      }

      interaction_idx = _lowRankBlocks[ block_idx ];

      nBlockRows = transpose ?
                      nCols
                    : node.offDiagonal()[ interaction_idx ]._rowList.size();

      Q = _basisStorage._dataPtrs[ block_idx ]._data;

      for ( int i = 0; i < 2; i++ )
      {
        MATRIX::gemv( Q, baseData, baseOutputData,
                      _basisStorage._dataPtrs[ block_idx ]._nRows, nBlockRows,
                      false /* do not transpose */ );

        MATRIX::gemv( Q, baseOutputData, baseData,
                      _basisStorage._dataPtrs[ block_idx ]._nRows, nBlockRows,
                      true, /* transpose */
                      -1.0, 1.0 /* Subtract */ );
      }

      baseData += nBlockRows;
    }

    // Multiply by A'
    if ( transpose ) {
      interactionMultiply( A, node, useFactor,
                           _vectorWork2, _vectorWork1, 1 );
    }
    else {
      interactionTransMultiply( A, node, useFactor,
                                _vectorWork2, _vectorWork1, 1 );
    }

    if ( iteration == powerIterations - 1 )
    {
      baseData = _vectorWork1;

      for ( int block_idx = 0; block_idx < _lowRankBlocks.size(); block_idx++ )
      {
        if ( !_lowRankBlockActive[ block_idx ] )
        {
          continue;
        }

        vectorSize = transpose ?
            node.offDiagonal()[ _lowRankBlocks[ block_idx ] ]._rowList.size()
          : nCols;

        vecNorm = VECTOR::norm2( baseData, vectorSize );

        vecNorm *= _blockErrors[ block_idx ][ 1 ];
        vecNorm = sqrt( vecNorm );

        // We are treating this as an upper bound.  It is not exact,
        // so we must multiply by the given boundMultiplier, which
        // is given as 10.0 in
        // [Liberty et al. 2007] "Randomized algorithms for the low-rank
        // approximation of matrices"
        _blockErrors[ block_idx ][ 1 ] = vecNorm * boundMultiplier;

        baseData += vectorSize;
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////
// Error estimator used when compressing blocks in a node's diagonal
//
// Note: This function assumes that the current basis for this block
// is stored in _basisStorage._dataPtrs[ 0 ], and that errors should
// be stored in _blockErrors[ 0 ]
//////////////////////////////////////////////////////////////////////
void FactorManager::diagonalBlockErrorEstimate(
                                const SPARSE_MATRIX::SparseColumnMatrix &A,
                                Supernode &node, int block_idx,
                                int powerIterations, Real boundMultiplier,
                                bool transpose, bool useFactor )
{
  int                        nRows;
  int                        nCols;
  Real                      *Q;
  Real                       vecNorm;

  const DenseBlock &block = node.diagonalLowRankBlocks()[ block_idx ];

  //////////////////////////////////////////////////////////////////////
  // Step 1: Estimate the norm of the matrix itself
  //////////////////////////////////////////////////////////////////////

  nRows = transpose ? block.numColumns() : block.numRows();
  nCols = transpose ? block.numRows() : block.numColumns();

  buildRandomTestVectors( node, nCols, true /* single vector */ );

  for ( int iteration = 0; iteration < powerIterations; iteration++ )
  {
    if ( iteration == powerIterations - 1 )
    {
      vecNorm = VECTOR::norm2( _vectorWork1, nCols );

      _blockErrors[ 0 ][ 0 ] = 1.0 / vecNorm;
    }

    if ( transpose ) {
      diagonalBlockTransMultiply( A, node, useFactor, block_idx, 1,
                                  // Override workspaces
                                  _vectorWork2, _vectorWork1 );

      diagonalBlockMultiply( A, node, useFactor, block_idx, 1,
                             // Override workspaces
                             _vectorWork1, _vectorWork2 );
    } else {
      diagonalBlockMultiply( A, node, useFactor, block_idx, 1,
                             // Override workspaces
                             _vectorWork1, _vectorWork2 );

      diagonalBlockTransMultiply( A, node, useFactor, block_idx, 1,
                                  // Override workspaces
                                  _vectorWork2, _vectorWork1 );
    }

    if ( iteration == powerIterations - 1 )
    {
      vecNorm = VECTOR::norm2( _vectorWork1, nCols );

      vecNorm *= _blockErrors[ 0 ][ 0 ];
      vecNorm = sqrt( vecNorm );

      _blockErrors[ 0 ][ 0 ] = vecNorm;
    }
  }

  //////////////////////////////////////////////////////////////////////
  // Step 2: Estimate the norm of the matrix less its low rank part
  //         A - Q * Q' * A
  //////////////////////////////////////////////////////////////////////
  buildRandomTestVectors( node, nCols, 1 );

  // For each iteration, multiply by
  //    A' * ( I - Q * Q' ) * ( I - Q * Q' ) * A
  for ( int iteration = 0; iteration < powerIterations; iteration++ )
  {
    if ( iteration == powerIterations - 1 )
    {
      vecNorm = VECTOR::norm2( _vectorWork1, nCols );

      _blockErrors[ 0 ][ 1 ] = 1.0 / vecNorm;
    }

    // Multiply by A
    if ( transpose ) {
      diagonalBlockTransMultiply( A, node, useFactor, block_idx, 1,
                                  // Override workspaces
                                  _vectorWork1, _vectorWork2 );
    } else {
      diagonalBlockMultiply( A, node, useFactor, block_idx, 1,
                             // Override workspaces
                             _vectorWork1, _vectorWork2 );
    }

#if 0
    baseData = _vectorWork2;
    baseOutputData = _vectorWork1;

    nBlockRows = node.offDiagonal()[ interaction_idx ]._rowList.size();
#endif

    Q = _basisStorage._dataPtrs[ 0 ]._data;

    for ( int i = 0; i < 2; i++ )
    {
      MATRIX::gemv( Q, _vectorWork2, _vectorWork1,
                    _basisStorage._dataPtrs[ 0 ]._nRows, nRows,
                    false /* do not transpose */ );

      MATRIX::gemv( Q, _vectorWork1, _vectorWork2,
                    _basisStorage._dataPtrs[ 0 ]._nRows, nRows,
                    true, /* transpose */
                    -1.0, 1.0 /* Subtract */ );
    }

    // Multiply by A'
    if ( transpose ) {
      diagonalBlockMultiply( A, node, useFactor, block_idx, 1,
                             // Overwrite workspaces
                             _vectorWork2, _vectorWork1 );
    } else {
      diagonalBlockTransMultiply( A, node, useFactor, block_idx, 1,
                                  // Overwrite workspaces
                                  _vectorWork2, _vectorWork1 );
    }

    if ( iteration == powerIterations - 1 )
    {
      vecNorm = VECTOR::norm2( _vectorWork1, nCols );

      vecNorm *= _blockErrors[ 0 ][ 1 ];
      vecNorm = sqrt( vecNorm );

      // We are treating this as an upper bound.  It is not exact,
      // so we must multiply by the given boundMultiplier, which
      // is given as 10.0 in
      // [Liberty et al. 2007] "Randomized algorithms for the low-rank
      // approximation of matrices"
      _blockErrors[ 0 ][ 1 ] = vecNorm * boundMultiplier;
    }
  }
}

//////////////////////////////////////////////////////////////////////
// Checks error for each block currently being decomposed in node.
// Sets the active flag (in _lowRankBlockActive) to false for any
// block A satisfying the error bound
//    || A - Q Q' A ||
//   ------------------ < tolerance
//        || A ||
//
// Also, for those nodes, copy the contents of _decompWorkspace
// to _basisWorkspace to prepare for the next iteration
//////////////////////////////////////////////////////////////////////
bool FactorManager::updateActiveBlocks(
                            const SPARSE_MATRIX::SparseColumnMatrix &A,
                            Supernode &node,
                            int blockSize, Real tolerance,
                            bool transpose,
                            bool useFactor )
{
  int                        interaction_idx;
  int                        nRows;
  int                        nCols = node.numColumns();
  int                        currentBlockSize;
  const Real                *inputData;
  Real                      *outputData = _basisWorkspace;
  bool                       blockActive = false;

  inputData = transpose ? _decompMultTransWorkspace : _decompWorkspace;

  // For each block, get a lower bound on the 2-norm of the block
  // itself, and an upper bound on the norm of the difference
  // between the block and its current low rank decomposition.
#ifdef DO_TIMING
  _timers[ ERROR_ESTIMATE ].tick();
#endif
  blockErrorEstimate( A, node, ERROR_POWER_ITERATIONS, _errorBoundMultiplier,
                      transpose, useFactor );
#ifdef DO_TIMING
  _timers[ ERROR_ESTIMATE ].tock();
#endif

#ifdef DO_TIMING
  _timers[ UPDATE_ACTIVE_BLOCKS ].tick();
#endif
  for ( int block_idx = 0; block_idx < _lowRankBlocks.size(); block_idx++ )
  {
    if ( !_lowRankBlockActive[ block_idx ] )
    {
      continue;
    }

    currentBlockSize = ( blockSize > 0 ) ? blockSize : _blockRanks[ block_idx ];

    interaction_idx = _lowRankBlocks[ block_idx ];

    nRows = transpose ? nCols
                      : node.offDiagonal()[ interaction_idx ]._rowList.size();

    const VEC3F             &blockErrors = _blockErrors[ block_idx ];

    // Check to see if the error tolerance has been met for
    // this block
    // TODO: Try an absolute tolerance here
    if ( blockErrors[ 1 ] / blockErrors[ 0 ] < tolerance
      || _basisStorage._dataPtrs[ block_idx ]._full )
    //if ( blockErrors[ 1 ] < tolerance )
    {
      // Remove this from the active list
      _lowRankBlockActive[ block_idx ] = false;
    }
    else
    {
      // Keep this block, which means we should copy its column
      // data from _decompWorkspace to _basisWorkspace
      MATRIX::copy( outputData, inputData, nRows, currentBlockSize );

      outputData += nRows * currentBlockSize;

      blockActive = true;
    }

    inputData += nRows * currentBlockSize;
  }
#ifdef DO_TIMING
  _timers[ UPDATE_ACTIVE_BLOCKS ].tock();
#endif

  return blockActive;
}

//////////////////////////////////////////////////////////////////////
// Similar to the above, but only acts on a single low rank block
// from the node's diagonal
//////////////////////////////////////////////////////////////////////
bool FactorManager::updateActiveDiagonalBlock(
                                const SPARSE_MATRIX::SparseColumnMatrix &A,
                                Supernode &node, int block_idx,
                                int blockSize, Real tolerance,
                                bool transpose, bool useFactor )
{
  const DenseBlock &block = node.diagonalLowRankBlocks()[ block_idx ];

  int                        nRows = transpose ? block.numColumns()
                                               : block.numRows();

  const VEC3F               &blockErrors = _blockErrors[ 0 ];

#ifdef DO_TIMING
  _timers[ DIAGONAL_ERROR_ESTIMATE ].tick();
#endif
  diagonalBlockErrorEstimate( A, node, block_idx, ERROR_POWER_ITERATIONS,
                              _errorBoundMultiplier,
                              transpose, useFactor );
#ifdef DO_TIMING
  _timers[ DIAGONAL_ERROR_ESTIMATE ].tock();
#endif

#ifdef DO_TIMING
  _timers[ DIAGONAL_UPDATE_ACTIVE_BLOCK ].tick();
#endif

  // Check to see if the error tolerance has been met
  if ( blockErrors[ 1 ] / blockErrors[ 0 ] < tolerance )
  {
    // We're done
#ifdef DO_TIMING
    _timers[ DIAGONAL_UPDATE_ACTIVE_BLOCK ].tock();
#endif
    return false;
  }

  // Keep this block, which means we should copy its column data from
  // _decompWorkspace to _basisWorkspace
  MATRIX::copy( _basisWorkspace,
                // Location of the data depends on whether or not we are
                // decomposing the transpose
                transpose ? _decompMultTransWorkspace : _decompWorkspace,
                nRows, blockSize );

#ifdef DO_TIMING
  _timers[ DIAGONAL_UPDATE_ACTIVE_BLOCK ].tock();
#endif

  return true;
}

//////////////////////////////////////////////////////////////////////
// Takes random Gaussian distributed values from _testMatrix
// and forms a set of normalized vectors (according to column
// sizes in low rank blocks) in _vectorWork1
//////////////////////////////////////////////////////////////////////
void FactorManager::buildRandomTestVectors( Supernode &node, int nColsTotal,
                                            bool singleVector )
{
  Real                      *baseData;
  int                        interaction_idx;
  Real                       vecNorm;
  int                        nCols = node.numColumns();

#if 0
  TRACE_ASSERT( nColsTotal <= _testMatrix.rows() * _testMatrix.cols(),
                "Random matrix is too small" );
#endif

  //////////////////////////////////////////////////////////////////////
  // Step 2: Estimate the norm of the matrix less its low rank part
  //         A - Q * Q' * A
  //////////////////////////////////////////////////////////////////////
  //MATRIX::copy( _vectorWork1, _testMatrix.data(), nColsTotal, 1 );

  // FIXME:
  TRACE_ASSERT( nColsTotal > 0 && nColsTotal <= _vectorWorkSz,
                "Invalid vector column size" );
  MathUtil::randomGaussianMatrix( nColsTotal, 1, _vectorWork1 );

  if ( singleVector )
  {
    vecNorm = VECTOR::norm2( _vectorWork1, nColsTotal );

    MATRIX::scale( _vectorWork1, nColsTotal, 1, 1.0 / vecNorm );
  }
  else
  {
    // Normalize the vector for each block
    baseData = _vectorWork1;

    for ( int block_idx = 0; block_idx < _lowRankBlocks.size(); block_idx++ )
    {
      if ( !_lowRankBlockActive[ block_idx ] )
      {
        continue;
      }

      interaction_idx = _lowRankBlocks[ block_idx ];

      vecNorm = VECTOR::norm2( baseData, nCols );

      MATRIX::scale( baseData, nCols, 1, 1.0 / vecNorm );

      baseData += nCols;
    }
  }
}

//////////////////////////////////////////////////////////////////////
// Assigns workspace data to each of a node's low rank blocks.
// This is used to incrementally build a basis for the column
// space of the block.
//////////////////////////////////////////////////////////////////////
void FactorManager::assignBasisWorkspace( Supernode &node,
                                          bool transpose,
                                          OffDiagonalCompressionType basisType )
{
  int                        remainingSize = _basisStorage._dataSz;
  int                        interaction_idx;
  int                        rank;
  int                        nRows;
  Real                      *baseData = _basisStorage._data;

  long int                   totalSize = 0;

  switch ( basisType ) {
    case COMPRESS_FULL_OFF_DIAGONAL:
    {
      // Only one basis to store
      _basisStorage._dataPtrs.resize( 1 );

      nRows = transpose ? node.numColumns() : node.countLowRankRows();

      rank = maxBlockRank( node.countLowRankRows(), node.numColumns() );

      BlockStorageData      &storageBlock = _basisStorage._dataPtrs[ 0 ];

      // We store the transpose of the basis, due to row-major ordering
      storageBlock._nCols = nRows;
      storageBlock._maxSize = rank * nRows;
      storageBlock._data = baseData;
      storageBlock._nRows = 0;

      totalSize = storageBlock._maxSize;

      break;
    }
    case COMPRESS_INDIVIDUAL_INTERACTIONS:
    default:
    {
      _basisStorage._dataPtrs.resize( _lowRankBlocks.size() );

      for ( int block_idx = 0; block_idx < _lowRankBlocks.size(); block_idx++ )
      {
        interaction_idx = _lowRankBlocks[ block_idx ];
        //rank = _blockRanks[ block_idx ];

        nRows = transpose ? node.numColumns()
                          : node.offDiagonal()[ interaction_idx ]._rowList.size();

        rank = maxBlockRank( node.offDiagonal()[ interaction_idx ]._rowList.size(),
                             node.numColumns() );

        BlockStorageData        &storageBlock = _basisStorage._dataPtrs[ block_idx ];

        // We store the transpose of the basis, due to row-major ordering
        storageBlock._nCols = nRows;
        storageBlock._maxSize = rank * nRows;
        storageBlock._data = baseData;
        storageBlock._nRows = 0;

        baseData += storageBlock._maxSize;

        totalSize += storageBlock._maxSize;
      }

      break;
    }
  }

  if ( totalSize > _basisStorage._dataSz ) {
    printf( "What the fuck is going on?  Node %d of %d\n",
            node.nodeID(), (int)_factor.size() );
    cout << SDUMP( totalSize ) << endl;
  }
  TRACE_ASSERT( totalSize <= _basisStorage._dataSz,
                "Basis storage capacity exceeded" );
}

//////////////////////////////////////////////////////////////////////
// Generates a block of random numbers and multiplies each
// low rank block by this block (possibly using power iteration).
//
// That is, forms the matrix (A * A')^q * A * G for each block to
// be decomposed, where q is the number of power iterations
//
// blockSize: The number of columns in G
//////////////////////////////////////////////////////////////////////
void FactorManager::blockRandomMultiply(
                          const SPARSE_MATRIX::SparseColumnMatrix &A,
                          Supernode &node,
                          int blockSize, int powerIterations,
                          bool transpose, bool useFactor )
{
  int                        nCols = node.numColumns();
  int                        nRows;

#ifdef DO_TIMING
  _timers[ RANDOM_MULTIPLY_GENERATE ].tick();
#endif

  int                        totalSize = 0;

  for ( int block_idx = 0; block_idx < _lowRankBlocks.size(); block_idx++ )
  {
    if ( !_lowRankBlockActive[ block_idx ] )
    {
      continue;
    }

    nRows = node.offDiagonal()[ _lowRankBlocks[ block_idx ] ]._rowList.size();

    if ( blockSize > 0 )
    {
      totalSize += transpose ? ( blockSize * nRows )
                             : ( blockSize * nCols );
    }
    else
    {
      totalSize += transpose ? ( _blockRanks[ block_idx ] * nRows )
                             : ( _blockRanks[ block_idx ] * nCols );
    }
  }

  if ( transpose ) {
    TRACE_ASSERT( totalSize <= _decompWorkspaceSz );

    MathUtil::randomGaussianMatrix( totalSize, 1, _decompWorkspace );
  }
  else {
    TRACE_ASSERT( totalSize <= _decompMultTransWorkspaceSz );

    MathUtil::randomGaussianMatrix( totalSize, 1, _decompMultTransWorkspace );
  }
#ifdef DO_TIMING
  _timers[ RANDOM_MULTIPLY_GENERATE ].tock();
#endif

#if 0
  // Make copies in _decompMultTransWorkspace
  _timers[ RANDOM_MULTIPLY_INIT ].tick();
  initMultiplyWorkspace( node, blockSize );
  _timers[ RANDOM_MULTIPLY_INIT ].tock();
#endif

#ifdef DO_TIMING
  _timers[ RANDOM_MULTIPLY_APPLY ].tick();
#endif
  if ( transpose ) {
    interactionTransMultiply( A, node, useFactor,
                              NULL, NULL, /* Use standard workspaces */
                              blockSize /* Use blockSize for all ranks */ );
  }
  else {
    interactionMultiply( A, node, useFactor,
                         NULL, NULL, /* Use standard workspaces */
                         blockSize /* Use blockSize for all ranks */ );
  }

  // Power iteration, if needed
  for ( int iteration = 0; iteration < powerIterations; iteration++ )
  {
    if ( transpose ) {
      interactionMultiply( A, node, useFactor, NULL, NULL, blockSize );
    }
    else {
      interactionTransMultiply( A, node, useFactor, NULL, NULL, blockSize );
    }

    if ( transpose ) {
      interactionTransMultiply( A, node, useFactor, NULL, NULL, blockSize );
    }
    else {
      interactionMultiply( A, node, useFactor, NULL, NULL, blockSize );
    }
  }
#ifdef DO_TIMING
  _timers[ RANDOM_MULTIPLY_APPLY ].tock();
#endif
}

//////////////////////////////////////////////////////////////////////
// Similar to the above, but used to decompose the full off-diagonal
// block of a node all at once
//////////////////////////////////////////////////////////////////////
void FactorManager::blockRandomMultiply_allRows(
                          const SPARSE_MATRIX::SparseColumnMatrix &A,
                          Supernode &node,
                          int blockSize, int powerIterations,
                          bool transpose, bool useFactor )
{
  int                        nCols = node.numColumns();
  int                        nRows = 0;
  int                        totalSize;

#ifdef DO_TIMING
  _timers[ RANDOM_MULTIPLY_GENERATE ].tick();
#endif

  nRows = node.countLowRankRows();

  TRACE_ASSERT( blockSize > 0 || _blockRanks.size() > 0 );
  if ( blockSize > 0 ) {
    totalSize = transpose ? ( blockSize * nRows )
                          : ( blockSize * nCols );
  } else {
    totalSize = transpose ? ( _blockRanks[ 0 ] * nRows )
                          : ( _blockRanks[ 0 ] * nCols );
  }

  if ( transpose ) {
    TRACE_ASSERT( totalSize <= _decompWorkspaceSz );

    MathUtil::randomGaussianMatrix( totalSize, 1, _decompWorkspace );

#if 0
    // FIXME: debugging
    MATRIX::write( _decompWorkspace, nRows,
                   ( blockSize > 0 ) ? blockSize : _blockRanks[ 0 ],
                   "multInput.matrix" );
#endif
  } else {
    TRACE_ASSERT( totalSize <= _decompMultTransWorkspaceSz );

    MathUtil::randomGaussianMatrix( totalSize, 1, _decompMultTransWorkspace );
  }
#ifdef DO_TIMING
  _timers[ RANDOM_MULTIPLY_GENERATE ].tock();
#endif

#ifdef DO_TIMING
  _timers[ RANDOM_MULTIPLY_APPLY ].tick();
#endif
  if ( transpose ) {
    interactionTransMultiply_allRows( A, node, useFactor,
                                      NULL, NULL, /* Use standard workspace */
                                      blockSize /* Use block size for rank */);

#if 0
    // FIXME: debugging
    MATRIX::write( _decompMultTransWorkspace, nCols,
                   ( blockSize > 0 ) ? blockSize : _blockRanks[ 0 ],
                   "multOutput1.matrix" );
    printf( "Done multiplying with node %d\n", node.nodeID() );
    abort();
#endif
  } else {
    interactionMultiply_allRows( A, node, useFactor,
                                 NULL, NULL, /* Use standard workspace */
                                 blockSize /* Use block size for rank */ );
  }

  // Power iteration, if needed
  for ( int iteration = 0; iteration < powerIterations; iteration++ ) {
    if ( transpose ) {
      interactionMultiply_allRows( A, node, useFactor,
                                   NULL, NULL, blockSize );

#if 0
      // FIXME: debugging
      MATRIX::write( _decompWorkspace, nRows,
                     ( blockSize > 0 ) ? blockSize : _blockRanks[ 0 ],
                     "multOutput2.matrix" );
#endif
    }
    else {
      interactionTransMultiply_allRows( A, node, useFactor,
                                        NULL, NULL, blockSize );
    }

    if ( transpose ) {
      interactionTransMultiply_allRows( A, node, useFactor,
                                        NULL, NULL, blockSize );

#if 0
      // FIXME: debugging
      MATRIX::write( _decompMultTransWorkspace, nCols,
                     ( blockSize > 0 ) ? blockSize : _blockRanks[ 0 ],
                     "multOutput3.matrix" );
#endif
    }
    else {
      interactionMultiply_allRows( A, node, useFactor,
                                   NULL, NULL, blockSize );
    }
  }
#ifdef DO_TIMING
  _timers[ RANDOM_MULTIPLY_APPLY ].tock();
#endif
}

//////////////////////////////////////////////////////////////////////
// Generates a block of random numbers and multiplies a single low
// rank block from the given node's diagonal by that matrix
//
// ie. forms the matrix (A * A')^q * A * G for the block to be
// decomposed, where q is the number of power iterations.
//
// blockSize: The number of columns in G
//////////////////////////////////////////////////////////////////////
void FactorManager::diagonalBlockRandomMultiply(
                          const SPARSE_MATRIX::SparseColumnMatrix &A,
                          Supernode &node, int block_idx,
                          int blockSize, int powerIterations,
                          bool transpose, bool useFactor )
{
  int                        nCols;

  const DenseBlock &block = node.diagonalLowRankBlocks()[ block_idx ];

  nCols = transpose ? block.numRows() : block.numColumns();

  // No need to put this in test matrix - just put it directly where
  // we need it to be
  if ( transpose ) {
    TRACE_ASSERT( nCols * blockSize <= _decompWorkspaceSz );

    MathUtil::randomGaussianMatrix( blockSize, nCols, _decompWorkspace );
  } else {
    TRACE_ASSERT( nCols * blockSize <= _decompMultTransWorkspaceSz );

    MathUtil::randomGaussianMatrix( nCols, blockSize,
                                    _decompMultTransWorkspace );
  }

  if ( transpose ) {
    diagonalBlockTransMultiply( A, node, useFactor,
                                block_idx, blockSize, NULL, NULL );
  } else {
    diagonalBlockMultiply( A, node, useFactor,
                           block_idx, blockSize, NULL, NULL );
  }

  // Power iteration, if needed
  for ( int iteration = 0; iteration < powerIterations; iteration++ )
  {
    if ( transpose ) {
      diagonalBlockMultiply( A, node, useFactor,
                             block_idx, blockSize, NULL, NULL );
      diagonalBlockTransMultiply( A, node, useFactor,
                                  block_idx, blockSize, NULL, NULL );
    } else {
      diagonalBlockTransMultiply( A, node, useFactor,
                                  block_idx, blockSize, NULL, NULL );
      diagonalBlockMultiply( A, node, useFactor,
                             block_idx, blockSize, NULL, NULL );
    }
  }
}

//////////////////////////////////////////////////////////////////////
// Counts the number of active low rank rows (that is, rows that
// are part of blocks that still need columns added to their
// bases).
//////////////////////////////////////////////////////////////////////
int FactorManager::countActiveLowRankRows( Supernode &node )
{
  int                        nRowsTotal = 0;
  int                        interaction_idx;

  for ( int block_idx = 0; block_idx < _lowRankBlocks.size(); block_idx++ )
  {
    if ( !_lowRankBlockActive[ block_idx ] )
    {
      continue;
    }

    interaction_idx = _lowRankBlocks[ block_idx ];

    nRowsTotal += (int)node.offDiagonal()[ interaction_idx ]._rowList.size();
  }

  return nRowsTotal;
}

//////////////////////////////////////////////////////////////////////
// Initialization step for adaptive low rank decomposition
//////////////////////////////////////////////////////////////////////
void FactorManager::initOffDiagonalDecomposition(
                            const SPARSE_MATRIX::SparseColumnMatrix &A,
                            Supernode &node,
                            int powerIterations,
                            bool adaptive,
                            bool transpose,
                            bool useFactor,
                            OffDiagonalCompressionType basisType )
{
  int                        currentBlockSize;
  Real                      *baseOutputData = _basisWorkspace;
  const Real                *baseInputData;
  
  baseInputData = transpose ? _decompMultTransWorkspace
                            : _decompWorkspace;

  setDecompositionBlockSizes( node, 0 /* iteration 0 */, adaptive, basisType );

  // Copy to block rank cache
  _blockRanksCache.clear();
  for ( int block_idx = 0; block_idx < _lowRankBlocks.size(); block_idx++ )
  {
    _blockRanksCache.push_back( _blockRanks[ block_idx ] );
  }

  assignBasisWorkspace( node, transpose, basisType );

  switch ( basisType ) {
    case COMPRESS_FULL_OFF_DIAGONAL: {
      _blockErrors.resize( _lowRankBlocks.size() );
      _blockErrors[ 0 ][ 0 ] = 1.0;
      _blockErrors[ 0 ][ 1 ] = 1.0;
      // Ignore the rest

      break;
    }
    case COMPRESS_INDIVIDUAL_INTERACTIONS:
    default: {
      _blockErrors.resize( _lowRankBlocks.size() );

      for ( int block_idx = 0; block_idx < _lowRankBlocks.size(); block_idx++ ) {
        _blockErrors[ block_idx ][ 0 ] = 1.0;
        _blockErrors[ block_idx ][ 1 ] = 1.0;
      }

      break;
    }
  }

  // Multiply each block (with the desired number of power iterations)
  // by a random Gaussian matrix, and put the results in _basisWorkspace
  if ( basisType == COMPRESS_INDIVIDUAL_INTERACTIONS ) {
    blockRandomMultiply( A, node, -1, /* Use block size in _blockRanks */
                         powerIterations, transpose, useFactor );
  } else {
    blockRandomMultiply_allRows(
                         A, node, -1, /* Use block size in _blockRanks */
                         powerIterations, transpose, useFactor );
  }

#if 0
#ifdef DO_TIMING
  _timers[ BASIS_COPY ].tick();
#endif
  long int copySize = 0;
  int blockRows;
  for ( int block_idx = 0; block_idx < _lowRankBlocks.size(); block_idx++ )
  {
    const SupernodeInteraction &interaction
                          = node.offDiagonal()[ _lowRankBlocks[ block_idx ] ];

    currentBlockSize = _blockRanks[ block_idx ];

    blockRows = transpose ? node.numColumns()
                          : interaction._rowList.size();

    MATRIX::copy( baseOutputData, baseInputData, blockRows, currentBlockSize );

    baseInputData += blockRows * currentBlockSize;
    baseOutputData += blockRows * currentBlockSize;

    copySize += blockRows * currentBlockSize;
  }
  TRACE_ASSERT( copySize <= _basisWorkspaceSz,
                "Basis workspace capacity exceeded" );
#ifdef DO_TIMING
  _timers[ BASIS_COPY ].tock();
#endif
#endif
  copyAllBases( node, baseInputData, baseOutputData, transpose, basisType );

#if 0
    // FIXME: debugging 
    MATRIX::write( _basisWorkspace, node.numColumns(), _blockRanks[ 0 ],
                   "basisWorkspace.matrix" );
#endif

  // If we aren't doing adaptive decomposition then we're done here (ie.
  // there's no need to set up a new block)
  if ( !adaptive ) {
    return;
  }

  // Set up another premultiplied random block
  if ( basisType == COMPRESS_INDIVIDUAL_INTERACTIONS ) {
    blockRandomMultiply( A, node, -1, /* Use block size in _blockRanks */
                         powerIterations, transpose, useFactor );
  } else {
    blockRandomMultiply_allRows(
                         A, node, -1, /* Use block size in _blockRanks */
                         powerIterations, transpose, useFactor );
  }
}

//////////////////////////////////////////////////////////////////////
// Helper function for initOffDiagonalDecomposition
//////////////////////////////////////////////////////////////////////
void FactorManager::copyAllBases( const Supernode &node,
                                  const Real *baseInputData,
                                  Real *baseOutputData,
                                  bool transpose,
                                  OffDiagonalCompressionType basisType )
{
  long int                   copySize = 0;
  int                        blockRows;
  int                        currentBlockSize;

#ifdef DO_TIMING
  _timers[ BASIS_COPY ].tick();
#endif

  switch ( basisType ) {
    case COMPRESS_FULL_OFF_DIAGONAL: {
      currentBlockSize = _blockRanks[ 0 ];

      blockRows = transpose ? node.numColumns() : node.countLowRankRows();

      MATRIX::copy( baseOutputData, baseInputData, blockRows, currentBlockSize );

      copySize = blockRows * currentBlockSize;

      break;
    }
    case COMPRESS_INDIVIDUAL_INTERACTIONS: {
      for ( int block_idx = 0; block_idx < _lowRankBlocks.size(); block_idx++ )
      {
        const SupernodeInteraction &interaction
                              = node.offDiagonal()[ _lowRankBlocks[ block_idx ] ];

        currentBlockSize = _blockRanks[ block_idx ];

        blockRows = transpose ? node.numColumns()
                              : interaction._rowList.size();

        MATRIX::copy( baseOutputData, baseInputData, blockRows, currentBlockSize );

        baseInputData += blockRows * currentBlockSize;
        baseOutputData += blockRows * currentBlockSize;

        copySize += blockRows * currentBlockSize;
      }

      break;
    }
  }

  TRACE_ASSERT( copySize <= _basisWorkspaceSz,
                "Basis workspace capacity exceeded" );

#ifdef DO_TIMING
  _timers[ BASIS_COPY ].tock();
#endif
}

//////////////////////////////////////////////////////////////////////
// Initialization for adaptive decomposition of a diagonal low rank block
//////////////////////////////////////////////////////////////////////
void FactorManager::initDiagonalDecomposition(
                                const SPARSE_MATRIX::SparseColumnMatrix &A,
                                Supernode &node, int block_idx,
                                int powerIterations,
                                int &blockSize,
                                bool adaptive,
                                bool transpose,
                                bool useFactor )
{
  const int                  SINGLE_BLOCK = 1;

  int                        maxRank;
  int                        nRows;
  int                        nonZeroColumns;

  const DenseBlock &block = node.diagonalLowRankBlocks()[ block_idx ];

  // FIXME: just use the size of the block for now
  //nonZeroColumns = node.diagonalLowRankBlockNonZeroColumns( block_idx );
  nonZeroColumns
    = (int)sqrt(
        3.0 * (Real)node.diagonalLowRankBlocks()[ block_idx ].numColumns() );

  // FIXME: Set the block size in some intelligent way
  if ( adaptive )
  {
    blockSize = decompositionBlockSize( nonZeroColumns, 0, /* iteration 0 */
                                        _currentDescendents.size() > 0,
                                        true /* diagonal block */ );
  }
  else
  {
    blockSize = fixedDecompositionBlockSize( nonZeroColumns,
                                             block.numRows(), block.numColumns(),
                                             _currentDescendents.size() > 0,
                                             true /* diagonal block */ );
  }

  nRows = transpose ? block.numColumns() : block.numRows();

  // Set up basis workspace
  _basisStorage._dataPtrs.resize( SINGLE_BLOCK );

  maxRank = diagonalMaxBlockRank( block.numRows(), block.numColumns() );

  BlockStorageData          &storageBlock = _basisStorage._dataPtrs[ 0 ];

  // We store the transpose of the basis
  storageBlock._nCols = nRows;
  storageBlock._maxSize = maxRank * nRows;
  storageBlock._data = _basisStorage._data;
  storageBlock._nRows = 0;

  // Error initialization
  _blockErrors.resize( SINGLE_BLOCK );

  _blockErrors[ 0 ][ 0 ] = 1.0;
  _blockErrors[ 0 ][ 1 ] = 1.0;

  _blockRanks.resize( SINGLE_BLOCK );
  _blockRanks[ 0 ] = maxRank;

  // Multiply the block by a random Gaussian matrix and store the
  // output in _basisWorkspace
#ifdef DO_TIMING
  _timers[ DIAGONAL_RANDOM_MULTIPLY ].tick();
#endif
  diagonalBlockRandomMultiply( A, node, block_idx, blockSize, powerIterations,
                               transpose, useFactor );
#ifdef DO_TIMING
  _timers[ DIAGONAL_RANDOM_MULTIPLY ].tock();
#endif

#ifdef DO_TIMING
  _timers[ DIAGONAL_BASIS_COPY ].tick();
#endif
  MATRIX::copy( _basisWorkspace,
                // Location of multiply output depends on whether or not
                // we are decomposing the transpose
                transpose ? _decompMultTransWorkspace : _decompWorkspace,
                nRows, blockSize );
#ifdef DO_TIMING
  _timers[ DIAGONAL_BASIS_COPY ].tock();
#endif

  // No need to set up another block if we are using fixed-rank decomposition
  if ( !adaptive )
  {
    return;
  }

  // Set up another premultiplied random block
#ifdef DO_TIMING
  _timers[ DIAGONAL_RANDOM_MULTIPLY ].tick();
#endif
  diagonalBlockRandomMultiply( A, node, block_idx, blockSize, powerIterations,
                               transpose, useFactor );
#ifdef DO_TIMING
  _timers[ DIAGONAL_RANDOM_MULTIPLY ].tock();
#endif
}

//////////////////////////////////////////////////////////////////////
// Adaptive decomposition iteration
//
// Roughly follows algorithm 4.2 of [Halko et al. 2009], "Finding Structure
// with Randomness: Stochastic Algorithms for Constructing Approximate
// Matrix Decompositions"
//////////////////////////////////////////////////////////////////////
void FactorManager::offDiagonalDecomposition(
                            const SPARSE_MATRIX::SparseColumnMatrix &A,
                            Supernode &node,
                            int powerIterations,
                            Real tolerance,
                            bool adaptive,
                            bool transpose,
                            bool useFactor,
                            OffDiagonalCompressionType basisType )
{
  initOffDiagonalDecomposition( A, node, powerIterations, adaptive,
                                transpose, useFactor,
                                basisType );

  if ( adaptive ) {
    adaptiveOffDiagonalDecomposition( A, node, powerIterations, tolerance,
                                      transpose, useFactor,
                                      basisType );
  } else {
    fixedOffDiagonalDecomposition( A, node, powerIterations,
                                   transpose, useFactor,
                                   basisType );
  }

  // Get everything in to _decompWorkspace.
  //
  // This requires an additional copy/transpose, but hopefully
  // that isn't too big of a deal.
#ifdef DO_TIMING
  _timers[ BASIS_COPY ].tick();
#endif
  if ( basisType == COMPRESS_FULL_OFF_DIAGONAL ) {
    // Only one basis to copy
    copyFinalBases( true, transpose );
  } else {
    copyFinalBases( false, transpose );
  }
#ifdef DO_TIMING
  _timers[ BASIS_COPY ].tock();
#endif
}

//////////////////////////////////////////////////////////////////////
// Fixed-rank off-diagonal decomposition
//////////////////////////////////////////////////////////////////////
void FactorManager::fixedOffDiagonalDecomposition(
                              const SPARSE_MATRIX::SparseColumnMatrix &A,
                              Supernode &node,
                              int powerIterations, bool transpose,
                              bool /* useFactor */,
                              OffDiagonalCompressionType basisType )
{
  // Compute a QR factorization for each interaction
#ifdef DO_TIMING
  _timers[ QR_FACTORIZATION ].tick();
#endif
  switch ( basisType ) {
    case COMPRESS_FULL_OFF_DIAGONAL: {
      node.lowRankQR( _basisWorkspace, _blockRanks[ 0 ],
                      _qrExtraData, _qrWorkspace, _qrExtraDataSz,
                      transpose );

      break;
    }
    case COMPRESS_INDIVIDUAL_INTERACTIONS:
    default: {
      node.lowRankQR( _lowRankBlocks, _lowRankBlockActive,
                      _basisWorkspace, _blockRanks,
                      _qrExtraData, _qrWorkspace, _qrExtraDataSz,
                      transpose );

      break;
    }
  }
#ifdef DO_TIMING
  _timers[ QR_FACTORIZATION ].tock();
#endif

  // Copy the new basis blocks in to the existing bases
#ifdef DO_TIMING
  _timers[ BASIS_COPY ].tick();
#endif
  switch ( basisType ) {
    case COMPRESS_FULL_OFF_DIAGONAL: {
      appendOffDiagonalBasis( node, -1 /* Use block size from _blockRanks */,
                              0 /* iteration 0 */, &A, transpose );

#if 0
      // FIXME: debugging
      MATRIX::write( _basisStorage._dataPtrs[ 0 ]._data,
                     _blockRanks[ 0 ], node.numColumns(),
                     "postStorageCopy.matrix" );
#endif

      break;
    }
    case COMPRESS_INDIVIDUAL_INTERACTIONS:
    default: {
      appendOffDiagonalBases( node, -1 /* Use block size from _blockRanks */,
                              0 /* iteration 0 */, &A, transpose );

      break;
    }
  }
#ifdef DO_TIMING
  _timers[ BASIS_COPY ].tock();
#endif
}

//////////////////////////////////////////////////////////////////////
// Adaptive off diagonal decomposition
//////////////////////////////////////////////////////////////////////
void FactorManager::adaptiveOffDiagonalDecomposition(
                              const SPARSE_MATRIX::SparseColumnMatrix &A,
                              Supernode &node,
                              int powerIterations,
                              Real tolerance,
                              bool transpose,
                              bool useFactor,
                              OffDiagonalCompressionType basisType )
{
  int                        iteration = 0;

  TRACE_ASSERT( basisType == COMPRESS_INDIVIDUAL_INTERACTIONS,
                "Adaptive decomposition of full off-diagonal not supported" );

  for ( ; ; )
  {
    iteration += 1;

    // Project the current contents of the basis workspace to
    // the complement space of each block's basis
#ifdef DO_TIMING
    _timers[ ORTHOGONAL_PROJECTION ].tick();
#endif
    basisOrthogonalProjection( node, -1 /* Use block size from _blockRanks */,
                               transpose );
#ifdef DO_TIMING
    _timers[ ORTHOGONAL_PROJECTION ].tock();
#endif

    // Compute a QR factorization for each interaction
#ifdef DO_TIMING
    _timers[ QR_FACTORIZATION ].tick();
#endif
    node.lowRankQR( _lowRankBlocks, _lowRankBlockActive,
                    _basisWorkspace, _blockRanks,
                    _qrExtraData, _qrWorkspace, _qrExtraDataSz, transpose );
#ifdef DO_TIMING
    _timers[ QR_FACTORIZATION ].tock();
#endif

    // Project this basis in to the complement space of the
    // existing basis.
#ifdef DO_TIMING
    _timers[ ORTHOGONAL_PROJECTION ].tick();
#endif
    basisOrthogonalProjection( node, -1 /* Use block size from _blockRanks */,
                               transpose );
#ifdef DO_TIMING
    _timers[ ORTHOGONAL_PROJECTION ].tock();
#endif

    // Copy the new basis blocks in to the existing bases
#ifdef DO_TIMING
    _timers[ BASIS_COPY ].tick();
#endif
    appendOffDiagonalBases( node, -1 /* Use block size from _blockRanks */,
                            iteration, &A, transpose );
#ifdef DO_TIMING
    _timers[ BASIS_COPY ].tock();
#endif

    cacheBlockRanks();

    // Check errors and break out of the loop if there are no blocks
    // remaining.  Otherwise copy data to _basisWorkspace to prepare
    // for the next iteration.
    if ( !updateActiveBlocks( A, node, -1, /* Use block size from _blockRanks */
                              tolerance, transpose, useFactor ) )
    {
      break;
    }

    cacheBlockRanks();

    setDecompositionBlockSizes( node, iteration, true, basisType );

    // Build a new block column of pre-multiplied random Gaussian data
    blockRandomMultiply( A, node, -1, /* Use block size from _blockRanks */
                         powerIterations, transpose, useFactor );

    cacheBlockRanks();
  }
}

//////////////////////////////////////////////////////////////////////
// Adaptive decomposition for low rank blocks on the main diagonal of
// a node
//////////////////////////////////////////////////////////////////////
void FactorManager::diagonalDecomposition(
                            const SPARSE_MATRIX::SparseColumnMatrix &A,
                            Supernode &node,
                            int powerIterations,
                            Real tolerance,
                            bool adaptive,
                            bool transpose,
                            bool useFactor )
{
  int                        numBlocks = node.diagonalLowRankBlocks().size();

  if ( numBlocks == 0 )
  {
    return;
  }

  // Set up diagonal block ranges for system multiplication
#ifdef DO_TIMING
  _timers[ FIND_BLOCK_RANGES ].tick();
#endif
  findDiagonalBlockRanges( node );
#ifdef DO_TIMING
  _timers[ FIND_BLOCK_RANGES ].tock();
#endif

  for ( int block_idx = 0; block_idx < numBlocks; block_idx++ )
  {
    diagonalDecomposition( A, node, block_idx, powerIterations, tolerance,
                           adaptive, transpose, useFactor );
  }
}

//////////////////////////////////////////////////////////////////////
// Factorization of a node's diagonal with in-place block compression
//////////////////////////////////////////////////////////////////////
SolverErrorPtr FactorManager::factorCompressedDiagonal(
                              const SPARSE_MATRIX::SparseColumnMatrix &A,
                              Supernode &node,
                              int powerIterations,
                              Real tolerance,
                              bool adaptive,
                              bool transpose,
                              bool useFactor )
{
  int                        numBlocks = node.diagonalLowRankBlocks().size();

  DiagonalBlock             &nodeDiagonal = node.nodeDiagonal();

  SolverErrorPtr             error;

  if ( numBlocks > 0 )
  {
    // Set up diagonal block ranges for system multiplication
#ifdef DO_TIMING
    _timers[ FIND_BLOCK_RANGES ].tick();
#endif
    findDiagonalBlockRanges( node );
#ifdef DO_TIMING
    _timers[ FIND_BLOCK_RANGES ].tock();
#endif
  }

  TRACE_ASSERT( useFactor, "useFactor == false not supported" );

#if 0
  // FIXME: debugging
  if ( nodeDiagonal.numBlocks() > 1 ) {
    printf( "Node %d's diagonal has %d blocks over %d columns\n",
            node.nodeID(), nodeDiagonal.numBlocks(), node.numColumns() );
  }
#endif

  // The DiagonalBlock data structure determines the correct order in
  // which to handle diagonal and off-diagonal blocks in the diagonal
  for ( int block_idx = 0; block_idx < nodeDiagonal.numBlocks(); block_idx++ )
  {
    DiagonalBlock::BlockEntry  blockEntry
                                    = nodeDiagonal.blockEntry( block_idx );

    if ( blockEntry._buildType == DiagonalBlock::EXPLICIT_BUILD ) {
      // Form this entire block explicitly, then compress its subblocks
      error = explicitDiagonalDecomposition( A, node, block_idx,
                                             powerIterations,
                                             transpose, useFactor );

      if ( error != SolverErrorPtr() ) {
        // Something went wrong - pass it to the caller
        return error;
      }
    } else if ( blockEntry._isDiagonal ) {
      TIMING_START( "factorCompressedDiagonal: update diagonal" );
      nodeDiagonal.updateDiagonal( blockEntry._index );
      TIMING_STOP( "factorCompressedDiagonal: update diagonal" );
      TIMING_START( "factorCompressedDiagonal: factor diagonal" );
      error = nodeDiagonal.factorDiagonal( blockEntry._index );
      TIMING_STOP( "factorCompressedDiagonal: factor diagonal" );

      if ( error != SolverErrorPtr() ) {
        // Something went wrong - pass it to the caller
        return error;
      }
    } else {
      // Form a low-rank decomposition for this block
      TIMING_START( "factorCompressedDiagonal: diagonalDecomposition" );
      diagonalDecomposition( A, node, blockEntry._index, powerIterations,
                             tolerance, adaptive, transpose, useFactor );
      TIMING_STOP( "factorCompressedDiagonal: diagonalDecomposition" );
    }
  }

#if 0
  // FIXME: debugging
  if ( nodeDiagonal.numDiagonalBlocks() > 1 ) {
    char buf[ 1024 ];
    sprintf( buf, "super_numeric/node_%05d_", node.nodeID() );
    nodeDiagonal.writeBlocks( buf );
  }
#endif

  return SolverErrorPtr();
}

//////////////////////////////////////////////////////////////////////
// Adaptive decomposition for a single low rank block in the node's
// main diagonal
//////////////////////////////////////////////////////////////////////
void FactorManager::diagonalDecomposition(
                            const SPARSE_MATRIX::SparseColumnMatrix &A,
                            Supernode &node,
                            int block_idx,
                            int powerIterations,
                            Real tolerance,
                            bool adaptive,
                            bool transpose,
                            bool useFactor )
{
  int                        blockSize;

  initDiagonalDecomposition( A, node, block_idx, powerIterations, blockSize,
                             adaptive, transpose, useFactor );

  if ( adaptive )
  {
    adaptiveDiagonalDecomposition( A, node, block_idx, powerIterations,
                                   tolerance, blockSize,
                                   transpose, useFactor );
  }
  else
  {
    fixedDiagonalDecomposition( A, node, block_idx, powerIterations,
                                blockSize,
                                transpose, useFactor );
  }

  // Copy the basis to _decompWorkspace
#ifdef DO_TIMING
  _timers[ DIAGONAL_BASIS_COPY ].tick();
#endif
  copyFinalBases( true /* only one basis */, transpose );
#ifdef DO_TIMING
  _timers[ DIAGONAL_BASIS_COPY ].tock();
#endif

  diagonalBasisProject( A, node, block_idx, false, transpose, useFactor );
}

//////////////////////////////////////////////////////////////////////
// Fixed-rank diagonal decomposition
//////////////////////////////////////////////////////////////////////
void FactorManager::fixedDiagonalDecomposition(
                                const SPARSE_MATRIX::SparseColumnMatrix &A,
                                Supernode &node,
                                int block_idx,
                                int powerIterations,
                                int blockSize,
                                bool transpose,
                                bool useFactor )
{
  // Compute a QR factorization for these columns
#ifdef DO_TIMING
  _timers[ DIAGONAL_QR_FACTORIZATION ].tick();
#endif
  node.lowRankQR( block_idx, _basisWorkspace, blockSize,
                  _qrExtraData, _qrWorkspace, _qrExtraDataSz, transpose );
#ifdef DO_TIMING
  _timers[ DIAGONAL_QR_FACTORIZATION ].tock();
#endif

#ifdef DO_TIMING
  _timers[ DIAGONAL_BASIS_COPY ].tick();
#endif
  appendDiagonalBasis( node, block_idx, blockSize, transpose );
#ifdef DO_TIMING
  _timers[ DIAGONAL_BASIS_COPY ].tock();
#endif
}

//////////////////////////////////////////////////////////////////////
// Adaptive diagonal decomposition
//////////////////////////////////////////////////////////////////////
void FactorManager::adaptiveDiagonalDecomposition(
                                const SPARSE_MATRIX::SparseColumnMatrix &A,
                                Supernode &node,
                                int block_idx,
                                int powerIterations,
                                Real tolerance,
                                int blockSize,
                                bool transpose,
                                bool useFactor )
{
  int                        iteration = 0;
  int                        oldBlockSize;
  int                        blockSizeTmp;
  int                        nonZeroColumns;

  oldBlockSize = blockSize;

  // FIXME: just use the size of the block for now
  //nonZeroColumns = node.diagonalLowRankBlockNonZeroColumns( block_idx );
  nonZeroColumns
    = (int)sqrt(
        3.0 * (Real)node.diagonalLowRankBlocks()[ block_idx ].numColumns() );

  for ( ; ; )
  {
    iteration += 1;

    // Project the current contents of the basis workspace to the
    // complement space of each block's basis
#ifdef DO_TIMING
    _timers[ DIAGONAL_ORTHOGONAL_PROJECTION ].tick();
#endif
    diagonalBasisOrthogonalProjection( node, block_idx, blockSize,
                                       transpose );
#ifdef DO_TIMING
    _timers[ DIAGONAL_ORTHOGONAL_PROJECTION ].tock();
#endif

    // Compute a QR factorization for these columns
#ifdef DO_TIMING
    _timers[ DIAGONAL_QR_FACTORIZATION ].tick();
#endif
    node.lowRankQR( block_idx, _basisWorkspace, blockSize,
                    _qrExtraData, _qrWorkspace, _qrExtraDataSz, transpose );
#ifdef DO_TIMING
    _timers[ DIAGONAL_QR_FACTORIZATION ].tock();
#endif

    // Project this basis in to the complement space of the
    // existing basis
#ifdef DO_TIMING
    _timers[ DIAGONAL_ORTHOGONAL_PROJECTION ].tick();
#endif
    diagonalBasisOrthogonalProjection( node, block_idx, blockSize,
                                       transpose );
#ifdef DO_TIMING
    _timers[ DIAGONAL_ORTHOGONAL_PROJECTION ].tock();
#endif

    // Copy the new basis blocks in to the existing bases
#ifdef DO_TIMING
    _timers[ DIAGONAL_BASIS_COPY ].tick();
#endif
    if ( !appendDiagonalBasis( node, block_idx, blockSize, transpose ) )
    {
      // We have run out of space without converging
      break;
    }
#ifdef DO_TIMING
    _timers[ DIAGONAL_BASIS_COPY ].tock();
#endif

    // Cache the old block size
    blockSizeTmp = oldBlockSize;
    oldBlockSize = blockSize;
    blockSize = blockSizeTmp;

    // Check the current error in this basis and break out of the
    // loop if we are done.  Otherwise, copy to _basisWorkspace to
    // prepare for the next iteration
    if ( !updateActiveDiagonalBlock( A, node, block_idx,
                                     blockSize, tolerance,
                                     transpose, useFactor ) )
    {
      break;
    }

    // Cache the old block size
    blockSizeTmp = oldBlockSize;
    oldBlockSize = blockSize;
    blockSize = blockSizeTmp;

    // Choose the next block size
    blockSize = decompositionBlockSize( nonZeroColumns, iteration,
                                        _currentDescendents.size() > 0,
                                        true /* diagonal block */ );

#if 0
    printf( "Node %d, diagonal block %d, iteration %d, block size = %d\n",
            node.nodeID(), block_idx, iteration, blockSize );
#endif

    // Build a new block column of pre-multiplied random Gaussian data
#ifdef DO_TIMING
    _timers[ DIAGONAL_RANDOM_MULTIPLY ].tick();
#endif
    diagonalBlockRandomMultiply( A, node, block_idx, blockSize,
                                 powerIterations,
                                 transpose, useFactor );
#ifdef DO_TIMING
    _timers[ DIAGONAL_RANDOM_MULTIPLY ].tock();
#endif

    // Cache the old block size
    blockSizeTmp = oldBlockSize;
    oldBlockSize = blockSize;
    blockSize = blockSizeTmp;
  }
}

//////////////////////////////////////////////////////////////////////
// Decomposes off-diagonal block block_idx, *AND* all of its children
// by ecplicitly forming the 2x2 diagonal block for which block[ block_idx ]
// is the off-diagonal part.  Next, this function directly forms low-rank
// approximations for pieces of this block.
//////////////////////////////////////////////////////////////////////
SolverErrorPtr FactorManager::explicitDiagonalDecomposition(
                          const SPARSE_MATRIX::SparseColumnMatrix &A,
                          Supernode &node,
                          int block_idx,
                          int powerIterations,
                          bool transpose,
                          bool useFactor )
{
  // For now, only the version that uses the factor is supported
  TRACE_ASSERT( useFactor, "Not implemented" );

  SolverErrorPtr             error;

  const DenseBlock          &fullBlock = node.nodeDiagonal().fullBlock(
                                                                  block_idx );

#if 0
  printf( "Running explicit decomposition of node %d, block %d, with %d "
          "columns\n", node.nodeID(), block_idx, fullBlock.numRows() );
#endif

  // Build a workspace for explicitly constructing this diagonal block
  MATRIX                     schurComplement( fullBlock.numRows(),
                                              fullBlock.numColumns() );

  // Copy matrix data in to this block
  node.copyDiagonalBlockMatrixData( block_idx, A, schurComplement.data() );

#if 0
  // FIXME: debugging
  schurComplement.write( "super_numeric/step1_sparseMatrix.matrix" );
#endif

  // Subtract contributions from all standard descendents
  for ( int i = 0; i < _currentDescendents.size(); i++ ) {
    int                      descendent_idx = _currentDescendents[ i ];

    if ( isInInteriorBlock(descendent_idx) ) {
      // Will be handled by the interior block loop
      continue;
    }

    constructNodeDiagonalUpdate( node, _factor[ descendent_idx ],
                                 fullBlock, schurComplement.data() );
  }

#if 0
  schurComplement.write( "super_numeric/step2_nodeUpdates.matrix" );
#endif

  // Subtract contributions from all interior block descendents
  for ( int i = 0; i < _currentBlockDescendents.size(); i++ ) {
    int                      block_descendent_idx;
    
    block_descendent_idx = _currentBlockDescendents[ i ];

    constructInteriorBlockDiagonalUpdate( A, node, fullBlock,
                                          block_descendent_idx,
                                          schurComplement.data() );
  }
  
#if 0
  schurComplement.write( "super_numeric/step3_blockUpdates.matrix" );
#endif

  // Pass this explicitly constructed block to the node's diagonal, which
  // will compress off-diagonal blocks and store them in the right places
  error = node.nodeDiagonal().compressExplicitBlock(
                                              block_idx,
                                              schurComplement.data(), transpose,
                                              _diagonalRankEstimator,
                                              _slackData,
                                              _slackDataOffset,
                                              _availableSlackData,
                                              powerIterations );

  TRACE_ASSERT( _availableSlackData >= 0, "Ran out of extra space" );

#if 0
  printf("......done\n");
#endif

  // Pass whatever error we received to the caller
  return error;
}

//////////////////////////////////////////////////////////////////////
// For each low rank block i to be decomposed, using the basis Q
// stored in _basisStorage[ i ], project the contents of _basisWorkspace
// in to the space orthogonal to Q.
//////////////////////////////////////////////////////////////////////
void FactorManager::basisOrthogonalProjection( Supernode &node, int blockSize,
                                               bool transpose )
{
  int                        interaction_idx;
  int                        nRows;
  int                        nCols = node.numColumns();
  int                        currentBlockSize;
  Real                      *baseData = _basisWorkspace;
  Real                      *Qtrans;

  for ( int block_idx = 0; block_idx < _lowRankBlocks.size(); block_idx++ )
  {
    if ( !_lowRankBlockActive[ block_idx ] )
    {
      continue;
    }

    currentBlockSize = ( blockSize > 0 ) ? blockSize : _blockRanks[ block_idx ];

    interaction_idx = _lowRankBlocks[ block_idx ];

    nRows = transpose ? nCols
                      : node.offDiagonal()[ interaction_idx ]._rowList.size();

    BlockStorageData &blockStorage = _basisStorage._dataPtrs[ block_idx ];

    TRACE_ASSERT( nRows == blockStorage._nCols, "Basis size mismatch" );

    Qtrans = blockStorage._data;

    RealWorkspace            workspace( _realWorkspaceManager,
                                        IndexPair( blockStorage._nRows,
                                                   currentBlockSize ) );

    // For this block, multiply Q' by the contents of _basisWorkspace
    MATRIX::gemm( Qtrans, baseData, workspace.workspaceData( 0 ),
                  blockStorage._nRows, blockStorage._nCols,
                  nRows, currentBlockSize,
                  false, false /* do not transpose */ );

    // Multiply Q by the matrix formed above and subtract it from
    // _basisStorage
    MATRIX::gemm( Qtrans, workspace.workspaceData( 0 ), baseData,
                  blockStorage._nRows, blockStorage._nCols,
                  blockStorage._nRows, currentBlockSize,
                  true, /* Transpose Qtrans to get Q */
                  false, /* Do not transpose the workspace */
                  -1.0, 1.0 /* Subtract product */ );

    baseData += nRows * currentBlockSize;
  }
}

//////////////////////////////////////////////////////////////////////
// Similar to the above, but only operates on a basis for a single
// low rank block from the node's diagonal
//////////////////////////////////////////////////////////////////////
void FactorManager::diagonalBasisOrthogonalProjection(
                                        Supernode &node, int block_idx,
                                        int blockSize,
                                        bool transpose )
{
  const DenseBlock &block = node.diagonalLowRankBlocks()[ block_idx ];

  int                        nRows = transpose ? block.numColumns()
                                               : block.numRows();
  Real                      *Qtrans;

  BlockStorageData          &blockStorage = _basisStorage._dataPtrs[ 0 ];

  TRACE_ASSERT( nRows == blockStorage._nCols, "Basis size mismatch" );

  Qtrans = blockStorage._data;

  RealWorkspace            workspace( _realWorkspaceManager,
                                      IndexPair( blockStorage._nRows,
                                                 blockSize ) );

  // For this block, multiply Q' by the contents of _basisWorkspace
  MATRIX::gemm( Qtrans, _basisWorkspace, workspace.workspaceData( 0 ),
                blockStorage._nRows, blockStorage._nCols,
                nRows, blockSize,
                false, false /* do not transpose */ );

  // Multiply Q by the matrix formed above and subtract it from
  // _basisStorage
  MATRIX::gemm( Qtrans, workspace.workspaceData( 0 ), _basisWorkspace,
                blockStorage._nRows, blockStorage._nCols,
                blockStorage._nRows, blockSize,
                true, /* Transpose Qtrans to get Q */
                false, /* Do not transpose the workspace */
                -1.0, 1.0 /* Subtract product */ );
}

//////////////////////////////////////////////////////////////////////
// Append bases stored in _basisStorage based on the contents
// of _basisWorkspace
//////////////////////////////////////////////////////////////////////
void FactorManager::appendOffDiagonalBases(
                              Supernode &node, int blockSize,
                              int iteration,
                              const SPARSE_MATRIX::SparseColumnMatrix *A,
                              bool transpose )
{
  int                        interaction_idx;
  int                        nRows;
  int                        nCols = node.numColumns();
  const Real                *baseInputData;
  Real                      *baseOutputData;
  int                        sz;
  int                        finalBlockSize;
  int                        remainingSize;
  int                        currentBlockSize;

  baseInputData = _basisWorkspace;

  for ( int block_idx = 0; block_idx < _lowRankBlocks.size(); block_idx++ )
  {
    if ( !_lowRankBlockActive[ block_idx ] )
    {
      continue;
    }

    bool                     copyBasis = true;

    currentBlockSize = ( blockSize > 0 ) ? blockSize : _blockRanks[ block_idx ];
    finalBlockSize = currentBlockSize;

    interaction_idx = _lowRankBlocks[ block_idx ];

    nRows = transpose ? nCols
                      : node.offDiagonal()[ interaction_idx ]._rowList.size();

    BlockStorageData &blockStorage = _basisStorage._dataPtrs[ block_idx ];

    // Note that we store the transpose of the basis, due
    // to row major ordering
    TRACE_ASSERT( nRows == blockStorage._nCols,
                  "Basis size mismatch" );

    sz = blockStorage._nRows * blockStorage._nCols;

    baseOutputData = blockStorage._data + sz;

    blockStorage._nRows += currentBlockSize;

    sz = blockStorage._nRows * blockStorage._nCols;

    // Error checking
    if ( sz > blockStorage._maxSize || blockStorage._nRows >= nRows )
    {
      printf( "Warning: Basis for block %d of node %d has exceeded "
              "maximum size.\n", interaction_idx, node.nodeID() );
      printf( "nRows = %d, nCols = %d, basis size = %d\n",
              nRows, node.numColumns(), blockStorage._nRows );
      printf( "Block size: %d\n", currentBlockSize );
      printf( "Iteration number: %d\n", iteration );
      printf( "Maximum basis size: %d x %d\n",
              blockStorage._maxSize / nRows, nRows );
      printf( "Current error: %f\n",
              _blockErrors[ block_idx ][ 1 ] / _blockErrors[ block_idx ][ 0 ] );

      // Flag the block as full
      blockStorage._full = true;

      // Roll back the changes and just take as much as we can.
      // If the basis size for some reason exceeds the number of rows in
      // the matrix, then replace the basis with the identity matrix.
      if ( blockStorage._nRows >= nRows )
      {
        blockStorage._nRows = nRows;

        // Replace the entire basis with the identity
        MATRIX::ident( blockStorage._data, nRows );

        printf( "Setting basis to %d x %d identity\n", nRows, nRows );

        copyBasis = false;
      }
      else
      {
        blockStorage._nRows -= currentBlockSize;

        // Take as many columns as we can
        remainingSize = blockStorage._maxSize
          - blockStorage._nRows * blockStorage._nCols;

        TRACE_ASSERT( remainingSize % blockStorage._nCols == 0 );

        // This is how many more basis columns we can accomodate
        remainingSize /= blockStorage._nCols;

        finalBlockSize = remainingSize;

        TRACE_ASSERT( finalBlockSize < currentBlockSize );

        printf( "Adding %d columns instead\n", finalBlockSize );

        blockStorage._nRows += finalBlockSize;

        printf( "Final basis size = %d\n", blockStorage._nRows );
      }

#if 0
      if ( A )
      {
        cout << "Writing" << endl;
        writeLowRankBlocks( *A, node );
      }

      abort();
#endif
    }

#if 0
    {
      char buf[ 1024 ];
      if ( _useInteriorBlocks ) {
        sprintf( buf, "node_%d_interaction_%d_basiscopy_int.matrix",
                 node.nodeID(), _lowRankBlocks[ block_idx ] );
      } else {
        sprintf( buf, "node_%d_interaction_%d_basiscopy.matrix",
                 node.nodeID(), _lowRankBlocks[ block_idx ] );
      }

      MATRIX::write( baseInputData, nRows, finalBlockSize, buf );
    }
#endif

    // Copy the transpose of the new block in to the basis
    if ( copyBasis )
    {
      MATRIX::transposeBLAS( baseOutputData, baseInputData, nRows,
                             finalBlockSize,
                             blockStorage._nCols, /* leading dim. for output */
                             currentBlockSize /* leading dim. for input */ );
    }

#if 0
    MATRIX tmp( blockStorage._nRows, blockStorage._nCols, blockStorage._data );
    tmp.write( "super_numeric/testBasis.matrix" );
#endif

    baseInputData += nRows * currentBlockSize;
  }
}

//////////////////////////////////////////////////////////////////////
// Appends basis data when building a single off-diagonal decomposition
//////////////////////////////////////////////////////////////////////
void FactorManager::appendOffDiagonalBasis(
                              Supernode &node, int blockSize,
                              int iteration,
                              const SPARSE_MATRIX::SparseColumnMatrix *A,
                              bool transpose )
{
  int                        nRows;
  int                        nCols = node.numColumns();
  const Real                *baseInputData;
  Real                      *baseOutputData;
  int                        sz;
  int                        finalBlockSize;
  int                        remainingSize;
  int                        currentBlockSize;

  bool                       copyBasis = true;

  baseInputData = _basisWorkspace;

  currentBlockSize = ( blockSize > 0 ) ? blockSize : _blockRanks[ 0 ];
  finalBlockSize = currentBlockSize;

  nRows = transpose ? nCols : node.countLowRankRows();

  BlockStorageData          &blockStorage = _basisStorage._dataPtrs[ 0 ];

  // Note that we store the transpose of the basis, due to row major
  // ordering
  TRACE_ASSERT( nRows == blockStorage._nCols );

  sz = blockStorage._nRows * blockStorage._nCols;

  baseOutputData = blockStorage._data + sz;

  blockStorage._nRows += currentBlockSize;

  sz = blockStorage._nRows * blockStorage._nCols;

  // Error checking
  if ( sz > blockStorage._maxSize || blockStorage._nRows >= nRows )
  {
    printf( "Warning: Off-diagonal basis for node %d has exceeded "
            "maximum size.\n", node.nodeID() );

    printf( "nRows = %d, nCols = %d, basis size = %d\n",
            nRows, node.numColumns(), blockStorage._nRows );
    printf( "Block size: %d\n", currentBlockSize );
    printf( "Iteration number: %d\n", iteration );
    printf( "Maximum basis size: %d x %d\n",
            blockStorage._maxSize / nRows, nRows );
    printf( "Current error: %f\n",
            _blockErrors[ 0 ][ 1 ] / _blockErrors[ 0 ][ 0 ] );

    // Flag the block as full
    blockStorage._full = true;

    // Roll back the changes and just take as much as we can.
    // If the basis size for some reason exceeds the number of rows in
    // the matrix, then replace the basis with the identity matrix.
    if ( blockStorage._nRows >= nRows ) {
      blockStorage._nRows = nRows;

      // Replace the entire basis with the identity
      MATRIX::ident( blockStorage._data, nRows );

      printf( "Setting basis to %d x %d identity\n", nRows, nRows );

      copyBasis = false;
    } else {
      blockStorage._nRows -= currentBlockSize;

      // Take as many columns as we can
      remainingSize = blockStorage._maxSize
        - blockStorage._nRows * blockStorage._nCols;

      TRACE_ASSERT( remainingSize % blockStorage._nCols == 0 );

      // This is how many more basis columns we can accomodate
      remainingSize /= blockStorage._nCols;

      finalBlockSize = remainingSize;

      TRACE_ASSERT( finalBlockSize < currentBlockSize );

      printf( "Adding %d columns instead\n", finalBlockSize );

      blockStorage._nRows += finalBlockSize;

      printf( "Final basis size = %d\n", blockStorage._nRows );
    }
  }

  // Copy the transpose of the new block in to the basis
  if ( copyBasis ) {
    MATRIX::transposeBLAS( baseOutputData, baseInputData, nRows,
                           finalBlockSize,
                           blockStorage._nCols, /* leading dim. for output */
                           currentBlockSize /* leading dim. for input */ );
  }
}

//////////////////////////////////////////////////////////////////////
// Appends to a single basis for a block on the given node's diagonal
//////////////////////////////////////////////////////////////////////
bool FactorManager::appendDiagonalBasis( Supernode &node, int block_idx,
                                         int blockSize,
                                         bool transpose )
{
  const DenseBlock &block = node.diagonalLowRankBlocks()[ block_idx ];

  int                        nRows = transpose ? block.numColumns()
                                               : block.numRows();
  int                        sz;
  int                        finalBlockSize = blockSize;
  int                        remainingSize;

  bool                       spaceRemaining = true;

  Real                      *baseOutputData;

  BlockStorageData          &blockStorage = _basisStorage._dataPtrs[ 0 ];

  TRACE_ASSERT( nRows == blockStorage._nCols, "Basis size mismatch" );

  sz = blockStorage._nRows * blockStorage._nCols;

  baseOutputData = blockStorage._data + sz;

  blockStorage._nRows += blockSize;

  sz = blockStorage._nRows * blockStorage._nCols;

  // Error checking
  if ( sz > blockStorage._maxSize )
  {
#if 0
    printf( "Basis for diagonal block %d of node %d has exceeded max size.\n",
            block_idx, node.nodeID() );
    printf( "nRows = %d, nCols = %d, basis size = %d\n",
            nRows, block.numColumns(), blockStorage._nRows );
    printf( "Maximum basis size: %d x %d\n",
            blockStorage._maxSize / nRows, nRows );
    printf( "Current error: %f\n",
            _blockErrors[ 0 ][ 1 ] / _blockErrors[ 0 ][ 0 ] );
    printf( "Factorization failed\n" );

#if 0
    // Dump out the bad basis
    blockStorage._nRows -= blockSize;
    MATRIX tmp( blockStorage._nRows, blockStorage._nCols, blockStorage._data );
    tmp.write( "super_numeric/badBasis.matrix" );
#endif

    abort();
#endif
    printf( "Warning: Basis for diagonal block %d of node %d "
            "has exceeded max size.\n",
            block_idx, node.nodeID() );
    printf( "nRows = %d, nCols = %d, basis size = %d\n",
            nRows, block.numColumns(), blockStorage._nRows );
    printf( "Previous basis size = %d\n",
            blockStorage._nRows - blockSize );
    printf( "Maximum basis size: %d x %d\n",
            blockStorage._maxSize / nRows, nRows );
    printf( "Current error: %f\n\n",
            _blockErrors[ 0 ][ 1 ] / _blockErrors[ 0 ][ 0 ] );

    // Roll back the changes and just take as much as we can
    blockStorage._nRows -= blockSize;

    remainingSize = blockStorage._maxSize
      - blockStorage._nRows * blockStorage._nCols;

    TRACE_ASSERT( remainingSize % blockStorage._nCols == 0 );

    // This is how many more basis columns we can accomodate
    remainingSize /= blockStorage._nCols;

    finalBlockSize = remainingSize;

    TRACE_ASSERT( finalBlockSize < blockSize );

    printf( "Adding %d columns instead\n", finalBlockSize );

    blockStorage._nRows += finalBlockSize;

    printf( "Final basis size = %d\n", blockStorage._nRows );

    spaceRemaining = false;

#if 0
    writeLowRankDiagonalBlock( _A, node, block_idx );
#endif
  }

  // Copy the transpose of the new block in to the basis
  MATRIX::transposeBLAS( baseOutputData, _basisWorkspace, nRows,
                         finalBlockSize,
                         blockStorage._nCols, /* leading dimension for output */
                         blockSize /* leading dimension for input */ );

  return spaceRemaining;
}

//////////////////////////////////////////////////////////////////////
// Copies finalized adaptive bases to _decompWorkspace in preparation
// for the second pass of low-rank decomposition
//////////////////////////////////////////////////////////////////////
void FactorManager::copyFinalBases( bool singleBasis, bool transpose )
{
  int                       remainingSize;
  Real                     *baseData;

  remainingSize = transpose ? _decompMultTransWorkspaceSz : _decompWorkspaceSz;
  baseData = transpose ? _decompMultTransWorkspace : _decompWorkspace;

  int                       numBases = singleBasis ? 1 : _lowRankBlocks.size();

  long int                  totalSize = 0;

  for ( int block_idx = 0; block_idx < numBases; block_idx++ )
  {
    const BlockStorageData &blockStorage = _basisStorage._dataPtrs[ block_idx ];

    remainingSize -= blockStorage._nRows * blockStorage._nCols;

    TRACE_ASSERT( remainingSize >= 0, "Out of room in decomp workspace" );

    MATRIX::transposeBLAS( baseData, blockStorage._data, 
                           blockStorage._nRows, blockStorage._nCols );

    baseData += blockStorage._nRows * blockStorage._nCols;

    totalSize += blockStorage._nRows * blockStorage._nCols;

    // Also, set the block rank to the rank of the adaptive basis
    _blockRanks[ block_idx ] = blockStorage._nRows;
  }

  TRACE_ASSERT( totalSize <= _decompWorkspaceSz,
                "Decomposition workspace capacity exceeded" );
}

//////////////////////////////////////////////////////////////////////
// For each descendent of the node currently being handled identify
// row/column ranges of the descendent that must be considered to
// form each low-rank block in the current node's main diagonal.
//////////////////////////////////////////////////////////////////////
void FactorManager::findDiagonalBlockRanges( const Supernode &currentNode )
{
  int                        descendent_idx;
  int                        start_idx;
  int                        row_idx, col_idx;
  int                        full_row_idx;

  int                        startRow, endRow;
  int                        startCol, endCol;

  const vector<DenseBlock> &denseBlocks = currentNode.diagonalLowRankBlocks();

  _inverseRowMap.resize( currentNode.numColumns() );

  for ( int descendent_ptr = 0; descendent_ptr < _currentDescendents.size();
        descendent_ptr++ )
  {
    descendent_idx = _currentDescendents[ descendent_ptr ];

    const Supernode         &descendent = _factor[ descendent_idx ];

    start_idx = _nextInteractionCache[ descendent_idx ];

    const SupernodeInteraction &interaction
                            = descendent.offDiagonal()[ start_idx ];

    vector<DenseBlock> &blockRanges = _diagonalBlockRanges[ descendent_ptr ];

    blockRanges.resize( denseBlocks.size() );

    for ( int block_idx = 0; block_idx < denseBlocks.size(); block_idx++ )
    {
      const DenseBlock      &block = denseBlocks[ block_idx ];

      DenseBlock            &descendentRanges = blockRanges[ block_idx ];

      findEntryRangeIntersection( interaction._rowList, block._rowRange,
                                  descendentRanges._rowRange );

      findEntryRangeIntersection( interaction._rowList, block._columnRange,
                                  descendentRanges._columnRange );
    }
  }
}

//////////////////////////////////////////////////////////////////////
// Gets maximum row/column counts, ranks, block sizes, etc.
// for all low rank diagonal blocks in the given node
//////////////////////////////////////////////////////////////////////
void FactorManager::diagonalBlockSizes( const Supernode &node,
                                        int &maxRows, int &maxColumns,
                                        int &maxRank )
{
  const vector<DenseBlock>  &blocks = node.diagonalLowRankBlocks();

  for ( int block_idx = 0; block_idx < blocks.size(); block_idx++ )
  {
    const DenseBlock        &block = blocks[ block_idx ];

    maxRows = max( maxRows, block.numRows() );
    maxColumns = max( maxColumns, block.numColumns() );
    maxRank = max( maxRank, diagonalMaxBlockRank( block.numRows(),
                                                  block.numColumns() ) );
  }
}

//////////////////////////////////////////////////////////////////////
// Yes, this is annoying
//////////////////////////////////////////////////////////////////////
void FactorManager::cacheBlockRanks()
{
  int                        rankTemp;

  if ( _blockRanksCache.size() < _blockRanks.size() )
  {
    _blockRanksCache.resize( _blockRanks.size(), 0 );
  }

  for ( int block_idx = 0; block_idx < _blockRanks.size(); block_idx++ )
  {
    rankTemp = _blockRanks[ block_idx ];
    _blockRanks[ block_idx ] = _blockRanksCache[ block_idx ];
    _blockRanksCache[ block_idx ] = rankTemp;
  }
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void FactorManager::initTimers()
{
  _timers.clear();

  for ( int timer_idx = 0; timer_idx < NUM_TIMERS; timer_idx++ )
  {
    _timers.push_back( Timer( TIMER_NAMES[ timer_idx ] ) );
  }
}

//////////////////////////////////////////////////////////////////////
// Assuming that sizes have been set in the extended part of the
// factorization, compute a new ordering for this part of the
// factorization using a weighted minimum degree ordering
//////////////////////////////////////////////////////////////////////
void FactorManager::buildExtendedPermutation( IntArray &permutation,
                                              IntArray &inversePermuation )
{
  // Track the fill in resulting from elimination of extended
  // nodes in the permuted order
  vector<set<int> >          extendedInteractions;
  IntArray                   extendedSizes;
  int                        numStandardNodes;
  
  numStandardNodes = _factor.size() - _numExtendedNodes;

  // Start by getting the basic structure of the Schur complement
  buildExtendedSchurComplement( extendedInteractions );

  for ( int node_idx = numStandardNodes; node_idx < _factor.size(); node_idx++ )
  {
    extendedSizes.push_back( _factor[ node_idx ].numColumns() );

    TRACE_ASSERT( extendedSizes.back() > 0, "Zero node size" );
  }

#if 0
  // Try reordering the Schur complement according to a minimum degree
  // reordering
  MinimumDegree::BlockMinimumDegree( extendedInteractions, extendedSizes,
                                     permutation, inversePermuation );

  MinimumDegree::PrintBlockOrderingStats( extendedInteractions,
                                          extendedSizes,
                                          permutation, inversePermuation );
#endif
}

//////////////////////////////////////////////////////////////////////
// Once all standard nodes have been handled, call this to
// build a symbolic representation of the Schur complement formed
// by eliminating just the standard nodes in the system.
//////////////////////////////////////////////////////////////////////
void FactorManager::buildExtendedSchurComplement(
                            vector<set<int> > &extendedInteractions )
{
  // Build an interaction list for each extended node
  int                        numStandardNodes;
  int                        startColumn, endColumn;

  extendedInteractions.clear();
  extendedInteractions.resize( _numExtendedNodes );
  
  numStandardNodes = _factor.size() - _numExtendedNodes;
  startColumn = _factor[ numStandardNodes ].startColumn();
  endColumn = _factor.back().endColumn();

  TRACE_ASSERT( startColumn >= 0, "Not all standard nodes have been factored" );
  TRACE_ASSERT( endColumn >= 0, "Not all standard nodes have been factored" );

  for ( int node_idx = 0; node_idx < numStandardNodes; node_idx++ )
  {
    const Supernode         &node = _factor[ node_idx ];

    for ( int interaction_idx = node.firstExtendedInteraction();
          interaction_idx < node.offDiagonal().size(); interaction_idx++ )
    {
      const SupernodeInteraction &interaction
        = node.offDiagonal()[ interaction_idx ];

      // Figure out what fill in this induces in the Schur complement
      int      nodeID = interaction._nodeID - numStandardNodes;

      for ( int next_interaction_idx = interaction_idx + 1;
            next_interaction_idx < node.offDiagonal().size();
            next_interaction_idx++ )
      {
        const SupernodeInteraction &nextInteraction
          = node.offDiagonal()[ next_interaction_idx ];

        int    nextNodeID = nextInteraction._nodeID - numStandardNodes;

        extendedInteractions[ nodeID ].insert( nextNodeID );
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////
// Identify "interior" block ranges; that is, ranges of supernodes
// whose off-diagonal will not be explicitly stored
//////////////////////////////////////////////////////////////////////
void FactorManager::identifyInteriorBlocks()
{
  int                        start_idx = 0;
  int                        end_idx = 0;

  _useInteriorBlocks = true;

  _interiorBlocks.clear();
  _interiorBlockMap.resize( _factor.size(), EMPTY );

  while ( start_idx < _factor.size() )
  {
    end_idx = start_idx;

    _interiorBlockMap[ end_idx ] = EMPTY;

    while ( end_idx < _factor.size()
         && _factor[ end_idx ].type() == STANDARD_NODE
         && !_factor[ end_idx ].compressOffDiagonal() )
    {
      _interiorBlockMap[ end_idx ] = (int)_interiorBlocks.size();
      end_idx++;
    }

    if ( end_idx > start_idx )
    {
      // Add a new node range
      _interiorBlocks.push_back( InteriorBlock( start_idx, end_idx - 1 ) );

      start_idx = end_idx;
    }
    else
    {
      start_idx += 1;
    }
  }

  _currentBlockDescendents.reserve( _interiorBlocks.size() );
  _blockDescendentFound.resize( _interiorBlocks.size(), false );

  buildInteriorBlockInteractions();

#ifdef FINALIZE_INTERIOR_BLOCKS
  // Go through each node in the interior block set and flag
  // interactions with nodes outside of that block as implicit (ie.
  // we will no longer form them explicitly
  for ( int block_idx = 0; block_idx < _interiorBlocks.size(); block_idx++ ) {
    const InteriorBlock     &block = _interiorBlocks[ block_idx ];

    for ( int node_idx = block._nodeRange.first;
          node_idx <= block._nodeRange.second; node_idx++ )
    {
      Supernode             &node = _factor[ node_idx ];

      for ( int interaction_idx = 0;
            interaction_idx < node.offDiagonal().size();
            interaction_idx++ )
      {
        if ( node.offDiagonal()[ interaction_idx ]._nodeID
                  > block._nodeRange.second )
        {
          node.flagImplicitInteraction( interaction_idx );
        }
      }
    }
  }
#endif

  // FIXME: debugging
  //
  verifyInteriorBlocks();
}

//////////////////////////////////////////////////////////////////////
// Builds interaction sets for interior blocks
//////////////////////////////////////////////////////////////////////
void FactorManager::buildInteriorBlockInteractions()
{
  for ( int interior_block_idx = 0;
        interior_block_idx < _interiorBlocks.size();
        interior_block_idx++ )
  {
    InteriorBlock           &block = _interiorBlocks[ interior_block_idx ];

    buildSingleInteriorBlockInteractions( block );
  }
}

//////////////////////////////////////////////////////////////////////
// Helper function for the above
//////////////////////////////////////////////////////////////////////
void FactorManager::buildSingleInteriorBlockInteractions( InteriorBlock &block )
{
  set<int>                   allInteractionBlocks;
  IntArray                   interactionBlockMap;
  int                        idx = 0;
  int                        lastColumn;

  // Figure other supernodes this block interacts with
  for ( int node_idx = block._nodeRange.first;
        node_idx <= block._nodeRange.second; node_idx++ )
  {
    const Supernode         &node = _factor[ node_idx ];

    for ( int interaction_idx = 0; interaction_idx < node.offDiagonal().size();
          interaction_idx++ )
    {
      const SupernodeInteraction &interaction
                                    = node.offDiagonal()[ interaction_idx ];

      // Ignore nodes in this block's range
      if ( interaction._nodeID > block._nodeRange.second )
      {
        allInteractionBlocks.insert( interaction._nodeID );
      }
    }
  }

  // Build an inverse map from the indices of nodes this block interacts
  // with to indices in this block's interaction set
  block._interactions.resize( allInteractionBlocks.size() );
  block._interactionMap.resize( _factor.size(), -1 );

  idx = 0;
  for ( set<int>::const_iterator iter = allInteractionBlocks.begin();
        iter != allInteractionBlocks.end(); iter++ )
  {
    block._interactionMap[ *iter ] = idx;

    idx++;
  }

  // Fill in all of the interaction blocks
  for ( int node_idx = block._nodeRange.first;
        node_idx <= block._nodeRange.second; node_idx++ )
  {
    const Supernode         &node = _factor[ node_idx ];

    for ( int interaction_idx = 0; interaction_idx < node.offDiagonal().size();
          interaction_idx++ )
    {
      const SupernodeInteraction &interaction
                                    = node.offDiagonal()[ interaction_idx ];

      // Ignore nodes in this block's range
      if ( interaction._nodeID > block._nodeRange.second )
      {
        int                  block_interaction_idx;

        block_interaction_idx = block._interactionMap[ interaction._nodeID ];

        TRACE_ASSERT( block_interaction_idx >= 0 );

        InteriorBlock::Interaction &interaction = 
                                block._interactions[ block_interaction_idx ];

        // Add the appropriate [node index, interaction_index] pair
        interaction.push_back( IndexPair( node_idx, interaction_idx ) );
      }
    }
  }

  // Identify off-diagonal row indices for every sparse matrix column in
  // this interior block
  block._firstOffDiagonalRows.clear();
  lastColumn = _factor[ block._nodeRange.second ].endColumn();

  for ( int node_idx = block._nodeRange.first;
        node_idx <= block._nodeRange.second; node_idx++ )
  {
    const Supernode         &node = _factor[ node_idx ];

    for ( int col_idx = node.startColumn(); col_idx <= node.endColumn();
          col_idx++ )
    {
      block._firstOffDiagonalRows.push_back( 0 );

      int                   &firstRow = block._firstOffDiagonalRows.back();

      for ( int row_ptr = _A._p[ col_idx ]; row_ptr < _A._p[ col_idx + 1 ];
            row_ptr++ )
      {
        if ( _A._i[ row_ptr ] <= lastColumn ) {
          firstRow += 1;
        } else {
          break;
        }
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////
// Debugging routine to make sure that the interior blocks are
// working the way they should
//////////////////////////////////////////////////////////////////////
void FactorManager::verifyInteriorBlocks() const
{
  printf( "System has %d interior blocks\n", (int)_interiorBlocks.size() );

  for ( int interior_block_idx = 0; interior_block_idx < _interiorBlocks.size();
        interior_block_idx++ )
  {
    const InteriorBlock      &block = _interiorBlocks[ interior_block_idx ];

    for ( int node_idx = block._nodeRange.first;
          node_idx <= block._nodeRange.second;
          node_idx++ )
    {
      const Supernode       &node = _factor[ node_idx ];

      for ( int interaction_idx = 0;
            interaction_idx < node.offDiagonal().size();
            interaction_idx++ )
      {
        const SupernodeInteraction &interaction
                                      = node.offDiagonal()[ interaction_idx ];

        // This interaction should either be with another node in the
        // interior block range, or it should be with a non-interior
        // node (eg. a compressed node, or an extended node)
        if ( interaction._nodeID > block._nodeRange.second
          && _interiorBlockMap[ interaction._nodeID ] != EMPTY )
        {
          printf( "Node %d from interior block range [%d, %d] "
                  "interacts with node %d from interior block %d\n",
                  node.nodeID(),
                  block._nodeRange.first, block._nodeRange.second,
                  interaction._nodeID,
                  _interiorBlockMap[ interaction._nodeID ] );

          TRACE_ASSERT( interaction._nodeID <= block._nodeRange.second
                     || _interiorBlockMap[ interaction._nodeID ] == EMPTY );
        }
#if 0
        else
        {
          printf( "Node %d from interior block range [%d, %d] "
                  "interacts with node %d from interior block %d\n",
                  node.nodeID(),
                  block._nodeRange.first, block._nodeRange.second,
                  interaction._nodeID,
                  _interiorBlockMap[ interaction._nodeID ] );


        }
#endif
      }
    }
  }
}

#if 0
//////////////////////////////////////////////////////////////////////
// Node list is the list of supernodes to use in this multiplication,
// relative to the entire matrix.
//
// subNodeList is the list of nodes (of the same size) relative to
// the particular interior block
//
// If A1 and A2 are the sparse matrix blocks associated with interaction
// index 1 and 2, respectively, then we form
//      A2 * ( L^{-T} * ( L^{-1} * ( A1^{T} * G ) ) )
//
// Performs the multiplication necessary to form a low rank block
// depending on the given interior block, without explicitly forming
// the interior block.
//////////////////////////////////////////////////////////////////////
void FactorManager::interiorBlockSchurMultiply(
                            const SPARSE_MATRIX::SparseColumnMatrix &A,
                            const IntArray &nodeList,
                            int interior_block_idx,
                            int interaction_idx1, int interaction_idx2,
                            const Real *G,
                            int nRowsG, int nColsG,
                            Real *multWorkspaceInitial,
                            Real *multWorkspaceFinal,
                            IndexRange rowRange,
                            IndexRange columnRange,
                            bool transpose, bool left )
{
  int                        fullColumnCount = 0;
  Real                      *baseInputData;
  Real                      *baseOutputData;

  int                        sparseRowStart, sparseColStart;
  int                        numSparseRows, numSparseCols;

  int                        base_node_idx;
  int                        interaction_node_idx;

#ifdef DO_TIMING
  _timers[ INTERIOR_BLOCK_MULT_SETUP ].tick();
#endif

  const InteriorBlock       &block = _interiorBlocks[ interior_block_idx ];

  // Choose which interaction is which based on whether or not
  // we want to compute the transpose
  const InteriorBlock::Interaction  &interaction1
      = transpose ? block._interactions[ interaction_idx2 ] :
                    block._interactions[ interaction_idx1 ];
  const InteriorBlock::Interaction  &interaction2
      = transpose ? block._interactions[ interaction_idx1 ] :
                    block._interactions[ interaction_idx2 ];

  // Swap the row and column sub-ranges if we are transposing
  if ( transpose )
  {
    IndexRange               tmp = rowRange;

    rowRange = columnRange;
    columnRange = tmp;
  }

#if 0
  printf( "Interior multiply over column range [ %d, %d ]\n",
          _factor[ block._nodeRange.first ].startColumn(),
          _factor[ block._nodeRange.second ].endColumn() );
#endif

  TRACE_ASSERT( interaction1.size() > 0 && interaction2.size() > 0 );

  // Index of the node in which we are building a low rank block, and
  // index of the node with which that low rank interaction is interacting.
  {
    const Supernode &firstNode = _factor[ interaction1[ 0 ].first ];
    const SupernodeInteraction &firstInteraction
                     = firstNode.offDiagonal().at( interaction1[ 0 ].second );

    base_node_idx = firstInteraction._nodeID;
  }

  {
    const Supernode &secondNode = _factor[ interaction2[ 0 ].first ];
    const SupernodeInteraction &secondInteraction
                     = secondNode.offDiagonal().at( interaction2[ 0 ].second );

    interaction_node_idx = secondInteraction._nodeID;
  }

  // Get the full number of columns necessary to store in workspaces
  for ( int node_idx = 0; node_idx < nodeList.size(); node_idx++ )
  {
    fullColumnCount += _factor.at( nodeList[ node_idx ] ).numColumns();
  }

  // Clear the workspace
  MATRIX::clear( multWorkspaceInitial, fullColumnCount, nColsG );

  // Transpose multiply the matrix A1 with G
  sparseRowStart = _factor[ base_node_idx ].startColumn();
  numSparseRows = _factor[ base_node_idx ].numColumns();

#ifdef DO_TIMING
  _timers[ INTERIOR_BLOCK_MULT_SETUP ].tock();
#endif

  // Restrict to the given column range, if it is valid
  if ( columnRange.first >= 0 && columnRange.second >= 0 )
  {
    sparseRowStart += columnRange.first;

    TRACE_ASSERT( columnRange.first + range_size( columnRange ) <= numSparseRows,
                  "Invalid row range" );

    numSparseRows = range_size( columnRange );
  }

  TRACE_ASSERT( nRowsG == numSparseRows );

#ifdef DO_TIMING
  _timers[ INTERIOR_BLOCK_MULT_SPARSE ].tick();
#endif
  sparseSubTransposeMultiply( A, nodeList, sparseRowStart, numSparseRows,
                              G, nColsG, multWorkspaceInitial, 1.0 );
#ifdef DO_TIMING
  _timers[ INTERIOR_BLOCK_MULT_SPARSE ].tock();
#endif

  // FIXME: debugging
  //
  // build our own workspace for now
#ifdef DO_TIMING
  _timers[ INTERIOR_BLOCK_MULT_WORKSPACE ].tick();
#endif
  MATRIX workspace( fullColumnCount, nColsG );
#ifdef DO_TIMING
  _timers[ INTERIOR_BLOCK_MULT_WORKSPACE ].tock();
#endif

  // Triangular solves
#ifdef DO_TIMING
  _timers[ INTERIOR_BLOCK_MULT_FORWARD_SOLVE ].tick();
#endif
  forwardSubSolve( nodeList, nColsG, multWorkspaceInitial, workspace.data() );
#ifdef DO_TIMING
  _timers[ INTERIOR_BLOCK_MULT_FORWARD_SOLVE ].tock();
#endif

#ifdef DO_TIMING
  _timers[ INTERIOR_BLOCK_MULT_BACKWARD_SOLVE ].tick();
#endif
  backwardSubSolve( nodeList, nColsG, multWorkspaceInitial, workspace.data() );
#ifdef DO_TIMING
  _timers[ INTERIOR_BLOCK_MULT_BACKWARD_SOLVE ].tock();
#endif

  // Standard multiply of matrix A2 with the result so far
  sparseRowStart = _factor[ interaction_node_idx ].startColumn();
  numSparseRows = _factor[ interaction_node_idx ].numColumns();

  // Restrict to the given row range, if it is valid
  if ( rowRange.first >= 0 && rowRange.second >= 0 )
  {
    sparseRowStart += rowRange.first;

    TRACE_ASSERT( rowRange.first + range_size( rowRange ) <= numSparseRows,
                  "Invalid row range" );

    numSparseRows = range_size( rowRange );
  }

  MATRIX::clear( multWorkspaceFinal, numSparseRows, nColsG );

#ifdef DO_TIMING
  _timers[ INTERIOR_BLOCK_MULT_SPARSE ].tick();
#endif
  sparseSubMultiply( A, nodeList, sparseRowStart, numSparseRows,
                     multWorkspaceInitial, nColsG, multWorkspaceFinal, 1.0 );
#ifdef DO_TIMING
  _timers[ INTERIOR_BLOCK_MULT_SPARSE ].tock();
#endif
}
#endif

//////////////////////////////////////////////////////////////////////
// Node list is the list of supernodes to use in this multiplication,
// relative to the entire matrix.
//
// subNodeList is the list of nodes (of the same size) relative to
// the particular interior block
//
// Performs the multiplication necessary to form a low rank block
// depending on the given interior block, without explicitly forming
// the interior block.
//
// inverseRowList maps rows in the "full" matrix that we are forming to
// those in a matrix in which some of these rows may be excluded
//
// Optionally restrict to a row and column range (relative to the
// dimensions of the product matrix).  Note that if meaningful row/column
// ranges are provided, this will ignore the provided inverse row list and
// assume that all rows in the desired range are required.
//////////////////////////////////////////////////////////////////////
#if 0
void FactorManager::interiorBlockSchurMultiply(
                            const SPARSE_MATRIX::SparseColumnMatrix &A,
                            const IntArray &nodeList,
                            int interior_block_idx,
                            int interaction_idx1, int interaction_idx2,
                            const Real *G,
                            int nRowsG, int nColsG,
                            const IntArray &inverseRowList, int nRows,
                            Real *multWorkspaceInitial,
                            Real *outputData,
                            IndexRange rowRange,
                            IndexRange columnRange,
                            bool transpose, bool left )
{
  int                        fullColumnCount = 0;

  int                        sparseRowStart, sparseColStart;
  int                        numSparseRows, numSparseCols;

  int                        base_node_idx;
  int                        interaction_node_idx;

  bool                       useRowList = true;

#ifdef DO_TIMING
  _timers[ INTERIOR_BLOCK_MULT_SETUP ].tick();
#endif
  TIMING_START( "Interior block mult setup" );

  const InteriorBlock       &block = _interiorBlocks[ interior_block_idx ];

  // Choose which interaction is which based on whether or not
  // we want to compute the transpose
  const InteriorBlock::Interaction  &interaction1
      = transpose ? block._interactions[ interaction_idx2 ] :
                    block._interactions[ interaction_idx1 ];
  const InteriorBlock::Interaction  &interaction2
      = transpose ? block._interactions[ interaction_idx1 ] :
                    block._interactions[ interaction_idx2 ];

  // Swap the row and column sub-ranges if we are transposing
  if ( transpose )
  {
    IndexRange               tmp = rowRange;

    rowRange = columnRange;
    columnRange = tmp;
  }

  TRACE_ASSERT( interaction1.size() > 0 && interaction2.size() > 0 );

  // Index of the node in which we are building a low rank block, and
  // index of the node with which that low rank interaction is interacting.
  {
    const Supernode &firstNode = _factor[ interaction1[ 0 ].first ];
    const SupernodeInteraction &firstInteraction
                     = firstNode.offDiagonal().at( interaction1[ 0 ].second );

    base_node_idx = firstInteraction._nodeID;
  }

  {
    const Supernode &secondNode = _factor[ interaction2[ 0 ].first ];
    const SupernodeInteraction &secondInteraction
                     = secondNode.offDiagonal().at( interaction2[ 0 ].second );

    interaction_node_idx = secondInteraction._nodeID;
  }

  // Get the full number of columns necessary to store in workspaces
  for ( int node_idx = 0; node_idx < nodeList.size(); node_idx++ )
  {
    fullColumnCount += _factor.at( nodeList[ node_idx ] ).numColumns();
  }

  // Clear the workspace
  MATRIX::clear( multWorkspaceInitial, fullColumnCount, nColsG );

  // Transpose multiply the matrix A1 with G
  sparseRowStart = _factor[ base_node_idx ].startColumn();

  // The complete matrix that we are multiplying with has dimension 
  // nRows x nColsG.  If we are transposing, then we should set this
  // to nRows; otherwise it is just the number of columns in the
  // relevant factor node
  if ( transpose ) {
    numSparseRows = nRows;
  } else {
    numSparseRows = _factor[ base_node_idx ].numColumns();
  }

  TIMING_STOP( "Interior block mult setup" );
#ifdef DO_TIMING
  _timers[ INTERIOR_BLOCK_MULT_SETUP ].tock();
#endif

  // Restrict to the given column range, if it is valid
  if ( columnRange.first >= 0 && columnRange.second >= 0 )
  {
    TRACE_ASSERT( rowRange.first >= 0 && rowRange.second >= rowRange.first );

    sparseRowStart += columnRange.first;

#if 0
    TRACE_ASSERT( columnRange.first + range_size( columnRange ) <= numSparseRows,
                  "Invalid row range" );
#endif

    numSparseRows = range_size( columnRange );

    useRowList = false;
  }

  TRACE_ASSERT( nRowsG == numSparseRows );

#ifdef DO_TIMING
  _timers[ INTERIOR_BLOCK_MULT_SPARSE ].tick();
#endif
  TIMING_START( "Interior block mult sparse" );
  if ( useRowList && transpose ) {
    sparseSubTransposeMultiply( A, nodeList,
                                sparseRowStart, numSparseRows, inverseRowList,
                                G, nColsG, multWorkspaceInitial, 1.0 );
  } else {
    sparseSubTransposeMultiply( A, nodeList, sparseRowStart, numSparseRows,
                                G, nColsG, multWorkspaceInitial, 1.0 );
  }
  TIMING_STOP( "Interior block mult sparse" );
#ifdef DO_TIMING
  _timers[ INTERIOR_BLOCK_MULT_SPARSE ].tock();
#endif

  // FIXME: debugging
  //
  // build our own workspace for now
#ifdef DO_TIMING
  _timers[ INTERIOR_BLOCK_MULT_WORKSPACE ].tick();
#endif
  TIMING_START( "Interior block build workspace" );
  MATRIX workspace( fullColumnCount, nColsG );
  TIMING_STOP( "Interior block build workspace" );
#ifdef DO_TIMING
  _timers[ INTERIOR_BLOCK_MULT_WORKSPACE ].tock();
#endif

  // Triangular solves
#ifdef DO_TIMING
  _timers[ INTERIOR_BLOCK_MULT_FORWARD_SOLVE ].tick();
#endif
  TIMING_START( "Interior block forward solve" );
  forwardSubSolve( nodeList, nColsG, multWorkspaceInitial, workspace.data() );
  TIMING_STOP( "Interior block forward solve" );
#ifdef DO_TIMING
  _timers[ INTERIOR_BLOCK_MULT_FORWARD_SOLVE ].tock();
#endif

#ifdef DO_TIMING
  _timers[ INTERIOR_BLOCK_MULT_BACKWARD_SOLVE ].tick();
#endif
  TIMING_START( "Interior block backward solve" );
  backwardSubSolve( nodeList, nColsG, multWorkspaceInitial, workspace.data() );
  TIMING_STOP( "Interior block backward solve" );
#ifdef DO_TIMING
  _timers[ INTERIOR_BLOCK_MULT_BACKWARD_SOLVE ].tock();
#endif

  // Standard multiply of matrix A2 with the result so far
  sparseRowStart = _factor[ interaction_node_idx ].startColumn();

  // If we are not transposing, then this is the matrix over which
  // we need to restrict the row set
  if ( transpose ) {
    numSparseRows = _factor[ interaction_node_idx ].numColumns();
  } else {
    numSparseRows = nRows;
  }

  // Restrict to the given row range, if it is valid
  if ( rowRange.first >= 0 && rowRange.second >= 0 )
  {
    sparseRowStart += rowRange.first;

#if 0
    TRACE_ASSERT( rowRange.first + range_size( rowRange ) <= numSparseRows,
                  "Invalid row range" );
#endif

    numSparseRows = range_size( rowRange );
  }

  // Next, subtract directly from the output data matrix
#ifdef DO_TIMING
  _timers[ INTERIOR_BLOCK_MULT_SPARSE ].tick();
#endif
  TIMING_START( "Interior block mult sparse" );
  if ( !transpose && useRowList ) {
    sparseSubMultiply( A, nodeList,
                       sparseRowStart, numSparseRows, inverseRowList,
                       multWorkspaceInitial, nColsG, outputData,
                       -1.0 /* subtract */ );
  } else {
    sparseSubMultiply( A, nodeList, sparseRowStart, numSparseRows,
                       multWorkspaceInitial, nColsG, outputData,
                       -1.0 /* subtract */ );
  }
  TIMING_STOP( "Interior block mult sparse" );
#ifdef DO_TIMING
  _timers[ INTERIOR_BLOCK_MULT_SPARSE ].tock();
#endif
}
#endif

//////////////////////////////////////////////////////////////////////
// New version to work on
//////////////////////////////////////////////////////////////////////
void FactorManager::interiorBlockSchurMultiply(
                            const SPARSE_MATRIX::SparseColumnMatrix &A,
                            IntArray &nodeList,
                            int interior_block_idx,
                            int interaction_idx1, int interaction_idx2,
                            const Real *G,
                            int nRowsG, int nColsG,
                            const IntArray &inverseRowList, int nRows,
                            Real *multWorkspaceInitial,
                            Real *multWorkspaceIntermediate,
                            Real *outputData,
                            IndexRange rowRange,
                            IndexRange columnRange,
                            bool transpose, bool left )
{
  int                        fullColumnCount = 0;

  // Get the full number of columns necessary to store in workspaces
  for ( int node_idx = 0; node_idx < nodeList.size(); node_idx++ )
  {
    fullColumnCount += _factor.at( nodeList[ node_idx ] ).numColumns();
  }

  MATRIX                     workspace( fullColumnCount, nColsG );

  InteriorBlockMultType      multType = MULT_FULL;

  if ( rowRange.first >= 0 ) {
    TRACE_ASSERT( rowRange.second >= rowRange.first );
    TRACE_ASSERT( columnRange.first >= 0 );
    TRACE_ASSERT( columnRange.second >= columnRange.first );

    multType = MULT_ROW_RANGE;
  }

  if ( transpose ) {
#if 0
    nodeList.clear();
    buildInteriorBlockNodeList( interior_block_idx,
                                interaction_idx2, interaction_idx2, nodeList );

    // Multiply; Interaction1 * Interaction2' * G
    interiorBlockInteractionMultiply( A, nodeList,
                                      interior_block_idx, interaction_idx2,
                                      G, nRowsG, nColsG,
                                      multType == MULT_ROW_RANGE
                                        ? MULT_ROW_RANGE : MULT_ROW_SET,
                                      inverseRowList, nRows,
                                      NULL,
                                      workspace.data(),
                                      multWorkspaceInitial,
                                      1.0 /* alpha */,
                                      rowRange,
                                      true );

    // Copy
    transformInteriorBlockWorkspace( interior_block_idx,
                                     interaction_idx2, interaction_idx1,
                                     nColsG, multWorkspaceInitial,
                                     multWorkspaceIntermediate );

    // new node list
    nodeList.clear();
    buildInteriorBlockNodeList( interior_block_idx,
                                interaction_idx1, interaction_idx1, nodeList );

    interiorBlockInteractionMultiply( A, nodeList,
                                      interior_block_idx, interaction_idx1,
                                      multWorkspaceIntermediate,
                                      fullColumnCount, nColsG,
                                      multType == MULT_ROW_RANGE
                                        ? MULT_ROW_RANGE : MULT_FULL,
                                      inverseRowList, nRows,
                                      multWorkspaceInitial,
                                      workspace.data(),
                                      outputData,
                                      -1.0 /* alpha */,
                                      columnRange,
                                      false );
#endif
    interiorBlockInteractionMultiply( A, nodeList,
                                      interior_block_idx, interaction_idx2,
                                      G, nRowsG, nColsG,
                                      multType == MULT_ROW_RANGE
                                        ? MULT_ROW_RANGE : MULT_ROW_SET,
                                      inverseRowList, nRows,
                                      NULL,
                                      workspace.data(),
                                      multWorkspaceIntermediate,
                                      1.0 /* alpha */,
                                      rowRange,
                                      true );

    interiorBlockInteractionMultiply( A, nodeList,
                                      interior_block_idx, interaction_idx1,
                                      multWorkspaceIntermediate,
                                      fullColumnCount, nColsG,
                                      multType == MULT_ROW_RANGE
                                        ? MULT_ROW_RANGE : MULT_FULL,
                                      inverseRowList, nRows,
                                      multWorkspaceInitial,
                                      workspace.data(),
                                      outputData,
                                      -1.0 /* alpha */,
                                      columnRange,
                                      false );
  } else {
    // Multiply; Interaction2 * Interaction1' * G
    interiorBlockInteractionMultiply( A, nodeList,
                                      interior_block_idx, interaction_idx1,
                                      G, nRowsG, nColsG,
                                      multType == MULT_ROW_RANGE
                                        ? MULT_ROW_RANGE : MULT_FULL,
                                      inverseRowList, nRows,
                                      NULL,
                                      workspace.data(),
                                      multWorkspaceIntermediate,
                                      1.0 /* alpha */,
                                      columnRange,
                                      true );

    interiorBlockInteractionMultiply( A, nodeList,
                                      interior_block_idx, interaction_idx2,
                                      multWorkspaceIntermediate,
                                      fullColumnCount, nColsG,
                                      multType == MULT_ROW_RANGE
                                        ? MULT_ROW_RANGE : MULT_ROW_SET,
                                      inverseRowList, nRows,
                                      multWorkspaceInitial,
                                      workspace.data(),
                                      outputData,
                                      -1.0 /* alpha */,
                                      rowRange,
                                      false );
  }
}

//////////////////////////////////////////////////////////////////////
// Same as the function above, but assumes that we are multiplying
// with all interactions in the given block below ancestor_idx.
//
// Note: we don't need any row/column range stuff here, because
// this will not be used for decomposing diagonal blocks
//////////////////////////////////////////////////////////////////////
void FactorManager::interiorBlockSchurMultiply(
                            const SPARSE_MATRIX::SparseColumnMatrix &A,
                            IntArray &nodeList,
                            int interior_block_idx,
                            int ancestor_idx,
                            const Real *G,
                            int nRowsG, int nColsG,
                            const IntArray &inverseRowList, int nRows,
                            Real *multWorkspaceInitial,
                            Real *multWorkspaceIntermediate,
                            Real *outputData,
                            int sparseRowStartOverride,
                            bool transpose )
{
  int                        fullColumnCount = 0;

  // Get the full number of columns necessary to store in workspaces
  for ( int node_idx = 0; node_idx < nodeList.size(); node_idx++ ) {
    fullColumnCount += _factor.at( nodeList[ node_idx ] ).numColumns();
  }

  MATRIX                     workspace( fullColumnCount, nColsG );

  // Note that we are assuming here that inverseRowList has been built to
  // include all rows starting after ancestor_idx
  if ( transpose ) {
    interiorBlockInteractionMultiply( A, nodeList,
                                      interior_block_idx, ancestor_idx + 1,
                                      G, nRowsG, nColsG,
                                      MULT_ROW_SET,
                                      inverseRowList, nRows,
                                      NULL,
                                      workspace.data(),
                                      multWorkspaceIntermediate,
                                      1.0 /* alpha */,
                                      IndexRange( -1, -1 ) /* placeholder */,
                                      true,
                                      sparseRowStartOverride );

    interiorBlockInteractionMultiply( A, nodeList,
                                      interior_block_idx, ancestor_idx,
                                      multWorkspaceIntermediate,
                                      fullColumnCount, nColsG,
                                      MULT_FULL,
                                      inverseRowList, nRows,
                                      multWorkspaceInitial,
                                      workspace.data(),
                                      outputData,
                                      -1.0 /* alpha */,
                                      IndexRange( -1, -1 ) /* placeholder */,
                                      false );
  } else {
    // Multiply; Interaction2 * Interaction1' * G
    interiorBlockInteractionMultiply( A, nodeList,
                                      interior_block_idx, ancestor_idx,
                                      G, nRowsG, nColsG,
                                      MULT_FULL,
                                      inverseRowList, nRows,
                                      NULL,
                                      workspace.data(),
                                      multWorkspaceIntermediate,
                                      1.0 /* alpha */,
                                      IndexRange( -1, -1 ) /* placeholder */,
                                      true );

    interiorBlockInteractionMultiply( A, nodeList,
                                      interior_block_idx, ancestor_idx + 1,
                                      multWorkspaceIntermediate,
                                      fullColumnCount, nColsG,
                                      MULT_ROW_SET,
                                      inverseRowList, nRows,
                                      multWorkspaceInitial,
                                      workspace.data(),
                                      outputData,
                                      -1.0 /* alpha */,
                                      IndexRange( -1, -1 ) /* placeholder */,
                                      false,
                                      sparseRowStartOverride );
  }
}

//////////////////////////////////////////////////////////////////////
// Helper function for the above
//
//
//////////////////////////////////////////////////////////////////////
void FactorManager::transformInteriorBlockWorkspace( int interior_block_idx,
                                                     int interaction_idx1,
                                                     int interaction_idx2,
                                                     int nCols,
                                                     const Real *inputWorkspace,
                                                     Real *outputWorkspace )
{
  const InteriorBlock       &block = _interiorBlocks[ interior_block_idx ];

  const InteriorBlock::Interaction  &interaction1
                                    = block._interactions[ interaction_idx1 ];
  const InteriorBlock::Interaction  &interaction2
                                    = block._interactions[ interaction_idx2 ];

  int                        nRows1 = 0;
  int                        nRows2 = 0;

  int                        idx1 = 0;
  int                        idx2 = 0;

  for ( int i = 0; i < interaction1.size(); i++ ) {
    nRows1 += _factor[ interaction1[ i ].first ].numColumns();
  }

  for ( int i = 0; i < interaction2.size(); i++ ) {
    nRows2 += _factor[ interaction2[ i ].first ].numColumns();
  }

  while ( nRows1 > 0 && nRows2 > 0 ) {
    int                      node_id1 = interaction1[ idx1 ].first;
    int                      node_id2 = interaction2[ idx2 ].first;

    if ( node_id1 < node_id2 ) {
      // Just skip this block
      nRows1 -= _factor[ node_id1 ].numColumns();
      inputWorkspace += _factor[ node_id1 ].numColumns() * nCols;
      idx1++;
    } else if ( node_id1 == node_id2 ) {
      // Copy the block
      MATRIX::copy( outputWorkspace, inputWorkspace,
                    _factor[ node_id1 ].numColumns(), nCols );

      nRows1 -= _factor[ node_id1 ].numColumns();
      nRows2 -= _factor[ node_id1 ].numColumns();
      inputWorkspace += _factor[ node_id1 ].numColumns() * nCols;
      outputWorkspace += _factor[ node_id1 ].numColumns() * nCols;
      idx1++;
      idx2++;
    } else {
      // Add a zero block to the output workspace
      MATRIX::clear( outputWorkspace,
                     _factor[ node_id1 ].numColumns(), nCols );

      nRows2 -= _factor[ node_id2 ].numColumns();
      outputWorkspace += _factor[ node_id2 ].numColumns() * nCols;
      idx2++;
    }
  }

  if ( nRows2 > 0 ) {
    MATRIX::clear( outputWorkspace, nRows2, nCols );
  }
}

//////////////////////////////////////////////////////////////////////
// Performs multiplication with a single row block from an interior block
// in the factor
//////////////////////////////////////////////////////////////////////
void FactorManager::interiorBlockInteractionMultiply(
                            const SPARSE_MATRIX::SparseColumnMatrix &A,
                            const IntArray &nodeList,
                            int interior_block_idx, int interaction_idx,
                            const Real *G,
                            int nRowsG, int nColsG,
                            InteriorBlockMultType multType,
                            // In case we only want to use a certain
                            // subset of rows when doing sparse multiplication
                            // (ie. if multType == MULT_ROW_SET)
                            const IntArray &inverseRowList, int nRows,
                            Real *multWorkspaceInitial,
                            Real *solveWorkspace,
                            Real *outputData,
                            // Factor to multiply by when adding to
                            // the output
                            Real alpha,
                            // Optional row range in the sparse matrix
                            // (used if multType == MULT_ROW_RANGE)
                            IndexRange rowRange,
                            bool transpose,
                            int sparseRowStartOverride )
{
  TIMING_START( "interiorBlockInteractionMultiply preamble" );
  int                        fullColumnCount = 0;

  int                        sparseRowStart, sparseColStart;
  int                        numSparseRows, numSparseCols;

  int                        base_node_idx;

  int                        node_idx;

  const InteriorBlock       &block = _interiorBlocks[ interior_block_idx ];

  const InteriorBlock::Interaction  &interaction
                                      = block._interactions[ interaction_idx ];

  {
    const Supernode &interactionNode = _factor[ interaction[ 0 ].first ];
    const SupernodeInteraction &firstInteraction
                  = interactionNode.offDiagonal().at( interaction[ 0 ].second );

    base_node_idx = firstInteraction._nodeID;
  }

  // Get the full number of columns necessary to store in workspaces
  for ( int node_idx = 0; node_idx < nodeList.size(); node_idx++ ) {
    fullColumnCount += _factor.at( nodeList[ node_idx ] ).numColumns();
  }

  sparseRowStart = sparseRowStartOverride > 0 ?
                              sparseRowStartOverride
                            : _factor[ base_node_idx ].startColumn();

#if 0
  // FIXME: debugging
  if ( sparseRowStartOverride > 0 ) {
    printf( "Overriding sparseRowStart with %d\n", sparseRowStartOverride );
  }
#endif

  switch ( multType ) {
    case MULT_FULL:
    {
      numSparseRows = _factor[ base_node_idx ].numColumns();
      break;
    }
    case MULT_ROW_SET:
    {
      numSparseRows = nRows;
      break;
    }
    case MULT_ROW_RANGE:
    default:
    {
      sparseRowStart += rowRange.first;
      numSparseRows = range_size( rowRange );
      break;
    }
  }

  TIMING_STOP( "interiorBlockInteractionMultiply preamble" );

  if ( transpose ) {
    TRACE_ASSERT( numSparseRows == nRowsG );

    TIMING_START( "interiorBlockInteractionMultiply sparseTransMult" );
    // If we are multiplying by the transpose, we start by applying the
    // transpose of the desired spase matrix and accumulating the result
    // in multWorkspaceInitial
    MATRIX::clear( outputData, fullColumnCount, nColsG );

    if ( multType == MULT_ROW_SET ) {
      sparseSubTransposeMultiply( A, nodeList,
                                  sparseRowStart, numSparseRows, inverseRowList,
                                  G, nColsG, outputData, 1.0 );
    } else {
      sparseSubTransposeMultiply( A, nodeList,
                                  sparseRowStart, numSparseRows,
                                  G, nColsG, outputData, 1.0 );
    }
    TIMING_STOP( "interiorBlockInteractionMultiply sparseTransMult" );

    // Follow this with a forward solve
    TIMING_START( "interiorBlockInteractionMultiply forwardSolve" );
    forwardSubSolve( nodeList, nColsG, outputData, solveWorkspace );
    TIMING_STOP( "interiorBlockInteractionMultiply forwardSolve" );
  } else {
    // If we are multiplying by the untransposed matrix, we start by
    // doing a backward solve
    TIMING_START( "interiorBlockInteractionMultiply backwardSolve" );
    MATRIX::copy( multWorkspaceInitial, G, fullColumnCount, nColsG );

    backwardSubSolve( nodeList, nColsG, multWorkspaceInitial,
                      solveWorkspace );
    TIMING_STOP( "interiorBlockInteractionMultiply backwardSolve" );

    // Now, perform an untransposed sparse matrix multiply and add to the
    // output using alpha
    TIMING_START( "interiorBlockInteractionMultiply sparseMult" );
    if ( multType == MULT_ROW_SET ) {
      sparseSubMultiply( A, nodeList,
                         sparseRowStart, numSparseRows, inverseRowList,
                         multWorkspaceInitial, nColsG, outputData,
                         alpha );
    } else {
      sparseSubMultiply( A, nodeList,
                         sparseRowStart, numSparseRows,
                         multWorkspaceInitial, nColsG, outputData,
                         alpha );
    }
    TIMING_STOP( "interiorBlockInteractionMultiply sparseMult" );
  }
}

#if 0
//////////////////////////////////////////////////////////////////////
// Similar to the above function, but multiplies with all interactions
// *BELOW* ancestor_idx.
//
// Note: We don't need the rowRange parameter here because this
// will not be used for diagonal block decomposition
//////////////////////////////////////////////////////////////////////
void FactorManager::interiorBlockInteractionMultiply_allRows(
                            const SPARSE_MATRIX::SparseColumnMatrix &A,
                            const IntArray &nodeList,
                            int interior_block_idx, int ancestor_idx,
                            const Real *G,
                            int nRowsG, int nColsG,
                            InteriorBlockMultType multType,
                            // In case we only want to use a certain
                            // subset of rows when doing sparse multiplication
                            // (ie. if multType == MULT_ROW_SET)
                            const IntArray &inverseRowList, int nRows,
                            Real *multWorkspaceInitial,
                            Real *solveWorkspace,
                            Real *outputData,
                            // Factor to multiply by when adding to
                            // the output.  Note; alpha is ignored for
                            // transposed multiplications
                            Real alpha,
                            bool transpose )
{
  TIMING_START( "interiorBlockInteractionMultiply_allRows preamble" );
  int                        fullColumnCount = 0;

  int                        sparseRowStart, sparseColStart;
  int                        numSparseRows, numSparseCols;

  int                        base_node_idx;

  int                        node_idx;

  const InteriorBlock       &block = _interiorBlocks[ interior_block_idx ];

  const InteriorBlock::Interaction  &interaction
                                      = block._interactions[ interaction_idx ];

  {
    const Supernode &interactionNode = _factor[ interaction[ 0 ].first ];
    const SupernodeInteraction &firstInteraction
                  = interactionNode.offDiagonal().at( interaction[ 0 ].second );

    base_node_idx = firstInteraction._nodeID;
  }

  // Get the full number of columns necessary to store in workspaces
  for ( int node_idx = 0; node_idx < nodeList.size(); node_idx++ )
  {
    fullColumnCount += _factor.at( nodeList[ node_idx ] ).numColumns();
  }

  sparseRowStart = _factor[ base_node_idx ].startColumn();

  switch ( multType ) {
    case MULT_FULL:
    {
      numSparseRows = _factor[ base_node_idx ].numColumns();
      break;
    }
    case MULT_ROW_SET:
    {
      numSparseRows = nRows;
      break;
    }
    case MULT_ROW_RANGE:
    default:
    {
      TRACE_ASSERT( NULL, "Invalid multiplication type" );
      break;
    }
  }

  TIMING_STOP( "interiorBlockInteractionMultiply_allRows preamble" );

  if ( transpose ) {
    TRACE_ASSERT( numSparseRows == nRowsG );

    TIMING_START( "interiorBlockInteractionMultiply_allRows sparseTransMult" );
    // If we are multiplying by the transpose, we start by applying the
    // transpose of the desired sparse matrix and accumulating the result
    // in multWorkspaceInitial
    MATRIX::clear( outputData, fullColumnCount, nColsG );

    if ( multType == MULT_ROW_SET ) {
      sparseSubTransposeMultiply( A, nodeList,
                                  sparseRowStart, numSparseRows, inverseRowList,
                                  G, nColsG, outputData, 1.0 );
    } else { 
      sparseSubTransposeMultiply( A, nodeList,
                                  sparseRowStart, numSparseRows,
                                  G, nColsG, outputData, 1.0 );
    }

    TIMING_STOP( "interiorBlockInteractionMultiply_allRows sparseTransMult" );

    // Follow this with a forward solve
    TIMING_START( "interiorBlockInteractionMultiply_allRows forwardSolve" );
    forwardSubSolve( nodeList, nColsG, outputData, solveWorkspace );
    TIMING_STOP( "interiorBlockInteractionMultiply_allRows forwardSolve" );
  } else {
    // If we are multiplying by the untransposed matrix, we start by
    // doing a backward solve
    TIMING_START( "interiorBlockInteractionMultiply_allRows backwardSolve" );
    MATRIX::copy( multWorkspaceInitial, G, fullColumnCount, nColsG );

    backwardSubSolve( nodeList, nColsG, multWorkspaceInitial,
                      solveWorkspace );
    TIMING_STOP( "interiorBlockInteractionMultiply_allRows backwardSolve" );

    // Now perform an untransposed sparse matrix multiply and add to the
    // output using alpha
    TIMING_START( "interiorBlockInteractionMultiply sparseMult" );
    if ( multType == MULT_ROW_SET ) {
      sparseSubMultiply( A, nodeList,
                         sparseRowStart, numSparseRows, inverseRowList,
                         multWorkspaceInitial, nColsG, outputData,
                         alpha );
    } else {
      sparseSubMultiply( A, nodeList,
                         sparseRowStart, numSparseRows,
                         multWorkspaceInitial, nColsG, outputData,
                         alpha );
    }
    TIMING_STOP( "interiorBlockInteractionMultiply sparseMult" );
  }
}
#endif

//////////////////////////////////////////////////////////////////////
// Builds a node list for the function above
//////////////////////////////////////////////////////////////////////
void FactorManager::buildInteriorBlockNodeList( int interior_block_idx,
                                                int interaction_idx1,
                                                int interaction_idx2,
                                                IntArray &nodeList )
{
  const InteriorBlock        &block = _interiorBlocks[ interior_block_idx ];

  const InteriorBlock::Interaction  &interaction1
                                      = block._interactions[ interaction_idx1 ];
  const InteriorBlock::Interaction  &interaction2
                                      = block._interactions[ interaction_idx2 ];

  int                        idx1 = 0;
  int                        idx2 = 0;

  nodeList.clear();

  // Get the union of the two sorted lists
  while ( idx1 < interaction1.size() || idx2 < interaction2.size() )
  {
    if ( idx1 < interaction1.size() && idx2 < interaction2.size()
      && interaction1[ idx1 ].first <= interaction2[ idx2 ].first )
    {
      nodeList.push_back( interaction1[ idx1 ].first );

      if ( interaction1[ idx1 ].first == interaction2[ idx2 ].first )
      {
        idx2 += 1;
      }

      idx1 += 1;
    }
    else if ( idx1 < interaction1.size() && idx2 < interaction2.size()
           && interaction1[ idx1 ].first > interaction2[ idx2 ].first )
    {
      nodeList.push_back( interaction2[ idx2 ].first );

      idx2 += 1;
    }
    else if ( idx1 < interaction1.size() )
    {
      nodeList.push_back( interaction1[ idx1 ].first );

      idx1 += 1;
    }
    else
    {
      nodeList.push_back( interaction2[ idx2 ].first );

      idx2 += 1;
    }
  }
#if 0
  // FIXME: try taking the intersection instead
  while ( idx1 < interaction1.size() && idx2 < interaction2.size() )
  {
    if ( interaction1[ idx1 ].first < interaction2[ idx2 ].first )
    {
      idx1++;
    }
    else if ( interaction2[ idx2 ].first < interaction1[ idx1 ].first )
    {
      idx2++;
    }
    else
    {
      TRACE_ASSERT( interaction1[ idx1 ].first == interaction2[ idx2 ].first );

      nodeList.push_back( interaction1[ idx1 ].first );

      idx1++;
      idx2++;
    }
  }
#endif
#if 0
  for ( idx2 = 0; idx2 < interaction2.size(); idx2++ )
  {
    nodeList.push_back( interaction2[ idx2 ].first );
  }
#endif
}

//////////////////////////////////////////////////////////////////////
// Builds a node list for interior block multiplies assuming
// that we are actually multiplying with everything "below"
// ancestor_idx
//////////////////////////////////////////////////////////////////////
void FactorManager::buildInteriorBlockNodeList( int interior_block_idx,
                                                int ancestor_idx,
                                                IntArray &nodeList )
{
  const InteriorBlock       &block = _interiorBlocks[ interior_block_idx ];

  // FIXME: this may be kind of inefficient, but maybe not
  OrderedSet<int>::type      nodeSet;

  // Build the union set of all nodes involved in interactions
  // below ancestor_idx
  for ( int interaction_idx = ancestor_idx;
        interaction_idx < block._interactions.size(); interaction_idx++ )
  {
    const InteriorBlock::Interaction  &interaction
                                      = block._interactions[ interaction_idx ];

    for ( int i = 0; i < interaction.size(); i++ ) {
      nodeSet.insert( interaction[ i ].first );
    }
  }

  nodeList.clear();
  for ( OrderedSet<int>::type::const_iterator i = nodeSet.begin();
        i != nodeSet.end(); i++ )
  {
    nodeList.push_back( *i );
  }
#if 0
  for ( int node_idx = block._nodeRange.first;
        node_idx <= block._nodeRange.second; node_idx++ )
  {
    nodeList.push_back( node_idx );
  }
#endif
}

//////////////////////////////////////////////////////////////////////
// Multiplies the transpose of the sparse matrix over the given
// row range, with the set of node indices specifying the columns to use
//////////////////////////////////////////////////////////////////////
void FactorManager::sparseSubMultiply(
                        const SPARSE_MATRIX::SparseColumnMatrix &A,
                        const IntArray &nodeList,
                        int sparseRowStart, int numSparseRows,
                        const Real *G, int nColsG,
                        Real *outputMatrix,
                        Real alpha )
{
  const Real                *baseInputData;

  int                        sparseColStart;
  int                        numSparseCols;

  baseInputData = G;

  for ( int node_idx = 0; node_idx < nodeList.size(); node_idx++ )
  {
    sparseColStart = _factor.at( nodeList[ node_idx ] ).startColumn();
    numSparseCols = _factor.at( nodeList[ node_idx ] ).numColumns();

    // Standard (non-transposed) multiply using this sub-block of the
    // sparse matrix
    SPARSE_MATRIX::subMatrixMultiply( A, baseInputData, outputMatrix,
                                      sparseRowStart, sparseColStart,
                                      numSparseRows, numSparseCols, nColsG,
                                      alpha,
                                      false /* Don't clear */,
                                      false /* Don't transpose */ );
  
    // Move the input data to the next block
    baseInputData += numSparseCols * nColsG;
  }
}

//////////////////////////////////////////////////////////////////////
// Only multiplies using a subset of the matrix rows using the provided
// inverse row map
//////////////////////////////////////////////////////////////////////
void FactorManager::sparseSubMultiply(
                                const SPARSE_MATRIX::SparseColumnMatrix &A,
                                const IntArray &nodeList,
                                int sparseRowStart, int numSparseRows,
                                const IntArray &inverseRowList,
                                const Real *G, int nColsG,
                                Real *outputMatrix,
                                Real alpha )
{
  const Real                *baseInputData;

  int                        sparseColStart;
  int                        numSparseCols;

  baseInputData = G;

#if 0
  cout << "sparseSubMultiply ";
  cout << SDUMP( numSparseRows );
  cout << SDUMP( inverseRowList.size() ) << endl;
#endif

  for ( int node_idx = 0; node_idx < nodeList.size(); node_idx++ )
  {
    sparseColStart = _factor.at( nodeList[ node_idx ] ).startColumn();
    numSparseCols = _factor.at( nodeList[ node_idx ] ).numColumns();

    // Standard (non-transposed) multiply using this sub-block of the
    // sparse matrix
    SPARSE_MATRIX::subMatrixMultiply( A, baseInputData, outputMatrix,
                                      sparseRowStart, sparseColStart,
                                      numSparseRows, numSparseCols, nColsG,
                                      inverseRowList,
                                      alpha,
                                      false /* Don't clear */,
                                      false /* Don't transpose */ );
  
    // Move the input data to the next block
    baseInputData += numSparseCols * nColsG;
  }
}

#if 0
//////////////////////////////////////////////////////////////////////
// Same as the above, but G is applied on the left
//////////////////////////////////////////////////////////////////////
void FactorManager::sparseSubLeftMultiply(
                                const SPARSE_MATRIX::SparseColumnMatrix &A,
                                const IntArray &nodeList,
                                int fullColumnCount,
                                int sparseRowStart, int numSparseRows,
                                const Real *G, int nRowsG,
                                Real *outputMatrix )
{
  Real                      *baseOutputData;

  int                        sparseColStart;
  int                        numSparseCols;

  baseOutputData = outputMatrix;

  for ( int node_idx = 0; node_idx < nodeList.size(); node_idx++ )
  {
    sparseColStart = _factor.at( nodeList[ node_idx ] ).startColumn();
    numSparseCols = _factor.at( nodeList[ node_idx ] ).numColumns();

    // Standard (non-transposed) multiply using this sub-block of the
    // sparse matrix
    SPARSE_MATRIX::subMatrixLeftMultiply( A, G, baseOutputData,
                                          sparseRowStart, sparseColStart,
                                          numSparseRows, numSparseCols, nRowsG,
                                          false /* Don't clear */,
                                          false /* Don't transpose */,
                                          -1 /* Don't specify G's leading
                                                dimension */,
                                          fullColumnCount /* Need a leading
                                                             dimension for
                                                             output */ );
  
    // Move to the next column in the output
    baseOutputData += numSparseCols;
  }
}
#endif

//////////////////////////////////////////////////////////////////////
// Multiplies the sparse matrix over the given row range, with the
// set of node indices specifying the columns to use
//
// This is a helper function for interiorBlockSchurMultiply
//////////////////////////////////////////////////////////////////////
void FactorManager::sparseSubTransposeMultiply(
                                const SPARSE_MATRIX::SparseColumnMatrix &A,
                                const IntArray &nodeList,
                                int sparseRowStart, int numSparseRows,
                                const Real *G, int nColsG,
                                Real *outputMatrix,
                                Real alpha )
{
  Real                      *baseOutputData;

  int                        sparseColStart;
  int                        numSparseCols;

  // FIXME: debugging
  int                        totalNNZ = 0;

  baseOutputData = outputMatrix;

  for ( int node_idx = 0; node_idx < nodeList.size(); node_idx++ )
  {
    sparseColStart = _factor.at( nodeList[ node_idx ] ).startColumn();
    numSparseCols = _factor.at( nodeList[ node_idx ] ).numColumns();

    // Transpose multiply using this sub-block of the sparse matrix
    totalNNZ += SPARSE_MATRIX::subMatrixMultiply(
                                      A, G, baseOutputData,
                                      sparseRowStart, sparseColStart,
                                      numSparseRows, numSparseCols, nColsG,
                                      alpha,
                                      false /* Don't need to clear */,
                                      true /* Transpose the block of A */ );

    // Move the output data to the next block
    baseOutputData += numSparseCols * nColsG;
  }
}

//////////////////////////////////////////////////////////////////////
// Only multiplies using a subset of the matrix rows using the
// provided inverse row map
//////////////////////////////////////////////////////////////////////
void FactorManager::sparseSubTransposeMultiply(
                                const SPARSE_MATRIX::SparseColumnMatrix &A,
                                const IntArray &nodeList,
                                int sparseRowStart, int numSparseRows,
                                const IntArray &inverseRowList,
                                const Real *G, int nColsG,
                                Real *outputMatrix,
                                Real alpha )
{
  Real                      *baseOutputData;

  int                        sparseColStart;
  int                        numSparseCols;

  // FIXME: debugging
  int                        totalNNZ = 0;

  baseOutputData = outputMatrix;

#if 0
  cout << "sparseSubTransposeMultiply ";
  cout << SDUMP( numSparseRows );
  cout << SDUMP( inverseRowList.size() ) << endl;
#endif

  for ( int node_idx = 0; node_idx < nodeList.size(); node_idx++ )
  {
    sparseColStart = _factor.at( nodeList[ node_idx ] ).startColumn();
    numSparseCols = _factor.at( nodeList[ node_idx ] ).numColumns();

    // Transpose multiply using this sub-block of the sparse matrix
    totalNNZ += SPARSE_MATRIX::subMatrixMultiply(
                                      A, G, baseOutputData,
                                      sparseRowStart, sparseColStart,
                                      numSparseRows, numSparseCols, nColsG,
                                      inverseRowList,
                                      alpha,
                                      false /* Don't need to clear */,
                                      true /* Transpose the block of A */ );

    // Move the output data to the next block
    baseOutputData += numSparseCols * nColsG;
  }
}

#if 0
//////////////////////////////////////////////////////////////////////
// Same as the above, but G is applied on the left
//////////////////////////////////////////////////////////////////////
void FactorManager::sparseSubTransposeLeftMultiply(
                                const SPARSE_MATRIX::SparseColumnMatrix &A,
                                const IntArray &nodeList,
                                int fullColumnCount,
                                int sparseRowStart, int numSparseRows,
                                const Real *G, int nRowsG,
                                Real *outputMatrix )
{
  const Real                *baseInputData;

  int                        sparseColStart;
  int                        numSparseCols;

  baseInputData = G;

  for ( int node_idx = 0; node_idx < nodeList.size(); node_idx++ )
  {
    sparseColStart = _factor.at( nodeList[ node_idx ] ).startColumn();
    numSparseCols = _factor.at( nodeList[ node_idx ] ).numColumns();

    SPARSE_MATRIX::subMatrixLeftMultiply( A, baseInputData, outputMatrix,
                                          sparseRowStart, sparseColStart,
                                          numSparseRows, numSparseCols, nRowsG,
                                          false /* Don't need to clear */,
                                          true /* Transpose the block of A */,
                                          fullColumnCount /* Leading dimension
                                                             of the input */ );

    // Move the input to the next column
    G += numSparseCols;
  }
}
#endif

//////////////////////////////////////////////////////////////////////
// Given the list of current descendents for a node, constructs
// a list of interior blocks to which those descendents belong
//////////////////////////////////////////////////////////////////////
void FactorManager::findInteriorBlockDescendents()
{
  int                        descendent_idx;
  int                        interior_block_idx;

  // Clear anything left over from previous iterations
  clearInteriorBlockDescendents();

  // Look at all individual node descendents and see which interior
  // block they belong to
  for ( int i = 0; i < _currentDescendents.size(); i++ )
  {
    descendent_idx = _currentDescendents[ i ];
    interior_block_idx = _interiorBlockMap[ descendent_idx ];

    if ( interior_block_idx == EMPTY ) {
      continue;
    }

    TRACE_ASSERT( interior_block_idx >= 0,
                  "Descendent does not belong to an interior block" );

    // If we haven't already found this interior block in the
    // descendent list, add it to our block descendent list and flag
    // it as having been found.
    if ( !_blockDescendentFound[ interior_block_idx ] )
    {
      _currentBlockDescendents.push_back( interior_block_idx );
      _blockDescendentFound[ interior_block_idx ] = true;
    }
  }

#if 0
  printf( "Found %d interior block descendents\n",
          (int)_currentBlockDescendents.size() );
#endif
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void FactorManager::clearInteriorBlockDescendents()
{
  for ( int i = 0; i < _currentBlockDescendents.size(); i++ ) {
    _blockDescendentFound.at( _currentBlockDescendents[ i ] ) = false;
  }

  _currentBlockDescendents.clear();
}
