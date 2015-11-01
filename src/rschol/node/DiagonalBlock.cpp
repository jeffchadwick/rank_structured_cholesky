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
// DiagonalBlock.cpp: Implementation
//
//////////////////////////////////////////////////////////////////////

#include "DiagonalBlock.h"

#include <rschol/node/Supernode.h>

#include <rschol/util/MathUtil.h>

#include <map>
#include <stack>

using namespace std;

//////////////////////////////////////////////////////////////////////
// Constructor
//////////////////////////////////////////////////////////////////////
DiagonalBlock::DiagonalBlock()
  : _node( NULL ),
    _root( NULL ),
    _explicitBlocks( false )
{
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
DiagonalBlock::DiagonalBlock( const DiagonalBlock &block )
  : _node( NULL ),
    _root( NULL )
{
  copy( block );
}

//////////////////////////////////////////////////////////////////////
// Destructor
//////////////////////////////////////////////////////////////////////
DiagonalBlock::~DiagonalBlock()
{
  clear();
#if 0
  delete _root;
#endif
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
DiagonalBlock &DiagonalBlock::operator=( const DiagonalBlock &block )
{
  clear();
  copy( block );
  return *this;
}

//////////////////////////////////////////////////////////////////////
// Initialize the block
//////////////////////////////////////////////////////////////////////
void DiagonalBlock::init( const vector<DenseBlock> &diagonalBlocks,
                          const vector<long int> &diagonalOffsets,
                          int explicitBlockThreshold )
{
  vector<BlockNode *>        currentLevel( diagonalBlocks.size() );

  // Build the bottom level of the tree
  for ( int block_idx = 0; block_idx < diagonalBlocks.size(); block_idx++ ) {
    currentLevel[ block_idx ] = new BlockNode(
                                    diagonalBlocks[ block_idx ]._rowRange,
                                    diagonalBlocks[ block_idx ]._columnRange );
  }

  // Assemble the tree from the bottom up
  while ( currentLevel.size() > 1 ) {
    TRACE_ASSERT( currentLevel.size() % 2 == 0,
                  "Number of diagonal blocks must be a power of two" );

    vector<BlockNode *>      nextLevel( currentLevel.size() / 2 );

    for ( int block_idx = 0; block_idx < nextLevel.size(); block_idx++ ) {
      // Assemble a new block out of two adjacent nodes
      nextLevel[ block_idx ] = new BlockNode();
      nextLevel[ block_idx ]->setChildren( currentLevel[ block_idx * 2 + 0 ],
                                           currentLevel[ block_idx * 2 + 1 ] );
    }

    currentLevel = nextLevel;
  }

  if ( currentLevel.size() != 1 ) {
    cout << "Current level size = " << currentLevel.size() << endl;
    cout << "Diagonal blocks size = " << diagonalBlocks.size() << endl;
  }
  TRACE_ASSERT( currentLevel.size() == 1 );

  _root = currentLevel[ 0 ];

  TRACE_ASSERT( _root->_parent == NULL );

  _nodeTraversal.clear();
  buildTraversal( _root );

  // Build lists of diagonal and off-diagonal blocks
  _diagonalBlocks.clear();
  _offDiagonalBlocks.clear();
  for ( int block_idx = 0; block_idx < _nodeTraversal.size(); block_idx++ ) {
    if ( _nodeTraversal[ block_idx ]->isDiagonal() ) {
      _diagonalBlocks.push_back( _nodeTraversal[ block_idx ] );
    } else {
      _offDiagonalBlocks.push_back( _nodeTraversal[ block_idx ] );
    }
  }

  if ( _diagonalBlocks.size() > 1 && explicitBlockThreshold > 0 ) {
    // Rebuild the node traversal so that it doesn't include children of
    // explicit blocks
    _nodeTraversal.clear();
    buildTraversal( _root, explicitBlockThreshold );
  } else {
    _explicitBlocks = false;
  }

  buildBlockList();

  // FIXME
  if ( _diagonalBlocks.size() > 1 && explicitBlockThreshold > 0 ) {
    for ( int i = 0; i < _blockList.size(); i++ ) {
      if ( _blockList[ i ]._buildType == EXPLICIT_BUILD ) {
        printf( "Block %d has %d columns; type EXPLICIT\n", i,
            _blockList[ i ]._isDiagonal ?
              _diagonalBlocks[_blockList[i]._index]->_fullBlock.numRows() :
              _offDiagonalBlocks[_blockList[i]._index]->_fullBlock.numRows() );
      } else {
        printf( "Block %d has %d columns; type RECURSIVE\n", i,
            _blockList[ i ]._isDiagonal ?
              _diagonalBlocks[_blockList[i]._index]->_fullBlock.numRows() :
              _offDiagonalBlocks[_blockList[i]._index]->_fullBlock.numRows() );
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////
// Initialize diagonal block data entries
//////////////////////////////////////////////////////////////////////
void DiagonalBlock::initData( Supernode &node,
                              const vector<long int> &diagonalOffsets )
{
  TRACE_ASSERT( _diagonalBlocks.size() == diagonalOffsets.size() );

  _node = &node;

  for ( int block_idx = 0; block_idx < _diagonalBlocks.size();
        block_idx++ )
  {
    // Some error checking
    if ( block_idx < _diagonalBlocks.size() - 1 ) {
      const DenseBlock      &block = _diagonalBlocks[ block_idx ]->_block;

      TRACE_ASSERT(
        block.numRows() * block.numRows()
          == diagonalOffsets[ block_idx + 1 ] - diagonalOffsets[ block_idx ] );
    }

    _diagonalBlocks[ block_idx ]->_data
      = _node->data() + diagonalOffsets[ block_idx ];
  }
}

//////////////////////////////////////////////////////////////////////
// Update the given diagonal block using previously computed
// off-diagonal blocks
//////////////////////////////////////////////////////////////////////
void DiagonalBlock::updateDiagonal( int block_idx )
{
  BlockNode                 *diagonalNode = _diagonalBlocks[ block_idx ];

#if 0
  BlockNode                 *ancestorNode = diagonalNode->_parent;
  BlockNode                 *currentNode = diagonalNode;
#endif

  const DenseBlock          &block = diagonalNode->_block;

  updateDiagonal( diagonalNode, block, diagonalNode->_data );
#if 0
  int                        nRows = block.numRows();

#if 0
  // FIXME: debugging
  if ( diagonalNode->_parent != NULL || nRows == 1144 ) {
    char buf[ 1024 ];
    cout << "Writing!!!!!" << endl;
    sprintf( buf, "super_numeric/block_%d_preupdate.matrix", block_idx );

    MATRIX::write( diagonalNode->_data, nRows, nRows, buf );
  }
#endif

  TRACE_ASSERT( diagonalNode->_data != NULL );

  // Walk up the tree to the root.
  for ( ; ancestorNode != NULL;
        currentNode = ancestorNode, ancestorNode = ancestorNode->_parent )
  {
    // Only ancestors to the right of this node in the tree
    // influence its contents
    if ( ancestorNode->_right != currentNode ) {
      TRACE_ASSERT( ancestorNode->_left == currentNode, "Invalid tree" );
      continue;
    }

    const DenseBlock        &ancestorBlock = ancestorNode->_block;
    int                      baseRow;

    const Real              *V;
    const Real              *U;

    // Workspace for forming the matrix update
    MATRIX                   workspace( nRows, ancestorBlock.numColumns() );

    TRACE_ASSERT( block._rowRange.first >= ancestorBlock._rowRange.first
               && block._rowRange.second <= ancestorBlock._rowRange.second,
               "Invalid block row range" );

    baseRow = block._rowRange.first - ancestorBlock._rowRange.first;

    // Skip to the desired row of the ancestor's low-rank decomposition
    TRACE_ASSERT( ancestorNode->_V && ancestorNode->_U,
                  "Undefined low-rank decomposition" );
    TRACE_ASSERT( ancestorNode->_rank > 0 );

    V = ancestorNode->_V + baseRow * ancestorNode->_rank;
    U = ancestorNode->_U;

    // Intermediate matrix
    MATRIX::gemm( V, U, workspace.data(),
                  // V and U dimensions
                  nRows, ancestorNode->_rank,
                  ancestorBlock.numColumns(), ancestorNode->_rank,
                  // Transpose U, but not V
                  false, true );

    // Symmetric update to the diagonal block
    MATRIX::syrk( workspace.data(), diagonalNode->_data,
                  nRows, ancestorBlock.numColumns(),
                  // No transpose
                  false,
                  // Subtract
                  -1.0, 1.0 );
  }
#endif
}

//////////////////////////////////////////////////////////////////////
// Updates an explicit schur complement formed for the given
// block.
//
// NOTE: block_entry_idx is an index in to _blockList
//////////////////////////////////////////////////////////////////////
void DiagonalBlock::updateDiagonalSchurComplement( int block_entry_idx,
                                                   Real *schurComplement )
{
  BlockNode                 *node;

  if ( _blockList[ block_entry_idx ]._isDiagonal ) {
    node = _diagonalBlocks[ _blockList[ block_entry_idx ]._index ];
  } else {
    node = _offDiagonalBlocks[ _blockList[ block_entry_idx ]._index ];
  }

  const DenseBlock          &block = node->_fullBlock;

  updateDiagonal( node, block, schurComplement );
}

#if 0
// FIXME: debugging
static Real DIAGONAL_COPY[ 1000000 ];
#endif

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
SolverErrorPtr DiagonalBlock::factorDiagonal( int block_idx )
{
  BlockNode                 *diagonalNode = _diagonalBlocks[ block_idx ];

  int                        nRows = diagonalNode->_block.numRows();
  int                        info;

#if 0
  // FIXME: debugging
  // Copy to diagonal copy
  MATRIX::copy( DIAGONAL_COPY, diagonalNode->_data, nRows, nRows );
#endif

#if 0
  // FIXME: debugging
  if ( diagonalNode->_parent != NULL || nRows == 1144 ) {
    char buf[ 1024 ];
    sprintf( buf, "super_numeric/block_%d_prefactor.matrix", block_idx );

    MATRIX::write( diagonalNode->_data, nRows, nRows, buf );
  }
#endif

  info = MATRIX::cholesky( diagonalNode->_data, nRows );

#if 0
  // FIXME: debugging
  if ( _node->compressOffDiagonal() ) {
    char buf[ 1024 ];
    sprintf( buf, "test/node_%d_block_%d.matrix", _node->nodeID(), block_idx );

    MATRIX::write( diagonalNode->_data, nRows, nRows, buf );
  }
#endif

#if 0
  // FIXME: debugging
  while ( info > 0 ) {
    printf( "Warning: indefinite matrix in row %d\n", info - 1 );

    // Bump this entry in the diagonal copy
    Real badEntry = MATRIX::access( diagonalNode->_data, nRows, nRows,
                                    info - 1, info - 1 );

    MATRIX::access( DIAGONAL_COPY, nRows, nRows, info - 1, info - 1 )
      += 11.0 * abs( badEntry );

    MATRIX::copy( diagonalNode->_data, DIAGONAL_COPY, nRows, nRows );

    info = MATRIX::cholesky( diagonalNode->_data, nRows );
  }
#endif

  if ( info < 0 ) {
    MATRIX::write( diagonalNode->_data, nRows, nRows, "invalidDiagonal.matrix" );
    TRACE_ASSERT( NULL, "ERROR: Invalid input to Cholesky solver" );
  } else if ( info > 0 ) {
    cerr << endl << "Block error in node " << _node->nodeID() << endl;

    SolverErrorPtr  error( new IndefiniteDiagonalError() );
    return error;

    MATRIX::write( diagonalNode->_data, nRows, nRows, "badDiagonal.matrix" );
    TRACE_ASSERT( NULL, "ERROR: Matrix not positive definite" );
    abort();
  }

  return SolverErrorPtr();
}


//////////////////////////////////////////////////////////////////////
// Given a workspace storing a Schur complement for the
// full block rooted at the given block index, factor this
// Schur complement and compress it's low-rank off-diagonal
// blocks
//////////////////////////////////////////////////////////////////////
SolverErrorPtr DiagonalBlock::compressExplicitBlock( int block_idx,
                                                     Real *schurComplement,
                                                     // Whether to decompose
                                                     // blocks directly or
                                                     // decompose their
                                                     // transposes
                                                     bool transpose,
                                                     // To determine ranks for
                                                     // off-diagonal blocks
                                                     RankEstimator rankEstimator,
                                                     // For storing compressed
                                                     // off-diagonal data
                                                     Real *extraData,
                                                     long int &offset,
                                                     long int &remainingSize,
                                                     int powerIterations )
{
#if 0
  printf( "      DiagonalBlock::compressExplicitBlock\n" );
#endif
  const BlockEntry          &entry = _blockList[ block_idx ];
  BlockNode                 *node;
  int                        info;

  if ( entry._isDiagonal ) {
    node = _diagonalBlocks[ entry._index ];
  } else {
    node = _offDiagonalBlocks[ entry._index ];
  }

  const DenseBlock          &block = node->_fullBlock;

  // Start by factoring the block
  info = MATRIX::cholesky( schurComplement, block.numRows() );
  if ( info < 0 ) {
    TRACE_ASSERT( NULL, "ERROR: Invalid input to Cholesky solver" );
  } else if ( info > 0 ) {
    printf( "Indefinite matrix found in node %d\n", _node->nodeID() );
    printf( "    Block range: [ %d, %d ]\n", block._rowRange.first,
                                             block._rowRange.second );

    SolverErrorPtr   error( new IndefiniteDiagonalError() );
    return error;

    MATRIX tmp( block.numRows(), block.numRows(), schurComplement );
    tmp.write( "badDiagonal.matrix" );

    TRACE_ASSERT( NULL, "Matrix is not positive definite" );
    abort();
  }

  // Pass this to the block to do compression and copy necessary data
  node->compressExplicitBlock( schurComplement, transpose, block,
                               rankEstimator, extraData,
                               offset, remainingSize, powerIterations );
#if 0
  printf( "      ........ done\n\n" );
#endif
  
  return SolverErrorPtr();
}

//////////////////////////////////////////////////////////////////////
// Multiply an off-diagonal block with a low-rank basis
//
// Note; this only accounts for contributions due to other off-diagonal
// blocks in this diagonal block.  Any other contributions are assumed
// to have been computed previously
//////////////////////////////////////////////////////////////////////
void DiagonalBlock::offDiagonalMultiply( int block_idx,
                                         const Real *G,
                                         int rowsG, int colsG,
                                         Real *output,
                                         bool left, bool transpose )
{
  BlockNode                 *offDiagonalNode = _offDiagonalBlocks[ block_idx ];

  BlockNode                 *ancestorNode = offDiagonalNode->_parent;
  BlockNode                 *currentNode = offDiagonalNode;

  const DenseBlock          &block = offDiagonalNode->_block;
  int                        nRows = block.numRows();
  int                        nCols = block.numColumns();

  MATRIX                     transposeWorkspace;
  MATRIX                     transposeOutputWorkspace;

  if ( left && transpose ) {
    TRACE_ASSERT( rowsG == block.numRows() );
  } else if ( left && !transpose ) {
    TRACE_ASSERT( rowsG == block.numColumns() );
  } else if ( !left && transpose ) {
    TRACE_ASSERT( colsG == block.numColumns() );
  } else {
    TRACE_ASSERT( colsG == block.numRows() );
  }

  // We at least 2 workspaces for multiplication of the form
  //
  // V * (U' * (U * (V' * G)))
  MATRIX                     workspaces[2];

  // Walk up the tree to the root.
  for ( ; ancestorNode != NULL;
        currentNode = ancestorNode, ancestorNode = ancestorNode->_parent )
  {
    // Only ancestors to the right of this node in the tree
    // influence its contents
    if ( ancestorNode->_right != currentNode ) {
      TRACE_ASSERT( ancestorNode->_left == currentNode, "Invalid tree" );
      continue;
    }

    const DenseBlock        &ancestorBlock = ancestorNode->_block;

    TRACE_ASSERT( block._rowRange.first >= ancestorBlock._rowRange.first
               && block._rowRange.second <= ancestorBlock._rowRange.second,
               "Invalid block row range" );
    TRACE_ASSERT( block._columnRange.first >= ancestorBlock._rowRange.first
               && block._columnRange.second <= ancestorBlock._rowRange.second,
               "Invalid block row range" );

    int                      baseRow;
    int                      baseColumn;

    const Real              *V1;
    const Real              *V2;
    const Real              *U;

    int                      V1rows = block.numRows();
    int                      V2rows = block.numColumns();
    int                      Urows = ancestorBlock.numColumns();

    int                      rank = ancestorNode->_rank;

    baseRow = block._rowRange.first - ancestorBlock._rowRange.first;
    baseColumn = block._columnRange.first - ancestorBlock._rowRange.first;

    V1 = ancestorNode->_V + baseRow * ancestorNode->_rank;
    V2 = ancestorNode->_V + baseColumn * ancestorNode->_rank;
    U = ancestorNode->_U;

    //////////////////////
    // Step 1
    //////////////////////
    if ( left ) {
      workspaces[0].resizeAndWipe( rank, colsG );
    } else {
      workspaces[0].resizeAndWipe( rowsG, rank );
    }

    if ( left && transpose ) {
      MATRIX::gemm( V1, G, workspaces[0].data(),
                    V1rows, rank, /* V1 dimensions */
                    rowsG, colsG, /* G dimensions */
                    true, false /* Transpose V1 */ );
    } else if ( left && !transpose ) {
      MATRIX::gemm( V2, G, workspaces[0].data(),
                    V2rows, rank, /* V2 dimensions */
                    rowsG, colsG, /* G dimensions */
                    true, false /* Transpose V2 */ );
    } else if ( !left && transpose ) {
      MATRIX::gemm( G, V2, workspaces[0].data(),
                    rowsG, colsG, /* G dimensions */
                    V2rows, rank, /* V2 dimensions */
                    false, false /* No transpose */ );
    } else {
      MATRIX::gemm( G, V1, workspaces[0].data(),
                    rowsG, colsG, /* G dimensions */
                    V1rows, rank, /* V2 dimensions */
                    false, false /* No transpose */ );
    }

    //////////////////////
    // Step 2
    //////////////////////
    if ( left ) {
      workspaces[1].resizeAndWipe( Urows, colsG );
    } else {
      workspaces[1].resizeAndWipe( rowsG, Urows );
    }

    if ( left ) {
      MATRIX::gemm( U, workspaces[0].data(), workspaces[1].data(),
                    Urows, rank, /* U dimensions */
                    workspaces[0].rows(), workspaces[0].cols(),
                    false, false /* No transpose */ );
    } else {
      MATRIX::gemm( workspaces[0].data(), U, workspaces[1].data(),
                    workspaces[0].rows(), workspaces[0].cols(),
                    Urows, rank, /* U dimensions */
                    false, true /* Transpose U */ );
    }

    //////////////////////
    // Step 3
    //////////////////////
    if ( left ) {
      workspaces[0].resizeAndWipe( rank, colsG );
    } else {
      workspaces[0].resizeAndWipe( rowsG, rank );
    }

    if ( left ) {
      MATRIX::gemm( U, workspaces[1].data(), workspaces[0].data(),
                    Urows, rank, /* U dimensions */
                    workspaces[1].rows(), workspaces[1].cols(),
                    true, false /* Transpose U */ );
    } else {
      MATRIX::gemm( workspaces[1].data(), U, workspaces[0].data(),
                    workspaces[1].rows(), workspaces[1].cols(),
                    Urows, rank, /* U dimensions */
                    false, false /* No transpose */ );
    }

    //////////////////////
    // Step 4
    //
    // Subtract final
    // updates from
    // output
    //////////////////////
    if ( left && transpose ) {
      MATRIX::gemm( V2, workspaces[0].data(), output,
                    V2rows, rank, /* V2 dimensions */
                    workspaces[0].rows(), workspaces[0].cols(),
                    false, false, /* No transpose */
                    -1.0, 1.0 /* Subtract */ );
    } else if ( left && !transpose ) {
      MATRIX::gemm( V1, workspaces[0].data(), output,
                    V1rows, rank, /* V1 dimensions */
                    workspaces[0].rows(), workspaces[0].cols(),
                    false, false, /* No transpose */
                    -1.0, 1.0 /* Subtract */ );
    } else if ( !left && transpose ) {
      MATRIX::gemm( workspaces[0].data(), V1, output,
                    workspaces[0].rows(), workspaces[0].cols(),
                    V1rows, rank, /* V1 dimensions */
                    false, true, /* Transpose V1 */
                    -1.0, 1.0 /* Subtract */ );
    } else {
      MATRIX::gemm( workspaces[0].data(), V2, output,
                    workspaces[0].rows(), workspaces[0].cols(),
                    V2rows, rank, /* V1 dimensions */
                    false, true, /* Transpose V2 */
                    -1.0, 1.0 /* Subtract */ );
    }
  }
}

//////////////////////////////////////////////////////////////////////
// Performs the lower-triangular solve associated with this
// off-diagonal block (determined by the tree structure)
//////////////////////////////////////////////////////////////////////
void DiagonalBlock::offDiagonalSolve( int block_idx,
                                      Real *G, int nRHS,
                                      bool transpose )
{
  BlockNode                 *offDiagonalNode = _offDiagonalBlocks[ block_idx ];

  // Do a solve with our left child, since this is what appears
  // "above" this block in the matrix
  diagonalSolve( offDiagonalNode->_left, G, nRHS, true /* left side */,
                 transpose );
}

//////////////////////////////////////////////////////////////////////
// Assigns low-rank decomposition data to an off-diagonal block
//////////////////////////////////////////////////////////////////////
void DiagonalBlock::assignOffDiagonalBlock( int block_idx,
                                            const Real *V, const Real *Utrans,
                                            int rank,
                                            Real *extraData,
                                            long int &offset,
                                            long int &remainingSize )
{
  BlockNode                 *offDiagonalNode = _offDiagonalBlocks[ block_idx ];

  const DenseBlock          &block = offDiagonalNode->_block;

  offDiagonalNode->_rank = rank;
  offDiagonalNode->_V = extraData + offset;
  offset += block.numRows() * rank;
  offDiagonalNode->_U = extraData + offset;
  offset += block.numColumns() * rank;

  remainingSize -= rank * ( block.numRows() + block.numColumns() );

  // Copy
  MATRIX::copy( offDiagonalNode->_V, V, block.numRows(), rank );

  MATRIX::transposeBLAS( offDiagonalNode->_U, Utrans,
                         rank, block.numColumns() );
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void DiagonalBlock::writeBlocks( const char *prefix )
{
  for ( int block_idx = 0; block_idx < _diagonalBlocks.size(); block_idx++ )
  {
    const DenseBlock &block = _diagonalBlocks[ block_idx ]->_block;

    if ( _diagonalBlocks[ block_idx ]->_data == NULL ) {
      continue;
    }

    char buf[ 1024 ];
    sprintf( buf, "%sdiag_block_%d.matrix", prefix, block_idx );

    MATRIX::write( _diagonalBlocks[ block_idx ]->_data,
                   block.numRows(), block.numColumns(), buf );
  }

  for ( int block_idx = 0; block_idx < _offDiagonalBlocks.size();
        block_idx++ )
  {
    const DenseBlock &block = _offDiagonalBlocks[ block_idx ]->_block;

    if ( _offDiagonalBlocks[ block_idx ]->_V == NULL
      || _offDiagonalBlocks[ block_idx ]->_U == NULL )
    {
      continue;
    }

    char buf[ 1024 ];
    sprintf( buf, "%soff_diag_block_V_%d.matrix", prefix, block_idx );

    MATRIX::write( _offDiagonalBlocks[ block_idx ]->_V,
                   block.numRows(), _offDiagonalBlocks[ block_idx ]->_rank,
                   buf );

    sprintf( buf, "%soff_diag_block_U_%d.matrix", prefix, block_idx );

    MATRIX::write( _offDiagonalBlocks[ block_idx ]->_U,
                   block.numColumns(), _offDiagonalBlocks[ block_idx ]->_rank,
                   buf );
  }
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void DiagonalBlock::buildTraversal( BlockNode *root,
                                    int explicitBlockThreshold )
{
  // Only recurse if this block's size is bigger than
  // the threshold
  if ( root->_left && root->_fullBlock.numRows() > explicitBlockThreshold ) {
    buildTraversal( root->_left, explicitBlockThreshold );
  }

  // If the size of block is small enough, identify it as an
  // EXPLICIT_BUILD block
  if ( root->_fullBlock.numRows() <= explicitBlockThreshold ) {
    _explicitBlocks = true;
    root->_buildType = EXPLICIT_BUILD;
  }

  _nodeTraversal.push_back( root );

  if ( root->_right && root->_fullBlock.numRows() > explicitBlockThreshold ) {
    buildTraversal( root->_right, explicitBlockThreshold );
  }
}

//////////////////////////////////////////////////////////////////////
// Assumes that buildTraversal has already been called
// and that _diagonalBlocks and _offDiagonalBlocks have been
// populated
//////////////////////////////////////////////////////////////////////
void DiagonalBlock::buildBlockList()
{
  map<BlockNode *, int>      indexMap;

  for ( int block_idx = 0; block_idx < _diagonalBlocks.size(); block_idx++ )
  {
    indexMap[ _diagonalBlocks[ block_idx ] ] = block_idx;
  }

  for ( int block_idx = 0; block_idx < _offDiagonalBlocks.size(); block_idx++ )
  {
    indexMap[ _offDiagonalBlocks[ block_idx ] ] = block_idx;
  }

  _blockList.resize( numBlocks() );
  for ( int block_idx = 0; block_idx < _nodeTraversal.size(); block_idx++ ) {
    _blockList[ block_idx ]._isDiagonal
          = _nodeTraversal[ block_idx ]->isDiagonal();

    _blockList[ block_idx ]._buildType
          = _nodeTraversal[ block_idx ]->buildType();

    TRACE_ASSERT(
      indexMap.find( _nodeTraversal[ block_idx ] ) != indexMap.end() );

    _blockList[ block_idx ]._index = indexMap[ _nodeTraversal[ block_idx ] ];
  }
}

//////////////////////////////////////////////////////////////////////
// Diagonal solve centered at this node
//////////////////////////////////////////////////////////////////////
void DiagonalBlock::diagonalSolve( BlockNode *root,
                                   Real *rhs, int nRHS,
                                   bool left, bool transpose ) const
{
  int                        ldRHS;

  if ( left ) {
    ldRHS = nRHS;
  } else {
    ldRHS = root->_fullBlock.numRows();
  }

#if 0
  // FIXME: debugging
  if ( root == _root && transpose && !left && _root->_left != NULL
    && _node->nodeID() == 9225 )
  {
    cout << "Doing debug write" << endl;

    // Dump some stuff out here...
    MATRIX::write( rhs, nRHS, ldRHS, "super_numeric/rhs_in.matrix" );

    for ( int block_idx = 0; block_idx < _diagonalBlocks.size(); block_idx++ )
    {
      const DenseBlock &block = _diagonalBlocks[ block_idx ]->_block;

      char buf[ 1024 ];
      sprintf( buf, "super_numeric/diag_block_%d.matrix", block_idx );

      MATRIX::write( _diagonalBlocks[ block_idx ]->_data,
                     block.numRows(), block.numColumns(), buf );
    }

    for ( int block_idx = 0; block_idx < _offDiagonalBlocks.size();
          block_idx++ )
    {
      const DenseBlock &block = _diagonalBlocks[ block_idx ]->_block;

      char buf[ 1024 ];
      sprintf( buf, "super_numeric/off_diag_block_V_%d.matrix", block_idx );

      MATRIX::write( _offDiagonalBlocks[ block_idx ]->_V,
                     block.numRows(), _offDiagonalBlocks[ block_idx ]->_rank,
                     buf );

      sprintf( buf, "super_numeric/off_diag_block_U_%d.matrix", block_idx );

      MATRIX::write( _offDiagonalBlocks[ block_idx ]->_U,
                     block.numColumns(), _offDiagonalBlocks[ block_idx ]->_rank,
                     buf );
    }
  }
#endif

  root->diagonalSolve( rhs, nRHS, ldRHS, left, transpose );

#if 0
  // FIXME: debugging
  if ( root == _root && transpose && !left && _root->_left != NULL
    && _node->nodeID() == 9225 )
  {
    // Dump some stuff out here...
    MATRIX::write( rhs, nRHS, ldRHS, "super_numeric/rhs_out.matrix" );
  }
#endif
}

// Diagonal solve with a vector right-hand side, centered at
// the given node.  This only allows solves from the left, and
// is meant to be used when applying the fully constructed
// factor matrix.
void DiagonalBlock::diagonalSolve( BlockNode *root, Real *rhs,
                                   bool transpose ) const
{
  root->diagonalSolve( rhs, transpose );
}

//////////////////////////////////////////////////////////////////////
// Deep copy of a diagonal block tree
//////////////////////////////////////////////////////////////////////
void DiagonalBlock::copy( const DiagonalBlock &block )
{
  map<const BlockNode *, BlockNode *>  parentMap;
  stack<const BlockNode *>             nodes;

  _node = block._node;

  parentMap[ NULL ] = NULL;

  _root = NULL;

  if ( block._root == NULL ) {
    return;
  }

  // Start by building the full node set, mapping existing nodes
  // to new nodes
  nodes.push( block._root );
  while ( !nodes.empty() ) {
    const BlockNode         *current = nodes.top();
    BlockNode               *newNode = new BlockNode();

    nodes.pop();

    newNode->_block = current->_block;
    newNode->_fullBlock = current->_fullBlock;
    newNode->_buildType = current->_buildType;

    // This data can be copied directly because block nodes aren't
    // responsible for it
    newNode->_data = current->_data;
    newNode->_V = current->_V;
    newNode->_U = current->_U;
    newNode->_rank = current->_rank;

    parentMap[ current ] = newNode;

    if ( current->_left != NULL ) {
      nodes.push( current->_left );
    }

    if ( current->_right != NULL ) {
      nodes.push( current->_right );
    }
  }

  // Next, build connections
  nodes.push( block._root );
  while ( !nodes.empty() ) {
    const BlockNode         *current = nodes.top();
    BlockNode               *newNode = parentMap[ current ];

    nodes.pop();

    newNode->_parent = parentMap[ current->_parent ];
    newNode->_left = parentMap[ current->_left ];
    newNode->_right = parentMap[ current->_right ];

    if ( current->_left != NULL ) {
      nodes.push( current->_left );
    }

    if ( current->_right != NULL ) {
      nodes.push( current->_right );
    }
  }

  _root = parentMap[ block._root ];

  // Finally, copy lists, node traversals etc.
  _diagonalBlocks.resize( block._diagonalBlocks.size() );
  _offDiagonalBlocks.resize( block._offDiagonalBlocks.size() );
  _nodeTraversal.resize( block._nodeTraversal.size() );

  for ( int i = 0; i < _diagonalBlocks.size(); i++ ) {
    _diagonalBlocks[ i ] = parentMap[ block._diagonalBlocks[ i ] ];
  }
  for ( int i = 0; i < _offDiagonalBlocks.size(); i++ ) {
    _offDiagonalBlocks[ i ] = parentMap[ block._offDiagonalBlocks[ i ] ];
  }
  for ( int i = 0; i < _nodeTraversal.size(); i++ ) {
    _nodeTraversal[ i ] = parentMap[ block._nodeTraversal[ i ] ];
  }

  _blockList = block._blockList;

  _explicitBlocks = block._explicitBlocks;
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void DiagonalBlock::clear()
{
  stack<BlockNode *>         nodes;

  if ( _root == NULL ) {
    return;
  }

  nodes.push( _root );
  while ( !nodes.empty() ) {
    BlockNode               *currentNode = nodes.top();

    nodes.pop();

    if ( currentNode->_left != NULL ) {
      nodes.push( currentNode->_left );
    }
    if ( currentNode->_right != NULL ) {
      nodes.push( currentNode->_right );
    }

    delete currentNode;
  }

  _root = NULL;
}

//////////////////////////////////////////////////////////////////////
// Called by update diagonal
//////////////////////////////////////////////////////////////////////
void DiagonalBlock::updateDiagonal( const BlockNode *node,
                                    const DenseBlock &block,
                                    Real *updateMatrix )
{
  const BlockNode           *ancestorNode = node->_parent;
  const BlockNode           *currentNode = node;

  int                        nRows = block.numRows();

#if 0
  // FIXME: debugging
  if ( diagonalNode->_parent != NULL || nRows == 1144 ) {
    char buf[ 1024 ];
    cout << "Writing!!!!!" << endl;
    sprintf( buf, "super_numeric/block_%d_preupdate.matrix", block_idx );

    MATRIX::write( diagonalNode->_data, nRows, nRows, buf );
  }
#endif

  // Walk up the tree to the root.
  for ( ; ancestorNode != NULL;
        currentNode = ancestorNode, ancestorNode = ancestorNode->_parent )
  {
    // Only ancestors to the right of this node in the tree
    // influence its contents
    if ( ancestorNode->_right != currentNode ) {
      TRACE_ASSERT( ancestorNode->_left == currentNode, "Invalid tree" );
      continue;
    }

    const DenseBlock        &ancestorBlock = ancestorNode->_block;
    int                      baseRow;

    const Real              *V;
    const Real              *U;

    // Workspace for forming the matrix update
    MATRIX                   workspace( nRows, ancestorBlock.numColumns() );

    TRACE_ASSERT( block._rowRange.first >= ancestorBlock._rowRange.first
               && block._rowRange.second <= ancestorBlock._rowRange.second,
               "Invalid block row range" );

    baseRow = block._rowRange.first - ancestorBlock._rowRange.first;

    // Skip to the desired row of the ancestor's low-rank decomposition
    TRACE_ASSERT( ancestorNode->_V && ancestorNode->_U,
                  "Undefined low-rank decomposition" );
    TRACE_ASSERT( ancestorNode->_rank > 0 );

    V = ancestorNode->_V + baseRow * ancestorNode->_rank;
    U = ancestorNode->_U;

    // Intermediate matrix
    MATRIX::gemm( V, U, workspace.data(),
                  // V and U dimensions
                  nRows, ancestorNode->_rank,
                  ancestorBlock.numColumns(), ancestorNode->_rank,
                  // Transpose U, but not V
                  false, true );

    // Symmetric update to the diagonal block
    MATRIX::syrk( workspace.data(), updateMatrix,
                  nRows, ancestorBlock.numColumns(),
                  // No transpose
                  false,
                  // Subtract
                  -1.0, 1.0 );
  }
}

//////////////////////////////////////////////////////////////////////
// BlockNode implementation
//////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////
// Constructor
//
// Empty block
//////////////////////////////////////////////////////////////////////
DiagonalBlock::BlockNode::BlockNode()
  : _parent( NULL ),
    _left( NULL ),
    _right( NULL ),
    _data( NULL ),
    _V( NULL ),
    _U( NULL ),
    _rank( 0 ),
    _buildType( RECURSIVE_BUILD )
{
}

//////////////////////////////////////////////////////////////////////
// Constructor
//
// Treat this as a leaf node, initially
//////////////////////////////////////////////////////////////////////
DiagonalBlock::BlockNode::BlockNode( const IndexRange &rowRange,
                                     const IndexRange &columnRange )
  : _block( rowRange, columnRange ),
    _fullBlock( rowRange, columnRange ),
    _parent( NULL ),
    _left( NULL ),
    _right( NULL ),
    _data( NULL ),
    _V( NULL ),
    _U( NULL ),
    _rank( 0 ),
    _buildType( RECURSIVE_BUILD )
{
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
DiagonalBlock::BlockNode::~BlockNode()
{
#if 0
  delete _left;
  delete _right;
#endif
}

//////////////////////////////////////////////////////////////////////
// Set children for this block
//////////////////////////////////////////////////////////////////////
void DiagonalBlock::BlockNode::setChildren( BlockNode *left, BlockNode *right )
{
  left->_parent = this;
  right->_parent = this;

  _left = left;
  _right = right;

  // Set the dense block for this node based on the children.
  //
  // Since this is not a leaf, it represents an off-diagonal block
  _block._rowRange.first = right->_fullBlock._rowRange.first;
  _block._rowRange.second = right->_fullBlock._rowRange.second;

  _block._columnRange.first = left->_fullBlock._columnRange.first;
  _block._columnRange.second = left->_fullBlock._columnRange.second;

  // Generate new row/column ranges for the full block represented by
  // this node
  _fullBlock._rowRange.first = left->_fullBlock._rowRange.first;
  _fullBlock._rowRange.second = right->_fullBlock._rowRange.second;

  _fullBlock._columnRange.first = left->_fullBlock._columnRange.first;
  _fullBlock._columnRange.second = right->_fullBlock._columnRange.second;
}

//////////////////////////////////////////////////////////////////////
// Diagonal solve centered at this node
//////////////////////////////////////////////////////////////////////
void DiagonalBlock::BlockNode::diagonalSolve( Real *rhs, int nRHS, int ldRHS,
                                              bool left, bool transpose ) const
{
  int                        nRows = _fullBlock.numRows();
#if 0
  int                        rowsRHS = left ? nRHS : nRows;
  int                        colsRHS = left ? nRows : nRHS;
#endif
  int                        rowsRHS = left ? nRows : nRHS;
  int                        colsRHS = left ? nRHS : nRows;

#if 0
  if ( nRHS == 1 ) {
    TRACE_ASSERT( left );
  }
#endif

  // If this is a leaf node, then we just need to do a straightforward
  // diagonal solve here
  if ( isDiagonal() ) {
    if ( nRHS == 1 && left ) {
      MATRIX::triangularSolve( _data, rhs, nRows,
                               true /* lower triangular */,
                               transpose,
                               false, /* Not unit diagonal */
                               1.0, ldRHS );
    } else {
      MATRIX::triangularSolve( _data, rhs, rowsRHS, colsRHS,
                               left, true /* lower triangular */,
                               transpose,
                               false, /* Not unit diagonal */
                               1.0, ldRHS );
    }

    return;
  }

  // The structure of the system is:
  //  /  L1       \
  //  \ V*U'  L2  /
  Real                      *Y1;
  Real                      *Y2;

  // Two parts of the RHS
  int                        dimY1 = _left->_fullBlock.numRows();
  int                        dimY2 = _right->_fullBlock.numRows();

  MATRIX                     workspace;

  if ( left ) {
    Y1 = rhs;

    // Align with correct row
    Y2 = rhs + dimY1 * colsRHS;
  } else {
    Y1 = rhs;

    // Align with correct column
    Y2 = rhs + dimY1;
  }

  if ( left && transpose ) {
    // Solve L2' X2 = Y2
    _right->diagonalSolve( Y2, nRHS, ldRHS, left, transpose );

    // Update Y1;   Y1 <-- Y1 - U * V' * X2
    workspace.resizeAndWipe( _rank, nRHS );
    MATRIX::gemm( _V, Y2, workspace.data(),
                  _block.numRows(), _rank, /* V dimensions */
                  dimY2, nRHS,
                  true, false, /* Transpose V */
                  1.0, 0.0, /* Overwrite */
                  -1, ldRHS, -1 /* Leading dimensions */ );

    MATRIX::gemm( _U, workspace.data(), Y1,
                  _block.numColumns(), _rank, /* U dimensions */
                  workspace.rows(), workspace.cols(),
                  false, false, /* No transpose */
                  -1.0, 1.0, /* Subtract */
                  -1, -1, ldRHS /* Leading dimensions */ );

    // Solve L1' X1 = Y1
    _left->diagonalSolve( Y1, nRHS, ldRHS, left, transpose );
  } else if ( left && !transpose ) {
    // Solve L1 X1 = Y1
    _left->diagonalSolve( Y1, nRHS, ldRHS, left, transpose );

    // Update Y2;   Y2 <-- Y2 - V * U' * X1
    workspace.resizeAndWipe( _rank, nRHS );
    MATRIX::gemm( _U, Y1, workspace.data(),
                  _block.numColumns(), _rank, /* U dimensions */
                  dimY1, nRHS,
                  true, false, /* Transpose U */
                  1.0, 0.0, /* Overwrite */
                  -1, ldRHS, -1 /* Leading dimensions */ );

    MATRIX::gemm( _V, workspace.data(), Y2,
                  _block.numRows(), _rank, /* V dimensions */
                  workspace.rows(), workspace.cols(),
                  false, false, /* No transpose */
                  -1.0, 1.0, /* Subtract */
                  -1, -1, ldRHS /* Leading dimensions */ );

    // Solve L2 X2 = Y2
    _right->diagonalSolve( Y2, nRHS, ldRHS, left, transpose );
  } else if ( !left && transpose ) {
    // Solve X1 L1' = Y1
    _left->diagonalSolve( Y1, nRHS, ldRHS, left, transpose );

    // Update Y2;   Y2 <-- Y2 - X1 * U * V'
    workspace.resizeAndWipe( nRHS, _rank );
    MATRIX::gemm( Y1, _U, workspace.data(),
                  nRHS, dimY1,
                  _block.numColumns(), _rank, /* U dimensions */
                  false, false, /* No transpose */
                  1.0, 0.0, /* Overwrite */
                  ldRHS, -1, -1 /* Leading dimensions */ );

    MATRIX::gemm( workspace.data(), _V, Y2,
                  workspace.rows(), workspace.cols(),
                  _block.numRows(), _rank, /* V dimensions */
                  false, true, /* Transpose V */
                  -1.0, 1.0, /* Subtract */
                  -1, -1, ldRHS /* Leading dimensions */ );

    // Solve X2 L2' = Y2
    _right->diagonalSolve( Y2, nRHS, ldRHS, left, transpose );
  } else {
    // Solve X2 L2 = Y2
    _right->diagonalSolve( Y2, nRHS, ldRHS, left, transpose );
    
    // Update Y1;   Y1 <-- Y1 - X2 * V * U'
    workspace.resizeAndWipe( nRHS, _rank );
    MATRIX::gemm( Y2, _V, workspace.data(),
                  nRHS, dimY2,
                  _block.numRows(), _rank, /* V dimensions */
                  false, false, /* No transpose */
                  1.0, 0.0,
                  ldRHS, -1, -1 /* Leading dimensions */ );

    MATRIX::gemm( workspace.data(), _U, Y1,
                  workspace.rows(), workspace.cols(),
                  _block.numColumns(), _rank, /* U dimensions */
                  false, true, /* Transpose U */
                  -1.0, 1.0, /* Subtract */
                  -1, -1, ldRHS /* Leading dimensions */ );

    // Solve X1 L1 = Y1
    _left->diagonalSolve( Y1, nRHS, ldRHS, left, transpose );
  }
}

// Diagonal solve (vector RHS) centered at this node.
// Only left-solves are allowed - this is to be used when
// solving systems with the factor matrix.
void DiagonalBlock::BlockNode::diagonalSolve( Real *rhs, bool transpose ) const
{
  int                        nRows = _fullBlock.numRows();
  int                        rowsRHS = nRows;
  int                        colsRHS = 1;

  // If this is a leaf node, then we just need to do a straightforward
  // diagonal solve here
  if ( isDiagonal() ) {
    MATRIX::triangularSolve( _data, rhs, nRows,
                             true /* lower triangular */,
                             transpose,
                             false, /* Not unit diagonal */
                             1.0, 1 );

    return;
  }

  // The structure of the system is:
  //  /  L1       \
  //  \ V*U'  L2  /
  Real                      *Y1;
  Real                      *Y2;

  // Two parts of the RHS
  int                        dimY1 = _left->_fullBlock.numRows();
  int                        dimY2 = _right->_fullBlock.numRows();

  MATRIX                     workspace;

  // Partition the right hand side
  Y1 = rhs;
  // Align with correct row
  Y2 = rhs + dimY1;

  if ( transpose ) {
    // Solve L2' X2 = Y2
    _right->diagonalSolve( Y2, transpose );

    // Update Y1;   Y1 <-- Y1 - U * V' * X2
    workspace.resizeAndWipe( _rank, 1 );
    MATRIX::gemv( _V, Y2, workspace.data(),
                  _block.numRows(), _rank, /* V dimensions */
                  true, /* Transpose V */
                  1.0, 0.0 /* Overwrite */ );

    MATRIX::gemv( _U, workspace.data(), Y1,
                  _block.numColumns(), _rank, /* U dimensions */
                  false, /* No transpose */
                  -1.0, 1.0 /* Subtract */ );

    // Solve L1' X1 = Y1
    _left->diagonalSolve( Y1, transpose );
  } else {
    // Solve L1 X1 = Y1
    _left->diagonalSolve( Y1, transpose );

    // Update Y2;   Y2 <-- Y2 - V * U' * X1
    workspace.resizeAndWipe( _rank, 1 );
    MATRIX::gemv( _U, Y1, workspace.data(),
                  _block.numColumns(), _rank, /* U dimensions */
                  true, /* Transpose U */
                  1.0, 0.0 /* Overwrite */ );

    MATRIX::gemv( _V, workspace.data(), Y2,
                  _block.numRows(), _rank, /* V dimensions */
                  false, /* No transpose */
                  -1.0, 1.0 /* Subtract */ );

    // Solve L2 X2 = Y2
    _right->diagonalSolve( Y2, transpose );
  }
}

//////////////////////////////////////////////////////////////////////
// Compress an explicitly provided matrix block in to this
// node's block structure
//////////////////////////////////////////////////////////////////////
void DiagonalBlock::BlockNode::compressExplicitBlock(
                                    Real *factor,
                                    // Whether to decompose blocks directly
                                    // or decompose their transposes
                                    bool transpose,
                                    // The matrix block stored in
                                    // factor
                                    const DenseBlock &factorBlock,
                                    // To determine ranks of off-diagonal
                                    // blocks
                                    RankEstimator rankEstimator,
                                    // For storing compressed off-diagonal data
                                    Real *extraData,
                                    long int &offset,
                                    long int &remainingSize,
                                    int powerIterations )
{
  int                        startRow;
  int                        startColumn;
  const Real                *factorData;
  
  startRow = _block._rowRange.first - factorBlock._rowRange.first;
  startColumn = _block._columnRange.first - factorBlock._columnRange.first;

  // The part of the factor data we need to work with
  factorData = factor + startRow * factorBlock.numColumns() + startColumn;

  // If this is a diagonal block, we can just copy the factored
  // data directly from factor
  if ( isDiagonal() ) {
    TRACE_ASSERT( startRow == startColumn );

    MATRIX::copy( _data, factorData, _block.numRows(), _block.numColumns(),
                  // Leading dimension for this data block
                  _block.numColumns(),
                  // Leading dimensions for the provided factor block
                  factorBlock.numColumns() );
  } else {
    int                      rank;
    MATRIX                   workspace;
    Real                    *workspace1;
    Real                    *workspace2;
    int                      info;
    FloatArray               qrExtraData( min( _block.numRows(),
                                               _block.numColumns() ) );

    // Compression rank for this block 
    rank = rankEstimator( _block.numRows(), _block.numColumns() );
    _rank = rank;

    workspace.resizeAndWipe( _block.numRows() + _block.numColumns(), rank );

    if ( transpose ) {
      workspace1 = workspace.data();
      workspace2 = workspace1 + _block.numRows() * rank;

      MathUtil::randomGaussianMatrix( _block.numRows(), rank, workspace1 );

#if 0
      // FIXME: debugging
      MATRIX decomposeBlock( _block.numRows(), _block.numColumns() );
      MATRIX::copy( decomposeBlock.data(), factorData,
                    _block.numRows(), _block.numColumns(),
                    // Leading dimensions
                    _block.numColumns(), factorBlock.numColumns() );
      decomposeBlock.write( "super_numeric/toBeDecomposed.matrix" );
#endif

#if 0
      // FIXME:
      MATRIX::write( workspace1, _block.numRows(), rank,
                     "super_numeric/multInput.matrix" );
#endif

      // Random multiplies with power iterations
      MATRIX::gemm( factorData, workspace1, workspace2,
                    // factorData and workspace dimensions
                    _block.numRows(), _block.numColumns(),
                    _block.numRows(), rank,
                    true, false, /* transpose factorData */
                    1.0, 0.0, /* overwrite */
                    factorBlock.numColumns() /* leading dimension */ );

#if 0
      // FIXME: debugging
      MATRIX::write( workspace2, _block.numColumns(), rank,
                     "super_numeric/multOutput1.matrix" );
#endif

      for ( int power_itr = 0; power_itr < powerIterations; power_itr++ ) {
        MATRIX::gemm( factorData, workspace2, workspace1,
                      // Dimensions
                      _block.numRows(), _block.numColumns(),
                      _block.numColumns(), rank,
                      false, false, /* No transpose */
                      1.0, 0.0, /* overwrite */
                      factorBlock.numColumns() /* leading dimension */ );

#if 0
        // FIXME: debugging
        MATRIX::write( workspace1, _block.numRows(), rank,
                       "super_numeric/multOutput2.matrix" );
#endif

        MATRIX::gemm( factorData, workspace1, workspace2,
                      // factorData and workspace dimensions
                      _block.numRows(), _block.numColumns(),
                      _block.numRows(), rank,
                      true, false, /* transpose factorData */
                      1.0, 0.0, /* overwrite */
                      factorBlock.numColumns() /* leading dimension */ );

#if 0
        // FIXME: debugging
        MATRIX::write( workspace2, _block.numColumns(), rank,
                       "super_numeric/multOutput3.matrix" );
#endif
      }

      // Turn workspace2 in to an orthonormal basis
      info = MATRIX::qr( workspace2, _block.numColumns(), rank,
                         qrExtraData.data(), NULL, -1 /* no workspace */ );
      TRACE_ASSERT( info == 0, "QR factorization failed" );

      info = MATRIX::extractQRfactor( workspace2, _block.numColumns(), rank,
                                      qrExtraData.data(), NULL, -1 );
      TRACE_ASSERT( info == 0, "QR extraction failed" );

#if 0
      // We are compressing the transpose, so get it, then compute
      // its rank-revealing QR factor
      MATRIX::transposeBLAS( workspace.data(), factorData,
                             _block.numRows(), _block.numColumns(),
                             -1 /* Ignore leading dimension */,
                             // Leading dimension for factor block input
                             factorBlock.numColumns() );

      MATRIX::qrPivot( workspace.data(),
                       NULL /* we don't care about pivots */,
                       _block.numColumns(), _block.numRows(),
                       qrExtraData.data(),
                       NULL, -1 /* no workspace needed */ );
      MATRIX::extractQRfactor( workspace.data(),
                               _block.numColumns(), _block.numRows(),
                               qrExtraData.data(),
                               NULL, -1 /* no workspace needed */ );
#endif

      // Copy the part of the basis that we're interested in
      // (the leading "rank" columns)
      _U = extraData + offset;
      offset += _block.numColumns() * rank;
#if 0
      MATRIX::copy( _U, workspace.data(), _block.numColumns(), rank,
                    rank /* Leading dimension of _U */,
                    _block.numRows() /* Leading dimension of workspace */ );
#endif
      MATRIX::copy( _U, workspace2, _block.numColumns(), rank );

      // Multiply to get the "V" part of the decomposition.  Only
      // consider the leading "rank" columns
      _V = extraData + offset;
      offset += _block.numRows() * rank;
      MATRIX::gemm( factorData, _U, _V,
                    // Dimensions of factor block
                    _block.numRows(), _block.numColumns(),
                    // Dimension of U
                    _block.numColumns(), rank,
                    // No transposition
                    false, false,
                    // Overwrite
                    1.0, 0.0,
                    // Leading dimension of factor block
                    factorBlock.numColumns() );
    } else {
      TRACE_ASSERT( NULL, "Not implemented yet" );
      // We are compressing the block itself, so get it, then compute
      // its rank-revealing QR factor
      MATRIX::copy( _data, factorData, _block.numRows(), _block.numColumns(),
                    -1 /* Ignore leading dimension */,
                    // Leading dimension for factor block input
                     factorBlock.numColumns() );

      MATRIX::qrPivot( workspace.data(),
                       NULL /* we don't care about pivots */,
                        _block.numRows(), _block.numColumns(),
                       qrExtraData.data(),
                       NULL, -1 /* no workspace needed */ );
      MATRIX::extractQRfactor( workspace.data(),
                               _block.numRows(), _block.numColumns(),
                               qrExtraData.data(),
                               NULL, -1 /* no workspace needed */ );

      // Copy the part of the basis that we're interested in
      // (the leading "rank" columns)
      _V = extraData + offset;
      offset += _block.numRows() * rank;
      MATRIX::copy( _V, workspace.data(), _block.numRows(), rank,
                    rank /* Leading dimension of V */,
                    _block.numColumns() /* Leading dimension of workspace */ );

      // Multiply to get the "U" part of the decomposition.
      _U = extraData + offset;
      offset += _block.numColumns() * rank;
      MATRIX::gemm( factorData, _V, _U,
                    // Dimensions of factor block
                    _block.numRows(), _block.numColumns(),
                    // Dimension of V
                    _block.numRows(), rank,
                    // Transpose factorBlock
                    true, false,
                    // Overwrite
                    1.0, 0.0,
                    // Leading dimension of factor block
                    factorBlock.numColumns() );
    }

    // Update the remaining size
    remainingSize -= rank * ( _block.numRows() + _block.numColumns() );

    // Recurse on our children
    //
    // TODO: try different tree orderings here, since it will affect
    // how things are laid out in memory.
    _left->compressExplicitBlock( factor, transpose, factorBlock,
                                  rankEstimator, extraData, offset,
                                  remainingSize, powerIterations );
    _right->compressExplicitBlock( factor, transpose, factorBlock,
                                   rankEstimator, extraData, offset,
                                   remainingSize, powerIterations );
  }
}
