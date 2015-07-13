//////////////////////////////////////////////////////////////////////
// Ordering.cpp: Implementation
//
//////////////////////////////////////////////////////////////////////

#include "Ordering.h"

#include <geometry/FiniteDifference.h>

#include <datastructure/KDTree.h>
#include <datastructure/SeparatorTree.h>

#if 0
#include <solver/SparseCholeskyFactor.h>
#endif

#include <util/IO.h>
#include <util/STLUtil.h>

#include <algorithm>
#include <set>

using namespace std;

//////////////////////////////////////////////////////////////////////
// Constructor
//////////////////////////////////////////////////////////////////////
Ordering::Ordering()
{
}

//////////////////////////////////////////////////////////////////////
// Destructor
//////////////////////////////////////////////////////////////////////
Ordering::~Ordering()
{
}

// Reorders a dense vector according to the given permutation
void Ordering::reorderVector( VECTOR &input, VECTOR &output,
                              const IntArray &permutation,
                              bool invert )
{
  if ( output.size() != input.size() )
  {
    output.resizeAndWipe( input.size() );
  }
  else
  {
    output.clear();
  }

  for ( int i = 0; i < input.size(); i++ )
  {
    if ( invert )
    {
      output( permutation[ i ] ) = input( i );
    }
    else
    {
      output( i ) = input( permutation[ i ] );
    }
  }
}

//////////////////////////////////////////////////////////////////////
// Reorders a dense square matrix according to the given permutation
//////////////////////////////////////////////////////////////////////
void Ordering::reorderMatrix( MATRIX &input, MATRIX &output,
                              const IntArray &permutation )
{
  TRACE_ASSERT( input.rows() == input.cols(),
                "Can't apply symmetric permutation to non-square matrix" );

  if ( output.rows() != input.rows() || output.cols() != input.cols() )
  {
    output.resizeAndWipe( input.rows(), input.cols() );
  }
  else
  {
    output.clear();
  }

  for ( int i = 0; i < input.rows(); i++ )
  {
    for ( int j = 0; j < input.cols(); j++ )
    {
      output( i, j ) = input( permutation[ i ], permutation[ j ] );
    }
  }
}

//////////////////////////////////////////////////////////////////////
// Reorders a sparse square matrix according to the given permutation
//////////////////////////////////////////////////////////////////////
void Ordering::reorderMatrix( SPARSE_MATRIX &input, SPARSE_MATRIX &output,
                              const IntArray &permutation )
{
  TRACE_ASSERT( input.rows() == input.cols(),
                "Can't apply symmetric permutation to non-square matrix" );

  IntArray permutationInverse( permutation.size() );

  for ( int i = 0; i < permutation.size(); i++ )
  {
    permutationInverse[ permutation[ i ] ] = i;
  }

  if ( output.rows() != input.rows() || output.cols() != input.cols() )
  {
    output.resize( input.rows(), input.cols() );
  }

  output.clearFull();

  for ( SPARSE_MATRIX::const_iterator i = input.begin(); i != input.end(); i++ )
  {
    int row_idx = permutationInverse[ i->first.first ];
    int col_idx = permutationInverse[ i->first.second ];

    int value = i->second;

    output( row_idx, col_idx ) = value;
  }
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void Ordering::reorderMatrix( SPARSE_MATRIX::SparseColumnMatrix &input,
                              SPARSE_MATRIX::SparseColumnMatrix &output,
                              const IntArray &permutation,
                              bool lowerOnly )
{
  TRACE_ASSERT( input._nrow == input._ncol,
                "Can't apply symmetric permutation to non-square matrix" );

  IntArray                   permutationInverse( permutation.size() );
  vector<pair<int, Real> >   columnEntries;
  int                        entryIndex;

  int                        old_col_idx;

  columnEntries.reserve( input._nrow );

  for ( int i = 0; i < permutation.size(); i++ )
  {
    permutationInverse[ permutation[ i ] ] = i;
  }

  // Allocates the output - we will overwrite all of its data
  if ( lowerOnly ) {
    output.clear();
    output._nrow = input._nrow;
    output._ncol = input._ncol;
    output._nzmax = input.symmetricNonZeros();
    output.allocate();
    output._p[ 0 ] = 0;
  } else {
    output = input;
    output._p[ 0 ] = 0;
  }

  entryIndex = 0;

  // Reorder data from the original matrix
  for ( int col_idx = 0; col_idx < input._ncol; col_idx++ )
  {
    old_col_idx = permutation[ col_idx ];

    // Get the entries from this column, with their new indices in the
    // matrix
    columnEntries.clear();
    for ( int row_ptr = input._p[ old_col_idx ];
          row_ptr < input._p[ old_col_idx + 1 ];
          row_ptr++ )
    {
      if ( !lowerOnly || permutationInverse[ input._i[ row_ptr ] ] >= col_idx )
      {
        columnEntries.push_back(
          pair<int, Real>( permutationInverse[ input._i[ row_ptr ] ],
                           input._x[ row_ptr ] ) );
      }
    }

    // Get the new column entries sorted in increasing order
    sort( columnEntries.begin(), columnEntries.end() );

    output._p[ col_idx + 1 ] = output._p[ col_idx ] + columnEntries.size();

    for ( int i = 0; i < columnEntries.size(); i++ )
    {
      output._i[ entryIndex ] = columnEntries[ i ].first;
      output._x[ entryIndex ] = columnEntries[ i ].second;

      entryIndex += 1;
    }
  }
}

//////////////////////////////////////////////////////////////////////
// Used to extract a set of submatrices from a sparse matrix
// using the provided range pairs.
// This will be pretty slow if we are trying to extract
// a lot of submatrices.  O(E*B*log(E)) to be specific, where
// E is the number of matrix entries and B is the number of
// submatrices to extract.
//////////////////////////////////////////////////////////////////////
void Ordering::extractSubMatrices( SPARSE_MATRIX &input,
                                   vector<SPARSE_MATRIX> &subMatrices,
                                   const vector<RangePair> &blocks,
                                   bool fullExtract )
{
  int                    n_rows, n_cols;
  int                    row_idx, col_idx;
  int                    sub_row_idx, sub_col_idx;
  Real                   value;

  // Set up matrix sizes
  subMatrices.resize( blocks.size() );

  for ( int block_num = 0; block_num < blocks.size(); block_num++ )
  {
    n_rows = blocks[ block_num ].first.second
            - blocks[ block_num ].first.first + 1;
    n_cols = blocks[ block_num ].second.second
            - blocks[ block_num ].second.first + 1;

    subMatrices[ block_num ].resize( n_rows, n_cols );
    subMatrices[ block_num ].clear();
  }

  int extractionErrors = 0;

  for ( SPARSE_MATRIX::const_iterator i = input.begin(); i != input.end(); i++ )
  {
    row_idx = i->first.first;
    col_idx = i->first.second;
    value = i->second;

    bool found = false;

    for ( int block_num = 0; block_num < blocks.size(); block_num++ )
    {
      const IndexRange  &row_range = blocks[ block_num ].first;
      const IndexRange  &col_range = blocks[ block_num ].second;
      SPARSE_MATRIX     &subMatrix = subMatrices[ block_num ];

      if ( row_idx >= row_range.first && row_idx <= row_range.second
        && col_idx >= col_range.first && col_idx <= col_range.second )
      {
        // Put this in the corresponding sub-matrix
        sub_row_idx = row_idx - row_range.first;
        sub_col_idx = col_idx - col_range.first;

        subMatrix( sub_row_idx, sub_col_idx ) = value;
        found = true;
      }
    }

    if ( fullExtract )
    {
      if ( !found && row_idx < col_idx )
      {
        cout << "ERROR: Nothing found for entry (" << row_idx << ", ";
        cout << col_idx << ") = " << value << endl;
        extractionErrors++;
      }
    }
  }

  if ( extractionErrors > 0 )
  {
    printf( "ERROR: %d out of %d entries not extracted.  Aborting.\n",
            extractionErrors, input.size() );
    abort();
  }
}

#if 0
//////////////////////////////////////////////////////////////////////
// Constructs a supernodal ordering on a matrix discretized
// using a finite difference grid.  We will make use of CHOLMOD's
// supernodal nested dissection ordering for leaf nodes in the
// nested dissection tree, and will impose the condition that
// each separator in our nested dissection tree corresponds to
// to supernode.
//
// Sparse matrix A is assumed to be discretized on a finite
// difference grid.  We are also making the assumption that
// the discrete differential operators used to construct A only
// involve immediate neighbours (eg. as in the 7 point Laplacian
// stencil) so that 
//
// gridDivisions stores the number of grid divisions along each axis
//////////////////////////////////////////////////////////////////////
SparseCholeskyFactor *Ordering::supernodalFiniteDifference( SPARSE_MATRIX &A,
                                                            Tuple3i gridDivisions,
                                                            int maxBlockSize,
                                                            int maxDiagonalBlock )
{
  NestedDissection::DissectionTree                  *orderTree;
  IntArray                                           permutation;
  IntArray                                           permutation_workspace;
  vector<FiniteDifference::GridPoint *>              separatorPermutation;
  vector<const NestedDissection::DissectionNode *>   nodeOrder;
  vector<const NestedDissection::DissectionNode *>   leafNodes;
  IntArray                                           numSeparatorLevels;
  vector<IntArray>                                   separatorNodeSizes;
  vector<IntArray>                                   separatorColumnStarts;
  IntArray                                           separatorStorageEstimate;
  const NestedDissection::DissectionNode            *treeNode;
  IndexRange                                         separatorIndices;
  IndexRange                                         columnRange;
  SPARSE_MATRIX                                      Areordered;
  vector<SPARSE_MATRIX>                              leafMatrices;
  vector<RangePair>                                  leafRanges;
  SPARSE_MATRIX::SparseColumnMatrix                  columnMatrix;
  SparseCholeskyFactor                              *factor = NULL;

  // List of the indices of extra rows we have to append to
  // each leaf supernode
  vector<vector<set<int > > >                        extraLeafRows;

  // List of indices of extra rows below the diagonal for
  // each separator block
  vector<set<int> >                                  extraSeparatorRows;

  // Symbolic factorizations of the leaf nodes
  vector<CHOLMOD_Environment::FactorWrapper>         leafFactors;

  // This is really hacky.  We'll have to figure out what to do here.
  // This basically controls how much extra space we allocate for
  // off-diagonal outer products in our separators.
  Real                                               storageConstant = 3.0;

  // We need a count of the number of supernodes, the number of subnodes
  // (which includes all subnodes in supernode trees), the number of
  // levels per node and the number of subnodes per node.
  // We also need the number of rows affected by each leaf node, and
  // each separator leaf node.
  int                                                nSuper = 0;
  int                                                nSubSuper = 0;
  IntArray                                           subNodesPerNode;
  vector<IntArray>                                   leafRowCounts;
  vector<IntArray>                                   separatorRowCounts;
  int                                                totalRowCount = 0;

  // maxLeafCount is the maximum row count for any subnode in a
  // given supernode tree
  IntArray                                           maxLeafCount;

  typedef SeparatorTree<FiniteDifference::GridPoint> PointTree;

  orderTree = NestedDissection::nestedDissectionFD( gridDivisions,
                                                    maxBlockSize,
                                                    10000, /* no max levels */
                                                    //true /* split largest */ );
                                                    false /* split smalles */ );

  // Build the high level permutation
  orderTree->buildPermutation( permutation );

  // Build sub-permutations for each separator in our nested
  // dissection tree
  orderTree->postOrder( nodeOrder );
  orderTree->getLeafNodes( leafNodes );

  for ( int i = 0; i < nodeOrder.size(); i++ )
  {
    treeNode = nodeOrder[ i ];

    // We are only interested in separators here, so skip the node if
    // it is a leaf.
    if ( treeNode->isLeaf() ) continue;

    separatorNodeSizes.push_back( IntArray() );

    separatorIndices.first = treeNode->right()->reorderedRange().second + 1;
    separatorIndices.second = treeNode->reorderedRange().second;

    separatorPermutation.clear();

    for ( int i = separatorIndices.first; i <= separatorIndices.second; i++ )
    {
      separatorPermutation.push_back(
        new FiniteDifference::GridPoint( permutation[ i ], gridDivisions ) );
    }

    // Build a KD tree to reorder the node.
    // We will just let the KD tree reorder the submatrix as much as
    // it wants, and pick what we want from the list
    // of leaf nodes later on.
    PointTree separatorTree( PointTree::LONGEST_AXIS,
                             1, /* largest leaf size */
                             PointTree::FIXED, -1.0 );

    separatorTree.build( separatorPermutation );

    separatorPermutation.clear();
    separatorTree.orderNodes( separatorPermutation );

    separatorTree.nodeSizes( separatorNodeSizes.back(), true );

    // TODO: Also need to get node sizes here so that we can use
    // them later on to remove off-diagonal blocks
    for ( int i = separatorIndices.first; i <= separatorIndices.second; i++ )
    {
      permutation[ i ]
        = separatorPermutation[ i - separatorIndices.first ]->index();
    }

    // This will delete all of the GridPoints we "new'd" earlier
    clearVectorContents( separatorPermutation );
  }

  // Reorder the matrix
  // TODO - it would be nice if we didn't have to do
  // anything like this, but it seems necessary for now.
  Ordering::reorderMatrix( A, Areordered, permutation );

  // Build a list of leaf blocks to extract from the matrix
  for ( int i = 0; i < leafNodes.size(); i++ )
  {
    treeNode = leafNodes[ i ];

    leafRanges.push_back(
        pair<IndexRange,IndexRange>( treeNode->reorderedRange(),
                                     treeNode->reorderedRange() ) );
  }

  // Extract all leaf submatrices from the main system
  extractSubMatrices( Areordered, leafMatrices, leafRanges );

  // Go through each leaf node and compute its symbolic, supernodal
  // factorization
  for ( int i = 0; i < leafMatrices.size(); i++ )
  {
    CHOLMOD_Environment solver;

    solver.setMatrix( leafMatrices[ i ] );
    solver.computeCholeskyFactor( false /* symbolic factorization only */ );

    leafFactors.push_back( CHOLMOD_Environment::FactorWrapper() );

    if ( !solver.getSupernodalFactor( leafFactors.back() ) )
    {
      cerr << "ERROR: Failed to get symbolic factorization" << endl;
      abort();
    }

    // Add to the total number of supernodes
    nSuper += leafFactors.back()._nsuper;
    nSubSuper += leafFactors.back()._nsuper;
  }

  TRACE_ASSERT( leafFactors.size() == leafNodes.size(),
                "Didn't get the right number of factors" );

  // Modify the current permutation according to the permutation
  // stored in each leaf's symbolic factorization
  for ( int i = 0; i < leafFactors.size(); i++ )
  {
    treeNode = leafNodes[ i ];

    int         *perm = leafFactors[ i ]._perm;
    IndexRange   idxRange = treeNode->reorderedRange();
    int          leafSize = idxRange.second - idxRange.first + 1;

    TRACE_ASSERT( perm, "No permutation found for leaf node" );

    permutation_workspace.clear();
    permutation_workspace.resize( leafSize );

    for ( int j = 0; j < leafSize; j++ )
    {
      permutation_workspace[ j ] = permutation[ idxRange.first + perm[ j ] ];
    }

    // Copy back to the full permutation
    for ( int j = 0; j < leafSize; j++ )
    {
      permutation[ idxRange.first + j ] = permutation_workspace[ j ];
    }
  }

  // Reorder the matrix (again... ugh)
  Ordering::reorderMatrix( A, Areordered, permutation );

  // Get the columns so that we can do row counts
  // FIXME: Check to make sure that this sorts rows in
  // ascending order.
  Areordered.constructSparseColumnCopy( columnMatrix );

  // Go through each leaf node and get a count of extra non-zero
  // row contributions from the lower triangular part of the matrix
  // (below the leaf block itself) which will contribute to super
  // nodes in the leaf.
  for ( int i = 0; i < leafFactors.size(); i++ )
  {
    extraLeafRows.push_back( vector<set<int> >() );
    extraLeafRows[ i ].resize( leafFactors[ i ]._nsuper );

    treeNode = leafNodes[ i ];

    IndexRange   idxRange = treeNode->reorderedRange();

    leafRowCounts.push_back( IntArray( leafFactors[ i ]._nsuper ) );

    // Iterate over each supernode in the leaf and figure out
    // which extra rows have to be appended
    for ( int j = 0; j < leafFactors[ i ]._nsuper; j++ )
    {
      int      col_start = leafFactors[ i ]._super[ j ] + idxRange.first;
      int      col_end = leafFactors[ i ]._super[ j + 1 ] + idxRange.first;
      int      last_row = idxRange.first + leafFactors[ i ]._s[
                                        leafFactors[ i ]._pi[ j + 1 ] - 1 ];

      // Check these columns in the original matrix to see
      // if there are any entries below this block that we
      // need to worry about
      for ( int col_idx = col_start; col_idx < col_end; col_idx++ )
      {
        for ( int row_p = columnMatrix._p[ col_idx ];
              row_p < columnMatrix._p[ col_idx + 1 ]; row_p++ )
        {
          int   row_idx = columnMatrix._i[ row_p ];

          //if ( row_idx >= col_end )
          if ( row_idx > last_row )
          {
            // We have an entry somewhere in the full matrix below
            // the leaf block in this supernode's column range.
            // We need to add this row's contribution to this supernode.
            extraLeafRows[ i ][ j ].insert( row_idx );
          }
        }
      }

      // Get a count of the total rows for this leaf supernode
      int nrows = leafFactors[ i ]._pi[ j + 1 ] - leafFactors[ i ]._pi[ j ];
      nrows += extraLeafRows[ i ][ j ].size();

      leafRowCounts[ i ][ j ] = nrows;
      totalRowCount += nrows;
    }
  }

  // Get the column range of each sub node resulting from the
  // sparsification of each separator
  numSeparatorLevels.resize( separatorNodeSizes.size() );
  subNodesPerNode.resize( separatorNodeSizes.size() );

  for ( int i = 0; i < separatorNodeSizes.size(); i++ )
  {
    int              level = 0;
    int              levelStart = 0;
    int              levelEnd = 1;
    int              levelSize = 1;
    IntArray        &nodeSizes = separatorNodeSizes[ i ];
    IntArray        &columnStarts = separatorColumnStarts[ i ];

    numSeparatorLevels[ i ] = -1;
    subNodesPerNode[ i ] = 0;

    while ( true )
    {
      bool           includeLevel = false;

      // We will include this level if it has any nodes larger
      // than the maximum diagonal size.

      // Ignore partial levels (this shouldn't happen for
      // reasonable amounts of sparsification anyways).
      if ( levelEnd > nodeSizes.size() ) break;

      for ( int j = levelStart; j < levelEnd && j < nodeSizes.size(); j++ )
      {
        if ( nodeSizes[ j ] >= maxDiagonalBlock )
        {
          includeLevel = true;
          break;
        }
      }

      if ( includeLevel )
      {
        // The first sub-node from each level starts at the
        // beginning of the parent supernode
        columnStarts.push_back( 0 );

        for ( int j = levelStart + 1; j < levelEnd; j++ )
        {
          columnStarts.push_back( columnStarts.back() + nodeSizes[ j - 1 ] );
        }

        numSeparatorLevels[ i ]++;
        subNodesPerNode[ i ] += levelSize;

        // Sanity check
        TRACE_ASSERT(
          columnStarts.back() + nodeSizes[ levelEnd - 1 ] == nodeSizes[ 0 ],
          "Invalid block sizes!" );
      }
      else
      {
        TRACE_ASSERT( level > 0,
                      "Schur complement smaller than max diagonal size" );

        // We're done
        break;
      }

      levelSize *= 2;
      levelStart = levelEnd;
      levelEnd = levelStart + levelSize;
      level++;
    }

    // Increment the number of supernodes and subnodes
    nSuper += 1;
    nSubSuper += subNodesPerNode[ i ];
  }

  // Figure out which extra rows are to be included below the
  int sepIndex = 0;

  for ( int i = 0; i < nodeOrder.size(); i++ )
  {
    treeNode = nodeOrder[ i ];

    // We are only interested in separators here, so skip the node if
    // it is a leaf.
    if ( treeNode->isLeaf() ) continue;

    separatorNodeSizes.push_back( IntArray() );

    separatorIndices.first = treeNode->right()->reorderedRange().second + 1;
    separatorIndices.second = treeNode->reorderedRange().second;

    // The ordering of separators in this array should
    // correspond to the ordering of separators from
    // the nested dissection tree, excluding the leaf nodes.
    int              nLevels = numSeparatorLevels[ sepIndex ];
    int              levelStart = 0;
    int              levelEnd = 1;
    int              levelSize = 1;
    IntArray        &columnStarts = separatorColumnStarts[ sepIndex ];
    IntArray        &nodeSizes = separatorNodeSizes[ sepIndex ];

    // Create an array of node counts for the separator supernode
    separatorRowCounts.push_back( IntArray() );
    maxLeafCount.push_back( 0 );

    // Get the start location for the bottom level
    //for ( int j = 0; j < nLevels - 1; j++ )
    for ( int j = 0; j < nLevels; j++ )
    {
      // Interior supernodes are empty
      for ( int levelNum = 0; levelNum < levelSize; levelNum++ )
      {
        separatorRowCounts[ sepIndex ].push_back( 0 );
      }

      levelSize *= 2;
      levelStart = levelEnd;
      levelEnd = levelStart + levelSize;
    }

    // For each leaf node, figure out the set of additional
    // rows affected by it
    for ( int j = levelStart; j < levelEnd; j++ )
    {
      int            col_start = columnStarts[ j ] + separatorIndices.first;
      int            col_end = col_start + nodeSizes[ j ];

      extraSeparatorRows.push_back( set<int>() );

      for ( int col_idx = col_start; col_idx < col_end; col_idx++ )
      {
        for ( int row_p = columnMatrix._p[ col_idx ];
              row_p < columnMatrix._p[ col_idx + 1 ]; row_p++ )
        {
          int        row_idx = columnMatrix._i[ row_p ];

          //if ( row_idx >= col_end )
          if ( row_idx > separatorIndices.second )
          {
            // We have an entry somewhere in the full matrix
            // below the leaf block in this sub-node's column range.
            extraSeparatorRows.back().insert( row_idx );
          }
        }
      }

      // Number of actual matrix rows for this leaf subnode
      int nrows = nodeSizes[ j ] + extraSeparatorRows.back().size();
      separatorRowCounts[ sepIndex ].push_back( nrows );

      maxLeafCount[ sepIndex ] = max( maxLeafCount[ sepIndex ], nrows );

      totalRowCount += nrows;
    }

    sepIndex++;
  }

  // Estimates of how many extra rows we want to allocate for
  // each separator supernode.
  for ( int i = 0; i < numSeparatorLevels.size(); i++ )
  {
    // Get the size of the separator in order to get a rough
    // estimate of how much storage we need
    int              baseSize = separatorNodeSizes[ i ][ 0 ];
    int              numRows = 0;

    Real             splitLength = sqrt( (Real)baseSize );

    for ( int j = 0; j < numSeparatorLevels[ i ]; j++ )
    {
      numRows += (int)( splitLength * storageConstant );

      splitLength /= 2.0;
    }

    separatorStorageEstimate.push_back( numRows );
  }

  // Allocate the factor
  factor = SparseCholeskyFactor::generateFactor();

  // Allocate the number of supernodes and subnodes, as well as
  // column offsets for subnodes, etc.
  factor->_nSuper = nSuper;
  factor->_nSubSuper = nSubSuper;
  factor->_super = (int *)malloc( (nSuper + 1) * sizeof( int ) );
  factor->_numLevels = (int *)malloc( nSuper * sizeof( int ) );
  factor->_numSubNodes = (int *)malloc( nSuper * sizeof( int ) );
  factor->_subSuper = (int *)malloc( nSubSuper * sizeof( int ) );
  factor->_subNodeSize = (int *)malloc( nSubSuper * sizeof( int ) );
  factor->_piIndex = (int *)malloc( nSuper * sizeof( int ) );
  factor->_pi = (int *)malloc( (nSubSuper + 1) * sizeof( int ) );
  factor->_s = (int *)malloc( totalRowCount * sizeof( int ) );
  factor->_px = (int *)malloc( nSuper * sizeof( int ) );
  factor->_nrows = (int *)malloc( nSuper * sizeof( int ) );
  factor->_columnOffset = (int *)malloc( nSubSuper * sizeof( int ) );
  factor->_rowOffset = (int *)malloc( nSubSuper * sizeof( int ) );
  factor->_mi = (int *)malloc( nSubSuper * sizeof( int ) );

  // Count up the number of subnodes in each supernode, multiplied
  // with the number of estimated storage for the supernode to get
  // a pessimistic upper bound on the amount of storage we need
  // for _rowOffset
  TRACE_ASSERT( subNodesPerNode.size() == separatorStorageEstimate.size(),
                "Data size mismatch" );

  int total_rowmap_rows = 0;
  for ( int i = 0; i < subNodesPerNode.size(); i++ )
  {
    total_rowmap_rows += subNodesPerNode[ i ] * separatorStorageEstimate[ i ];
  }

  factor->_rowMap = (int *)malloc( total_rowmap_rows * sizeof( int  ) );

  // Indices used to populate the factor structure
  int super_idx = 0;
  int sub_super_idx = 0;
  int leaf_factor_idx = 0;
  int separator_factor_idx = 0;
  int separator_leaf_idx = 0;
  int s_idx = 0;
  int s_idx_old;
  int x_size = 0;
  int row_map_idx = 0;

  // Populate these arrays
  for ( int i = 0; i < nodeOrder.size(); i++ )
  {
    treeNode = nodeOrder[ i ];

    // We are only interested in separators here, so skip the node if
    // it is a leaf.
    if ( treeNode->isLeaf() )
    {
      columnRange = treeNode->reorderedRange();

      CHOLMOD_Environment::FactorWrapper &leafFactor
        = leafFactors[ leaf_factor_idx ];

      addLeafSuperNodeSet( factor, columnRange, leafFactor, leafRowCounts,
                           extraLeafRows,
                           // These will be updated by addLeafSuperNodeSet
                           super_idx, sub_super_idx, leaf_factor_idx, s_idx,
                           x_size );

    }
    else
    {
      columnRange.first = treeNode->right()->reorderedRange().second + 1;
      columnRange.second = treeNode->reorderedRange().second;

      IntArray &nodeSizes = separatorNodeSizes[ separator_factor_idx ];
      IntArray &startColumns = separatorColumnStarts[ separator_factor_idx ];

      addSeparatorSuperNode( factor, columnRange, nodeSizes, startColumns,
                             numSeparatorLevels, subNodesPerNode,
                             separatorStorageEstimate, maxLeafCount,
                             separatorRowCounts, extraSeparatorRows,
                             // These will be updated by addSeparatorSuperNode
                             super_idx, sub_super_idx, separator_factor_idx,
                             separator_leaf_idx, s_idx, x_size );
    }
  }

  factor->_pi[ nSubSuper ] = totalRowCount;

  // TODO - still need to figure out how to setup
  // numRows, numFullRows, columnOffset, rowOffset, mi, rowMap and x
  // in factor

  return factor;
}
#endif

#if 0
//////////////////////////////////////////////////////////////////////
// Similar to the above function, but it just spits out a
// CHOLMOD_Environment object with a symbolic factor in which
// parts of Schur complements have been explicitly removed.
// This doesn't do anything to accomodate the sparsification
// of low diagonal blocks, etc.  It just removes them entirely.
//////////////////////////////////////////////////////////////////////
Ordering::SparseFactorInfo Ordering::supernodalFiniteDifferenceSparseCHOLMOD(
                                                 SPARSE_MATRIX &A,
                                                 Tuple3i gridDivisions,
                                                 int maxBlockSize,
                                                 int maxDiagonalBlock )
{
  NestedDissection::DissectionTree                  *orderTree;
  IntArray                                           permutation;
  IntArray                                           permutation_workspace;
  vector<FiniteDifference::GridPoint *>              separatorPermutation;
  vector<const NestedDissection::DissectionNode *>   nodeOrder;
  vector<const NestedDissection::DissectionNode *>   leafNodes;
  IntArray                                           numSeparatorLevels;
  vector<IntArray>                                   separatorNodeSizes;
  vector<IntArray>                                   separatorColumnStarts;
  const NestedDissection::DissectionNode            *treeNode;
  IndexRange                                         separatorIndices;
  IndexRange                                         columnRange;
  SPARSE_MATRIX                                      Areordered;
  vector<SPARSE_MATRIX>                              leafMatrices;
  vector<RangePair>                                  leafRanges;
  SPARSE_MATRIX::SparseColumnMatrix                  columnMatrix;

  // List of the indices of extra rows we have to append to
  // each leaf supernode
  vector<vector<set<int > > >                        extraLeafRows;

  // List of indices of extra rows below the diagonal for
  // each separator block
  vector<set<int> >                                  extraSeparatorRows;

  // Symbolic factorizations of the leaf nodes
  vector<CHOLMOD_Environment::FactorWrapper>         leafFactors;

  // We need a count of the number of supernodes, the number of subnodes
  // (which includes all subnodes in supernode trees), the number of
  // levels per node and the number of subnodes per node.
  // We also need the number of rows affected by each leaf node, and
  // each separator leaf node.
  int                                                nSuper = 0;
  int                                                nSubSuper = 0;
  IntArray                                           subNodesPerNode;
  vector<IntArray>                                   leafRowCounts;
  vector<IntArray>                                   separatorRowCounts;
  int                                                totalRowCount = 0;

  // maxLeafCount is the maximum row count for any subnode in a
  // given supernode tree
  IntArray                                           maxLeafCount;

  typedef SeparatorTree<FiniteDifference::GridPoint> PointTree;

  cout << "Building dissection tree..." << endl;

  orderTree = NestedDissection::nestedDissectionFD( gridDivisions,
                                                    maxBlockSize,
                                                    10000, /* no max levels */
                                                    true /* split largest */ );
                                                    //false /* split smallest */ );

  cout << "Building base permutation..." << endl;

  // Build the high level permutation
  orderTree->buildPermutation( permutation );

  // Build sub-permutations for each separator in our nested
  // dissection tree
  orderTree->postOrder( nodeOrder );
  orderTree->getLeafNodes( leafNodes );

  cout << "Building separator permutations..." << endl;

  for ( int i = 0; i < nodeOrder.size(); i++ )
  {
    treeNode = nodeOrder[ i ];

    // We are only interested in separators here, so skip the node if
    // it is a leaf.
    if ( treeNode->isLeaf() ) continue;

    separatorNodeSizes.push_back( IntArray() );

    separatorIndices.first = treeNode->right()->reorderedRange().second + 1;
    separatorIndices.second = treeNode->reorderedRange().second;

    printf( "Separator column range (%d, %d), size = %d\n",
            separatorIndices.first, separatorIndices.second,
            separatorIndices.second - separatorIndices.first + 1 );

    separatorPermutation.clear();

    for ( int i = separatorIndices.first; i <= separatorIndices.second; i++ )
    {
      separatorPermutation.push_back(
        new FiniteDifference::GridPoint( permutation[ i ], gridDivisions ) );
    }

    // Build a KD tree to reorder the node.
    // We will just let the KD tree reorder the submatrix as much as
    // it wants, and pick what we want from the list
    // of leaf nodes later on.
    //
    // FIXME - the distance (-0.1) that we are setting here
    // could cause problems for general grid spacing.
    PointTree separatorTree( PointTree::LONGEST_AXIS,
                             1, /* largest leaf size */
                             PointTree::FIXED, 0.0 );

    separatorTree.build( separatorPermutation );

    separatorPermutation.clear();
    separatorTree.orderNodes( separatorPermutation );

    separatorTree.nodeSizes( separatorNodeSizes.back(), true );

    // TODO: Also need to get node sizes here so that we can use
    // them later on to remove off-diagonal blocks
    for ( int i = separatorIndices.first; i <= separatorIndices.second; i++ )
    {
      permutation[ i ]
        = separatorPermutation[ i - separatorIndices.first ]->index();
    }

    // This will delete all of the GridPoints we "new'd" earlier
    clearVectorContents( separatorPermutation );
  }

  cout << "Reordering matrix..." << endl;

  // Reorder the matrix
  // TODO - it would be nice if we didn't have to do
  // anything like this, but it seems necessary for now.
  Ordering::reorderMatrix( A, Areordered, permutation );

  cout << "Getting leaf block index ranges..." << endl;

  // Build a list of leaf blocks to extract from the matrix
  for ( int i = 0; i < leafNodes.size(); i++ )
  {
    treeNode = leafNodes[ i ];

    leafRanges.push_back(
        pair<IndexRange,IndexRange>( treeNode->reorderedRange(),
                                     treeNode->reorderedRange() ) );
  }

  cout << "Extracting sub-matrices..." << endl;

  // Extract all leaf submatrices from the main system
  extractSubMatrices( Areordered, leafMatrices, leafRanges );

  cout << "Computing leaf block symbolic factors..." << endl;

#if 0
  for ( int i = 0; i < permutation.size(); i++ )
  {
    printf( "Original permutation[ %d ] = %d\n", i, permutation[ i ] );
  }
#endif

  // Go through each leaf node and compute its symbolic, supernodal
  // factorization
  for ( int i = 0; i < leafMatrices.size(); i++ )
  {
    CHOLMOD_Environment solver;

    solver.setMatrix( leafMatrices[ i ] );
    solver.computeCholeskyFactor( false /* symbolic factorization only */ );

    leafFactors.push_back( CHOLMOD_Environment::FactorWrapper() );

    if ( !solver.getSupernodalFactor( leafFactors.back() ) )
    {
      cerr << "ERROR: Failed to get symbolic factorization" << endl;
      abort();
    }

    // FIXME
    if ( i == 0 )
    {
      cout << "Computing leaf 0 Cholesky factor!" << endl;
      solver.computeCholeskyFactor( true, false );
      cout << endl;
    }

#if 0
    for ( int j = 0; j < leafFactors.back()._nsuper; j++ )
    {
      printf( "Leaf %d supernode %d starts at column %d\n", i, j, leafFactors.back()._super[ j ] );
    }
#endif

    // Add to the total number of supernodes
    nSuper += leafFactors.back()._nsuper;
    nSubSuper += leafFactors.back()._nsuper;
  }

  TRACE_ASSERT( leafFactors.size() == leafNodes.size(),
                "Didn't get the right number of factors" );

  cout << "Updating leaf permutations..." << endl;

  // Modify the current permutation according to the permutation
  // stored in each leaf's symbolic factorization
  for ( int i = 0; i < leafFactors.size(); i++ )
  {
    treeNode = leafNodes[ i ];

    int         *perm = leafFactors[ i ]._perm;
    IndexRange   idxRange = treeNode->reorderedRange();
    int          leafSize = idxRange.second - idxRange.first + 1;

    TRACE_ASSERT( perm, "No permutation found for leaf node" );

    permutation_workspace.clear();
    permutation_workspace.resize( leafSize );

    for ( int j = 0; j < leafSize; j++ )
    {
#if 0
      printf( "Leaf %d permutation[ %d ] = %d\n", i, j, perm[ j ] );
#endif
      permutation_workspace[ j ] = permutation[ idxRange.first + perm[ j ] ];
    }

    // Copy back to the full permutation
    for ( int j = 0; j < leafSize; j++ )
    {
      permutation[ idxRange.first + j ] = permutation_workspace[ j ];
    }
  }

#if 0
  for ( int i = 0; i < permutation.size(); i++ )
  {
    printf( "Permutation[ %d ] = %d\n", i, permutation[ i ] );
  }
#endif

  cout << "Reordering matrix..." << endl;

  // Reorder the matrix (again... ugh)
  Ordering::reorderMatrix( A, Areordered, permutation );

  cout << "Reformatting matrix..." << endl;

  // Get the columns so that we can do row counts
  // FIXME: Check to make sure that this sorts rows in
  // ascending order.
  Areordered.constructSparseColumnCopy( columnMatrix );

  cout << "Counting extra leaf supernode rows..." << endl;

  // Go through each leaf node and get a count of extra non-zero
  // row contributions from the lower triangular part of the matrix
  // (below the leaf block itself) which will contribute to super
  // nodes in the leaf.
  for ( int i = 0; i < leafFactors.size(); i++ )
  {
    extraLeafRows.push_back( vector<set<int> >() );
    extraLeafRows[ i ].resize( leafFactors[ i ]._nsuper );

    treeNode = leafNodes[ i ];

    IndexRange   idxRange = treeNode->reorderedRange();

    leafRowCounts.push_back( IntArray( leafFactors[ i ]._nsuper ) );

    // Iterate over each supernode in the leaf and figure out
    // which extra rows have to be appended
    for ( int j = 0; j < leafFactors[ i ]._nsuper; j++ )
    {
      int      col_start = leafFactors[ i ]._super[ j ] + idxRange.first;
      int      col_end = leafFactors[ i ]._super[ j + 1 ] + idxRange.first;
      int      last_row = idxRange.first + leafFactors[ i ]._s[
                                        leafFactors[ i ]._pi[ j + 1 ] - 1 ];

      // Check these columns in the original matrix to see
      // if there are any entries below this block that we
      // need to worry about
      for ( int col_idx = col_start; col_idx < col_end; col_idx++ )
      {
        for ( int row_p = columnMatrix._p[ col_idx ];
              row_p < columnMatrix._p[ col_idx + 1 ]; row_p++ )
        {
          int   row_idx = columnMatrix._i[ row_p ];

          //if ( row_idx >= col_end )
          if ( row_idx > last_row )
          {
            // We have an entry somewhere in the full matrix below
            // the leaf block in this supernode's column range.
            // We need to add this row's contribution to this supernode.
            extraLeafRows[ i ][ j ].insert( row_idx );
          }
        }
      }

      // Get a count of the total rows for this leaf supernode
      int nrows = leafFactors[ i ]._pi[ j + 1 ] - leafFactors[ i ]._pi[ j ];
      nrows += extraLeafRows[ i ][ j ].size();

      leafRowCounts[ i ][ j ] = nrows;
      totalRowCount += nrows;
    }
  }

  cout << "Computing subnode column ranges..." << endl;

  // Get the column range of each sub node resulting from the
  // sparsification of each separator
  numSeparatorLevels.resize( separatorNodeSizes.size() );
  subNodesPerNode.resize( separatorNodeSizes.size() );
  separatorColumnStarts.resize( separatorNodeSizes.size() );

  for ( int i = 0; i < separatorNodeSizes.size(); i++ )
  {
    int              level = 0;
    int              levelStart = 0;
    int              levelEnd = 1;
    int              levelSize = 1;
    IntArray        &nodeSizes = separatorNodeSizes[ i ];
    IntArray        &columnStarts = separatorColumnStarts[ i ];

    numSeparatorLevels[ i ] = -1;
    subNodesPerNode[ i ] = 0;

    while ( true )
    {
      bool           includeLevel = false;

      // We will include this level if it has any nodes larger
      // than the maximum diagonal size.

      // Ignore partial levels (this shouldn't happen for
      // reasonable amounts of sparsification anyways).
      if ( levelEnd > nodeSizes.size() ) break;

      for ( int j = levelStart; j < levelEnd && j < nodeSizes.size(); j++ )
      {
        if ( nodeSizes[ j ] >= maxDiagonalBlock )
        {
          printf( "Node at level %d has size = %d.  Including level.\n",
                  level, nodeSizes[ j ] );
          includeLevel = true;
          break;
        }
      }

      if ( includeLevel )
      {
        // The first sub-node from each level starts at the
        // beginning of the parent supernode
        columnStarts.push_back( 0 );

        for ( int j = levelStart + 1; j < levelEnd; j++ )
        {
          columnStarts.push_back( columnStarts.back() + nodeSizes[ j - 1 ] );
        }

        numSeparatorLevels[ i ]++;
        subNodesPerNode[ i ] += levelSize;

        cout << SDUMP( numSeparatorLevels[ i ] ) << endl;
        cout << SDUMP( subNodesPerNode[ i ] ) << endl;

        // Sanity check
        TRACE_ASSERT(
          columnStarts.back() + nodeSizes[ levelEnd - 1 ] == nodeSizes[ 0 ],
          "Invalid block sizes!" );
      }
      else
      {
        TRACE_ASSERT( level > 0,
                      "Schur complement smaller than max diagonal size" );

        // We're done
        break;
      }

      levelSize *= 2;
      levelStart = levelEnd;
      levelEnd = levelStart + levelSize;
      level++;
    }

    // Increment the number of supernodes and subnodes
    nSuper += 1;
    nSubSuper += subNodesPerNode[ i ];
  }

  cout << "Computing extra rows for diagonal blocks..." << endl;

  // Figure out which extra rows are to be included below the
  int sepIndex = 0;

  for ( int i = 0; i < nodeOrder.size(); i++ )
  {
    treeNode = nodeOrder[ i ];

    // We are only interested in separators here, so skip the node if
    // it is a leaf.
    if ( treeNode->isLeaf() ) continue;

    separatorNodeSizes.push_back( IntArray() );

    separatorIndices.first = treeNode->right()->reorderedRange().second + 1;
    separatorIndices.second = treeNode->reorderedRange().second;

    // The ordering of separators in this array should
    // correspond to the ordering of separators from
    // the nested dissection tree, excluding the leaf nodes.
    int              nLevels = numSeparatorLevels[ sepIndex ];
    int              levelStart = 0;
    int              levelEnd = 1;
    int              levelSize = 1;
    IntArray        &columnStarts = separatorColumnStarts[ sepIndex ];
    IntArray        &nodeSizes = separatorNodeSizes[ sepIndex ];

    // Create an array of node counts for the separator supernode
    separatorRowCounts.push_back( IntArray() );
    maxLeafCount.push_back( 0 );

    // Get the start location for the bottom level
    //for ( int j = 0; j < nLevels - 1; j++ )
    for ( int j = 0; j < nLevels; j++ )
    {
      // Interior supernodes are empty
      for ( int levelNum = 0; levelNum < levelSize; levelNum++ )
      {
        separatorRowCounts[ sepIndex ].push_back( 0 );
      }

      levelSize *= 2;
      levelStart = levelEnd;
      levelEnd = levelStart + levelSize;
    }

    // For each leaf node, figure out the set of additional
    // rows affected by it
    for ( int j = levelStart; j < levelEnd; j++ )
    {
      int            col_start = columnStarts[ j ] + separatorIndices.first;
      int            col_end = col_start + nodeSizes[ j ];

      extraSeparatorRows.push_back( set<int>() );

      for ( int col_idx = col_start; col_idx < col_end; col_idx++ )
      {
        for ( int row_p = columnMatrix._p[ col_idx ];
              row_p < columnMatrix._p[ col_idx + 1 ]; row_p++ )
        {
          int        row_idx = columnMatrix._i[ row_p ];

          //if ( row_idx >= col_end )
          if ( row_idx > separatorIndices.second )
          {
            // We have an entry somewhere in the full matrix
            // below the leaf block in this sub-node's column range.
            extraSeparatorRows.back().insert( row_idx );
          }
        }
      }

      // Number of actual matrix rows for this leaf subnode
      int nrows = nodeSizes[ j ] + extraSeparatorRows.back().size();
      separatorRowCounts[ sepIndex ].push_back( nrows );

      maxLeafCount[ sepIndex ] = max( maxLeafCount[ sepIndex ], nrows );

      totalRowCount += nrows;
    }

    sepIndex++;
  }

  cout << "Building sparse symbolic factor..." << endl;

  Ordering::SparseFactorInfo factorData;

  factorData = buildSparseCHOLMODFactor(
                                   A, permutation, nodeOrder, leafNodes,
                                   numSeparatorLevels, separatorNodeSizes,
                                   separatorColumnStarts, extraLeafRows,
                                   extraSeparatorRows, leafFactors );

  return factorData;
}
#endif

//////////////////////////////////////////////////////////////////////
// Computes fill-in for a block matrix and returns the number of
// non-zeros needed to represend the factor
//////////////////////////////////////////////////////////////////////
long int Ordering::BlockFillIn(
                  vector<set<int> > &blockInteractions,
                  const IntArray &blockSizes )
{
  long int                   totalNNZ = 0;

  for ( int block_idx = 0; block_idx < blockInteractions.size(); block_idx++ )
  {
    const set<int>          &interactions = blockInteractions.at( block_idx );
    int                      nCols = blockSizes.at( block_idx );

    //totalNNZ += nCols * ( nCols + 1 ) / 2;
    totalNNZ += nCols * nCols;

    for ( set<int>::iterator iter = interactions.begin();
          iter != interactions.end(); iter++ )
    {
      TRACE_ASSERT( *iter > block_idx );

      totalNNZ += nCols * blockSizes.at( *iter );

      set<int> &neighbourInteractions = blockInteractions.at( *iter );

      set<int>::iterator nextIter = iter;
      nextIter++;

      // Every other neighbour after the current one fill be filled
      // in the current neighbour
      for ( ; nextIter != interactions.end(); nextIter++ )
      {
        neighbourInteractions.insert( *nextIter );
      }
    }
  }

  return totalNNZ;
}

//////////////////////////////////////////////////////////////////////
// Builds a symbolic block schur complement given a set of interactions
// between previously computed blocks and the desired variable space
//////////////////////////////////////////////////////////////////////
void Ordering::BuildBlockSchurComplement(
                  const vector<set<int> > &blockInteractions,
                  vector<set<int> > &schurInteractions,
                  int start_block, int end_block )
{
  int                        numBlocks = end_block - start_block;

  schurInteractions.clear();
  schurInteractions.resize( numBlocks );

  printf( "Ordering::BuildBlockSchurComplement: start_block = %d "
          "end_block = %d\n", start_block, end_block );

  for ( int block_idx = 0; block_idx < blockInteractions.size(); block_idx++ )
  {
    const set<int>          &interactions = blockInteractions.at( block_idx );

    for ( set<int>::const_iterator iter = interactions.begin();
          iter != interactions.end(); iter++ )
    {
      if ( *iter < start_block || *iter >= end_block )
      {
        //cerr << "*** Warning *** out of range" << endl;
        continue;
      }

      int                    newBlockIdx = *iter - start_block;

      set<int>::const_iterator nextIter = iter;
      nextIter++;

      // Connect to all subsequent interactions
      for ( ; nextIter != interactions.end(); nextIter++ )
      {
        int                  interactionBlockIdx = *nextIter - start_block;

        if ( *nextIter < start_block || *nextIter >= end_block )
        {
          //cerr << "*** Warning *** out of range" << endl;
          continue;
        }

        //cout << "Actually inserting something" << endl;

        schurInteractions.at( newBlockIdx ).insert( interactionBlockIdx );
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////
// Builds a symbolic block schur complement given a set of interactions
// between previously computed blocks and the desired variable space
//////////////////////////////////////////////////////////////////////
void Ordering::BuildBlockSchurComplement(
                  const vector<IntArray> &blockInteractions,
                  vector<set<int> > &schurInteractions,
                  int start_block, int end_block,
                  const vector<vector<IndexRange> > *blockColumnRanges )
{
  int                        numBlocks = end_block - start_block;

  schurInteractions.clear();
  schurInteractions.resize( numBlocks );

  if ( blockColumnRanges )
  {
    TRACE_ASSERT( blockColumnRanges->size() == blockInteractions.size() );
  }

  printf( "Ordering::BuildBlockSchurComplement: start_block = %d "
          "end_block = %d\n", start_block, end_block );

  for ( int block_idx = 0; block_idx < blockInteractions.size(); block_idx++ )
  {
    const IntArray          &interactions = blockInteractions.at( block_idx );

    const vector<IndexRange> *columnRanges = NULL;

    if ( blockColumnRanges )
    {
      columnRanges = &blockColumnRanges->at( block_idx );

      TRACE_ASSERT( columnRanges->size() == interactions.size() );
    }

#if 0
    for ( set<int>::const_iterator iter = interactions.begin();
          iter != interactions.end(); iter++ )
#endif
#if 0
    for ( IntArray::const_iterator iter = interactions.begin();
          iter != interactions.end(); iter++ )
#endif
    for ( int interaction_idx = 0; interaction_idx < interactions.size();
          interaction_idx++ )
    {
      int                    ancestor_block = interactions[ interaction_idx ];
      int                    next_interaction_idx;

      if ( ancestor_block < start_block || ancestor_block >= end_block )
      {
        //cerr << "*** Warning *** out of range" << endl;
        continue;
      }

      int                    newBlockIdx = ancestor_block - start_block;

#if 0
      set<int>::const_iterator nextIter = iter;
#endif
#if 0
      IntArray::const_iterator nextIter = iter;
      nextIter++;
#endif
      next_interaction_idx = interaction_idx + 1;

      // Connect to all subsequent interactions
#if 0
      for ( ; nextIter != interactions.end(); nextIter++ )
#endif
      for ( ; next_interaction_idx < interactions.size();
            next_interaction_idx++ )
      {
        int         next_ancestor_block = interactions[ next_interaction_idx ];

        int         interactionBlockIdx = next_ancestor_block - start_block;

        if ( next_ancestor_block < start_block
          || next_ancestor_block >= end_block )
        {
          //cerr << "*** Warning *** out of range" << endl;
          continue;
        }

        //cout << "Actually inserting something" << endl;

        // If we have a column range list, and the two column ranges
        // do not intersect then skip adding this piece of fill in to the
        // Schur complement
        if ( columnRanges
          && !range_overlap( columnRanges->at( interaction_idx ),
                             columnRanges->at( next_interaction_idx ) ) )
        {
          continue;
        }

        schurInteractions.at( newBlockIdx ).insert( interactionBlockIdx );
      }
    }
  }
}

#if 0
//////////////////////////////////////////////////////////////////////
// Helper function for adding a leaf supernode set to
// our factor structure given some current indices, etc.
//////////////////////////////////////////////////////////////////////
void Ordering::addLeafSuperNodeSet(
                            SparseCholeskyFactor *factor,
                            IndexRange columnRange,
                            CHOLMOD_Environment::FactorWrapper &leafFactor,
                            vector<IntArray> &leafRowCounts,
                            vector<vector<set<int> > > &extraLeafRows,
                            // These keep track of the current factor
                            // assembly state
                            int &super_idx,
                            int &sub_super_idx,
                            int &leaf_factor_idx,
                            int &s_idx,
                            int &x_size )
{
  int s_idx_old;

  // Add all supernodes from this leaf
  for ( int j = 0; j < leafFactor._nsuper; j++ )
  {
    factor->_super[ super_idx ]
      =   columnRange.first /* Nested Dissection offset */
        + leafFactor._super[ j ]; /* Supernode offset */
    factor->_numLevels[ super_idx ] = 0;
    factor->_numSubNodes[ super_idx ] = 1;
    factor->_piIndex[ super_idx ] = sub_super_idx;
    factor->_nrows[ super_idx ] = leafRowCounts[ leaf_factor_idx ][ j ];

    factor->_subSuper[ sub_super_idx ] = factor->_super[ super_idx ];
    factor->_subNodeSize[ sub_super_idx ]
      =   leafFactor._super[ j + 1 ]
        - leafFactor._super[ j ];
    factor->_pi[ sub_super_idx ] = s_idx;

    // This is not needed for normal supernodes
    factor->_numFullRows[ sub_super_idx ] = 0;

    // Increment the index in to s according to how many
    // rows this node needs
    s_idx_old = s_idx;
    s_idx += leafRowCounts[ leaf_factor_idx ][ j ];

    // Increment the size of the data array
    x_size += leafRowCounts[ leaf_factor_idx ][ j ]
            * factor->_subNodeSize[ sub_super_idx ];

    // Fill in the non-zero row pattern for this subnode
    TRACE_ASSERT( s_idx - s_idx_old
      == extraLeafRows[ leaf_factor_idx ][ j ].size()
       + factor->_subNodeSize[ sub_super_idx ],
       "Error: Size mismatch in row pattern" );

    // Assign base diagonal block rows first
    for ( int row_idx = 0;
          row_idx < factor->_subNodeSize[ sub_super_idx ]; row_idx++ )
    {
      factor->_s[ s_idx_old + row_idx ]
        = factor->_subSuper[ sub_super_idx ] + row_idx;
    }
    s_idx_old += factor->_subNodeSize[ sub_super_idx ];

    // Next, assign extra rows
    for ( set<int>::iterator row_idx
            = extraLeafRows[ leaf_factor_idx ][ j ].begin();
          row_idx != extraLeafRows[ leaf_factor_idx ][ j ].end();
          row_idx++ )
    {
      factor->_s[ s_idx_old ] = *row_idx;
      s_idx_old++;
    }

    TRACE_ASSERT( s_idx_old == s_idx,
                  "Error: Size mismatch in row pattern" );

    super_idx++;
    sub_super_idx++;
  }

  leaf_factor_idx++;
}

//////////////////////////////////////////////////////////////////////
// Helper function for adding a separator supernode to our
// factor structure given some current indices, etc.
//////////////////////////////////////////////////////////////////////
void Ordering::addSeparatorSuperNode( 
                            SparseCholeskyFactor *factor,
                            IndexRange columnRange,
                            IntArray &nodeSizes,
                            IntArray &startColumns,
                            IntArray &numSeparatorLevels,
                            IntArray &subNodesPerNode,
                            IntArray &separatorStorageEstimate,
                            IntArray &maxLeafCount,
                            vector<IntArray> &separatorRowCounts,
                            vector<set<int> > &extraSeparatorRows,
                            // These keep track of the current factor
                            // assembly state
                            int &super_idx,
                            int &sub_super_idx,
                            int &separator_factor_idx,
                            int &separator_leaf_idx,
                            int &s_idx,
                            int &x_size )
{
  int s_idx_old;

  factor->_super[ super_idx ] = columnRange.first;
  factor->_numLevels[ super_idx ]
    = numSeparatorLevels[ separator_factor_idx ];
  factor->_numSubNodes[ super_idx ]
    = subNodesPerNode[ separator_factor_idx ];
  factor->_piIndex[ super_idx ] = sub_super_idx;
  factor->_nrows[ super_idx ]
    =   separatorStorageEstimate[ separator_factor_idx ]
      + maxLeafCount[ separator_factor_idx ];

  x_size += ( columnRange.second - columnRange.first + 1 )
          * factor->_nrows[ super_idx ];

  int max_full_rows = 0;

  // Set starting columns for each subnode
  for ( int subNode_idx = 0; subNode_idx < nodeSizes.size(); subNode_idx++ )
  {
    factor->_subSuper[ sub_super_idx ]
      = startColumns[ subNode_idx ] + columnRange.first;
    factor->_subNodeSize[ sub_super_idx ] = nodeSizes[ subNode_idx ];

    factor->_pi[ sub_super_idx ] = s_idx;
    s_idx_old = s_idx;

    factor->_columnOffset[ sub_super_idx ] = startColumns[ subNode_idx ];

    int rowCount = separatorRowCounts[ separator_factor_idx ][ subNode_idx ];

    if ( rowCount > 0 )
    {
      // Get the maximum off-diagonal contribution from any leaf node
      max_full_rows = max( max_full_rows,
                         (int)extraSeparatorRows[ separator_leaf_idx ].size() );

      // We must be at a leaf node, so add to factor->_s
      for ( int col_idx = 0;
            col_idx < factor->_subNodeSize[ sub_super_idx ];
            col_idx++ )
      {
        factor->_s[ s_idx ] = factor->_subSuper[ sub_super_idx ] + col_idx;
        s_idx++;
      }

      for ( set<int>::iterator row_idx
              = extraSeparatorRows[ separator_leaf_idx ].begin();
            row_idx != extraSeparatorRows[ separator_leaf_idx ].end();
            row_idx++ )
      {
        factor->_s[ s_idx ] = *row_idx;
        s_idx++;
      }

      TRACE_ASSERT( s_idx - s_idx_old == rowCount,
                    "Error: Size mismatch in subnode row pattern" );

      separator_leaf_idx++;
    }

    sub_super_idx++;
  }

  factor->_numFullRows[ super_idx ] = max_full_rows;

  super_idx++;
}

//////////////////////////////////////////////////////////////////////
// Generates a CHOLMOD factor, in which parts of Schur complements
// have been explicitly sparsified.
//////////////////////////////////////////////////////////////////////
Ordering::SparseFactorInfo Ordering::buildSparseCHOLMODFactor(
                SPARSE_MATRIX &A,
                IntArray &permutation,
                vector<const NestedDissection::DissectionNode *> &nodeOrder,
                vector<const NestedDissection::DissectionNode *> &leafNodes,
                IntArray &numSeparatorLevels,
                vector<IntArray> &separatorNodeSizes,
                vector<IntArray> &separatorColumnStarts,
                vector<vector<set<int> > > &extraLeafRows,
                vector<set<int> > &extraSeparatorRows,
                vector<CHOLMOD_Environment::FactorWrapper> &leafFactors )
{
  CHOLMOD_Environment                     *env = new CHOLMOD_Environment();
  CHOLMOD_Environment                     *env_full = new CHOLMOD_Environment();
  cholmod_factor                          *factor;
  cholmod_factor                          *factor_full;
  const NestedDissection::DissectionNode  *treeNode;
  int                                      leafIndex = 0;
  int                                      separatorIndex = 0;
  int                                      diagBlockIndex = 0;
  IndexRange                               nodeRange;
  int                                      start_idx;
  int                                     *Super;
  int                                     *Lpi;
  int                                     *Ls;
  int                                      ncols;

  // Used to figure out how much space needs to be allocated
  // in various parts of the supernodal factorization
  int                                      nSuper = 0;
  int                                      nSuper_full = 0;
  int                                      totalRows = 0;
  int                                      totalRows_full = 0;
  int                                      totalXSize = 0;
  int                                      totalXSize_full = 0;
  IntArray                                 startColumns;
  IntArray                                 startColumns_full;
  vector<set<int> >                        nodeRows;
  vector<set<int> >                        nodeRows_full;
  int                                      maxesize = 1;
  int                                      maxesize_full = 1;
  int                                      maxcsize = 1;
  int                                      maxcsize_full = 1;

  // This is used when figuring out what additional fill-in
  // rows to add to a supernode.  For each supernode, fill
  // in rows are only added if their index is larger than this.
  // The purpose of this is to make sure we don't fill in
  // precisely the parts of separators that we have sparsified
  // earlier.
  IntArray                                 minFillInRows;
  IntArray                                 minFillInRows_full;

  Ordering::SparseFactorInfo               factorData;

  factorData._env = env;
  factorData._envFull = env_full;

  // TODO: Still need to set maxcsize

  env->setMatrix( A );
  env->generateBlankFactor();

  factor = env->factor();

  env_full->setMatrix( A );
  env_full->generateBlankFactor();

  factor_full = env_full->factor();

  // Change the CHOLMOD factor to the desired type.
  // (Symbolic, supernodal, Cholesky)
  cholmod_change_factor( CHOLMOD_PATTERN, 1, 1, 0, 0,
                         factor, env->common() );
  cholmod_change_factor( CHOLMOD_PATTERN, 1, 1, 0, 0,
                         factor_full, env_full->common() );

#if 0
  cout << "NODE ORDER:" << endl;
#endif

  // Walk through each supernode and figure out which rows
  // it affects, etc. (supernodes in this case will be the
  // supernodes produced by CHOLMOD in the Leaf nodes, as well
  for ( int i = 0; i < nodeOrder.size(); i++ )
  {
    treeNode = nodeOrder[ i ];

    if ( treeNode->isLeaf() )
    {
      start_idx = treeNode->reorderedRange().first;

      printf( "Leaf node: (%d, %d)\n",
              treeNode->reorderedRange().first,
              treeNode->reorderedRange().second );

      // Grab the set of supernodes from this leaf (computed
      // earlier by CHOLMOD)
      CHOLMOD_Environment::FactorWrapper &wrapper = leafFactors[ leafIndex ];

      factorData._leafNNZ += wrapper._xsize;

      Super = wrapper._super;
      Lpi = wrapper._pi;
      Ls = wrapper._s;

      for ( int super_idx = 0; super_idx < wrapper._nsuper; super_idx++ )
      {
#if 0
        if ( i == 0 )
        {
          cout << "Introducing leaf supernode " << nSuper;
        }
#endif
        nSuper++;
        nSuper_full++;
        startColumns.push_back( start_idx + Super[ super_idx ] );
        startColumns_full.push_back( start_idx + Super[ super_idx ] );
#if 0
        if ( i == 0 )
        {
          cout << "   startColumn = " << startColumns.back();
          cout << "   lastCol = " << start_idx + Super[ super_idx + 1 ] << endl;
        }
#endif

        ncols = Super[ super_idx + 1 ] - Super[ super_idx ];

        minFillInRows.push_back( startColumns.back() + ncols );
        minFillInRows_full.push_back( startColumns_full.back() + ncols );

        nodeRows.push_back( set<int>() );
        nodeRows_full.push_back( set<int>() );

        set<int> &rowSet = nodeRows.back();
        set<int> &rowSet_full = nodeRows_full.back();

#if 0
        if ( i == 0 )
        {
          cout << "Base rows..." << endl;
        }
#endif

        for ( int row_idx = Lpi[ super_idx ];
              row_idx < Lpi[ super_idx + 1 ]; row_idx++ )
        {
          rowSet.insert( Ls[ row_idx ] + start_idx );
          rowSet_full.insert( Ls[ row_idx ] + start_idx );
#if 0
          if ( i == 0 )
          {
            cout << Ls[ row_idx ] + start_idx << endl;
          }
#endif
        }
#if 0
        if ( i == 0 )
        {
          cout << endl;
        }
#endif

        set<int> &extraRows = extraLeafRows[ leafIndex ][ super_idx ];

#if 0
        if ( i == 0 )
        {
          cout << "Extra rows..." << endl;
        }
#endif

        for ( set<int>::iterator extraRowIter = extraRows.begin();
              extraRowIter != extraRows.end(); extraRowIter++ )
        {
          rowSet.insert( *extraRowIter );
          rowSet_full.insert( *extraRowIter );
#if 0
          if ( i == 0 )
          {
            cout << *extraRowIter << endl;
          }
#endif
        }
#if 0
        if ( i == 0 )
        {
          cout << endl;
        }
#endif

#if 0
        maxesize = max( maxesize, (int)rowSet.size() - ncols );

        totalRows += rowSet.size();
        totalXSize += rowSet.size() * ncols;
#endif
      }

      leafIndex++;
    }
    else
    {
      int            nLevels = numSeparatorLevels[ separatorIndex ];
      int            levelStart = 0;
      int            levelEnd = 1;
      int            levelSize = 1;
      IntArray      &columnStarts = separatorColumnStarts[ separatorIndex ];
      IntArray      &nodeSizes = separatorNodeSizes[ separatorIndex ];

      IndexRange tempRange( treeNode->right()->reorderedRange().second + 1,
                            treeNode->reorderedRange().second );
      printf( "Separator node: (%d, %d)\n", tempRange.first, tempRange.second );

      factorData._separatorNNZ += ( tempRange.second - tempRange.first + 1 )
                                * ( tempRange.second - tempRange.first + 1 );

      factorData._separatorRanges.push_back( tempRange );

      factorData._finalSeparatorRange = tempRange;

      start_idx = treeNode->right()->reorderedRange().second + 1;

      // Add a single node in to the full symbolic factorization
      nSuper_full++;
      startColumns_full.push_back( start_idx );
      minFillInRows_full.push_back( tempRange.second + 1 );

      nodeRows_full.push_back( set<int>() );

      set<int> &rowSet_full = nodeRows_full.back();

      // Get to the start location of the bottom level
      //for ( int j = 0; j < nLevels - 1; j++ )
      for ( int j = 0; j < nLevels; j++ )
      {
        levelSize *= 2;
        levelStart = levelEnd;
        levelEnd = levelStart + levelSize;
      }

      // Create a supernode for each leaf node
      for ( int j = levelStart; j < levelEnd; j++ )
      {
        ncols = nodeSizes[ j ];

        factorData._sparseSeparatorNNZ += ncols * ncols;

        nodeRows.push_back( set<int>() );

        set<int> &rowSet = nodeRows.back();

        cout << "Introducing separator supernode " << nSuper;

        // Only fill in below the maximum separator index, NOT
        // the maximum index for this particular sub-block
        minFillInRows.push_back( treeNode->reorderedRange().second + 1 );

        nSuper++;
        startColumns.push_back( start_idx + columnStarts[ j ] );

        cout << "   startColumn = " << startColumns.back() << endl;

        for ( int row_idx = 0; row_idx < ncols; row_idx++ )
        {
          rowSet.insert( start_idx + columnStarts[ j ] + row_idx );
          rowSet_full.insert( start_idx + columnStarts[ j ] + row_idx );
        }

        set<int> &extraRows = extraSeparatorRows[ diagBlockIndex ];

        for ( set<int>::iterator extraRowIter = extraRows.begin();
              extraRowIter != extraRows.end(); extraRowIter++ )
        {
          rowSet.insert( *extraRowIter );
          rowSet_full.insert( *extraRowIter );
        }

#if 0
        maxesize = max( maxesize, (int)rowSet.size() - ncols );

        totalRows += rowSet.size();
        totalXSize += rowSet.size() * ncols;
#endif

        diagBlockIndex++;
      }

      separatorIndex++;
    }
  }

  startColumns.push_back( A.rows() );
  startColumns_full.push_back( A.rows() );

  // Figure out whatever additional fill in should exist
  cout << "Adding additional fill-in rows..." << endl;
  addFillInRows( startColumns, minFillInRows, nodeRows );
  addFillInRows( startColumns_full, minFillInRows_full, nodeRows_full );
  cout << "Done" << endl;

  // Compute a final count of rows, etc.
  for ( int i = 0; i < nodeRows.size(); i++ )
  {
    set<int> &rowSet = nodeRows[ i ];

    ncols = startColumns[ i + 1 ] - startColumns[ i ];

    totalRows += rowSet.size();
    totalXSize += rowSet.size() * ncols;

    maxesize = max( maxesize, (int)rowSet.size() - ncols );

    factorData._offDiagonalNNZ += ( rowSet.size() - ncols ) * ncols;
  }

  // Do the same thing for the full factor
  for ( int i = 0; i < nodeRows_full.size(); i++ )
  {
    set<int> &rowSet = nodeRows_full[ i ];

    ncols = startColumns_full[ i + 1 ] - startColumns_full[ i ];

    totalRows_full += rowSet.size();
    totalXSize_full += rowSet.size() * ncols;

    maxesize_full = max( maxesize_full, (int)rowSet.size() - ncols );
  }

  // Time to actually allocate the factor

  // Start by generating workspace for CHOLMOD.
  // TODO: I think that we don't need to worry about setting
  // any data in this workspace between calls, but I'm
  // not entirely sure of this.
  size_t w;
  int ok = 1;
  w = (size_t)A.rows() * 5;

  cholmod_allocate_work( (size_t)A.rows(), w, 0, env->common() );
  if ( env->common()->status < CHOLMOD_OK )
  {
    cerr << "Out of memory" << endl;
    return factorData;
  }

  // Do the same thing for the full factor
  cholmod_allocate_work( (size_t)A.rows(), w, 0, env_full->common() );
  if ( env_full->common()->status < CHOLMOD_OK )
  {
    cerr << "Out of memory" << endl;
    return factorData;
  }

  factor->n = (size_t)A.rows();
  factor->minor = (size_t)A.rows();

  CHOLMOD_Environment::clearCHOLMODFactor( factor );

  // Set up the factor permutation
  factor->ordering = CHOLMOD_GIVEN;
  factor->Perm = malloc( (size_t)A.rows() * sizeof( int ) );
  factor->ColCount = malloc( (size_t)A.cols() * sizeof( int ) );
  int *Perm = (int *)factor->Perm;
  int *ColCount = (int *)factor->ColCount;

  for ( int i = 0; i < permutation.size(); i++ )
  {
    Perm[ i ] = permutation[ i ];
  }

  // Set up the supernode list
  factor->nsuper = nSuper;
  factor->super = malloc( (nSuper + 1) * sizeof( int ) );
  int *super = (int *)factor->super;

  for ( int i = 0; i < startColumns.size(); i++ )
  {
    super[ i ] = startColumns[ i ];
  }

  // Set up pi, px and s arrays
  factor->ssize = totalRows;
  factor->xsize = totalXSize;

  factor->pi = malloc( (nSuper + 1) * sizeof( int ) );
  factor->px = malloc( (nSuper + 1) * sizeof( int ) );
  factor->s = malloc( factor->ssize * sizeof( int ) );
  
  Lpi = (int *)factor->pi;
  int *Lpx = (int *)factor->px;
  Ls = (int *)factor->s;

  // factor->x stays uninitialized until the actual
  // numerical factorization.

  int s_idx = 0;
  int x_idx = 0;

  for ( int i = 0; i < nodeRows.size(); i++ )
  {
    set<int> &rowSet = nodeRows[ i ];

    ncols = startColumns[ i + 1 ] - startColumns[ i ];

    Lpi[ i ] = s_idx;
    Lpx[ i ] = x_idx;

    // Initialize the s entry for this supernode.
    // TODO: Need to check that this orders rows in
    // ascending order.
    for ( set<int>::iterator rowSetIter = rowSet.begin();
          rowSetIter != rowSet.end(); rowSetIter++ )
    {
      Ls[ s_idx ] = *rowSetIter;
      s_idx++;
    }

    x_idx += rowSet.size() * ncols;
  }
  Lpi[ nodeRows.size() ] = s_idx;
  Lpx[ nodeRows.size() ] = x_idx;

  int *SuperMap = (int *)env->common()->Iwork;
  
  // Construct a supermap
  for ( int i = 0; i < nSuper; i++ )
  {
    for ( int j = startColumns[ i ]; j < startColumns[ i + 1 ]; j++ )
    {
      SuperMap[ j ] = i;
      ColCount[ j ] = Lpi[ i + 1 ] - Lpi[ i ];
    }

    // Maybe we should just set col counts for supernodes?
    //ColCount[ i ] = Lpi[ i + 1 ] - Lpi[ i ];
  }

  // Figure out maxcsize.
  // Copied (and modified) from cholmod_super_symbolic.c
  for (int i = 0; i < nSuper; i++)
  {
    ncols = startColumns[ i + 1 ] - startColumns[ i ];
    int p = Lpi[ i ] + ncols;
    int plast = p;
    int pend = Lpi[ i + 1 ];
    int slast = (p == pend) ? (-1) : (SuperMap [ Ls[p] ]);
    for ( ; p <= pend ; p++)
    {
      int s = (p == pend) ? (-1) : (SuperMap[ Ls[p] ]);
      if (s != slast)
      {
        /* row i is the start of a new supernode */
        int ndrow1 = p - plast;
        int ndrow2 = pend - plast;
        int csize = ndrow2 * ndrow1;
        maxcsize = max (maxcsize, csize);
        plast = p;
        slast = s;
      }
    }
  }

  factor->maxesize = maxesize;
  factor->maxcsize = maxcsize;

  if ( env->common()->status != CHOLMOD_OK )
  {
    cerr << "Error in CHOLMOD common!" << endl;
    abort();
  }

  ///////////////////////////////////////////////////////////////////////////
  // Repeat this whole procedure fo the full factor.
  // This is pretty bad in terms of code repetition...
  ///////////////////////////////////////////////////////////////////////////

  factor_full->n = (size_t)A.rows();
  factor_full->minor = (size_t)A.rows();

  CHOLMOD_Environment::clearCHOLMODFactor( factor_full );

  // Set up the factor permutation
  factor_full->ordering = CHOLMOD_GIVEN;
  factor_full->Perm = malloc( (size_t)A.rows() * sizeof( int ) );
  factor_full->ColCount = malloc( (size_t)A.cols() * sizeof( int ) );
  Perm = (int *)factor_full->Perm;
  ColCount = (int *)factor_full->ColCount;

  for ( int i = 0; i < permutation.size(); i++ )
  {
    Perm[ i ] = permutation[ i ];
  }

  // Set up the supernode list
  factor_full->nsuper = nSuper_full;
  factor_full->super = malloc( (nSuper_full + 1) * sizeof( int ) );
  super = (int *)factor_full->super;

  for ( int i = 0; i < startColumns_full.size(); i++ )
  {
    super[ i ] = startColumns_full[ i ];
  }

  // Set up pi, px and s arrays
  factor_full->ssize = totalRows_full;
  factor_full->xsize = totalXSize_full;

  factor_full->pi = malloc( (nSuper_full + 1) * sizeof( int ) );
  factor_full->px = malloc( (nSuper_full + 1) * sizeof( int ) );
  factor_full->s = malloc( factor_full->ssize * sizeof( int ) );
  
  Lpi = (int *)factor_full->pi;
  Lpx = (int *)factor_full->px;
  Ls = (int *)factor_full->s;

  // factor_full->x stays uninitialized until the actual
  // numerical factorization.

  s_idx = 0;
  x_idx = 0;

  for ( int i = 0; i < nodeRows_full.size(); i++ )
  {
    set<int> &rowSet = nodeRows_full[ i ];

    ncols = startColumns_full[ i + 1 ] - startColumns_full[ i ];

    Lpi[ i ] = s_idx;
    Lpx[ i ] = x_idx;

    // Initialize the s entry for this supernode.
    // TODO: Need to check that this orders rows in
    // ascending order.
    for ( set<int>::iterator rowSetIter = rowSet.begin();
          rowSetIter != rowSet.end(); rowSetIter++ )
    {
      Ls[ s_idx ] = *rowSetIter;
      s_idx++;
    }

    x_idx += rowSet.size() * ncols;
  }
  Lpi[ nodeRows_full.size() ] = s_idx;
  Lpx[ nodeRows_full.size() ] = x_idx;

  SuperMap = (int *)env_full->common()->Iwork;
  
  // Construct a supermap
  for ( int i = 0; i < nSuper_full; i++ )
  {
    for ( int j = startColumns_full[ i ]; j < startColumns_full[ i + 1 ]; j++ )
    {
      SuperMap[ j ] = i;
      ColCount[ j ] = Lpi[ i + 1 ] - Lpi[ i ];
    }

    // Maybe we should just set col counts for supernodes?
    //ColCount[ i ] = Lpi[ i + 1 ] - Lpi[ i ];
  }

  // Figure out maxcsize_full.
  // Copied (and modified) from cholmod_super_symbolic.c
  for (int i = 0; i < nSuper_full; i++)
  {
    ncols = startColumns_full[ i + 1 ] - startColumns_full[ i ];
    int p = Lpi[ i ] + ncols;
    int plast = p;
    int pend = Lpi[ i + 1 ];
    int slast = (p == pend) ? (-1) : (SuperMap [ Ls[p] ]);
    for ( ; p <= pend ; p++)
    {
      int s = (p == pend) ? (-1) : (SuperMap[ Ls[p] ]);
      if (s != slast)
      {
        /* row i is the start of a new supernode */
        int ndrow1 = p - plast;
        int ndrow2 = pend - plast;
        int csize = ndrow2 * ndrow1;
        maxcsize_full = max (maxcsize_full, csize);
        plast = p;
        slast = s;
      }
    }
  }

  factor_full->maxesize = maxesize_full;
  factor_full->maxcsize = maxcsize_full;

  if ( env_full->common()->status != CHOLMOD_OK )
  {
    cerr << "Error in CHOLMOD common!" << endl;
    abort();
  }

  return factorData;
}

//////////////////////////////////////////////////////////////////////
// Corrects the set of row indices associated with each super
// node to account for fill-in.
//
// This is quadratic in the number of supernodes, which is bad
//////////////////////////////////////////////////////////////////////
void Ordering::addFillInRows( vector<int> &startColumns,
                              vector<int> &minFillInRows,
                              vector<set<int> > &rowSets )
{
  int                      k2;
  int                      kd1, kd2;
  int                      d;
  vector<int>              supermap;
  int                      descendentSupernode;

  for ( int s = 0; s < startColumns.size() - 1; s++ )
  {
    for ( int col_idx = startColumns[s];
          col_idx < startColumns[s + 1]; col_idx++ )
    {
      supermap.push_back( s );
    }
  }

  for ( int s = 0; s < startColumns.size() - 1; s++ )
  {
    set<int>              &rowSet = rowSets[ s ];

    k2 = startColumns[ s + 1 ];

    descendentSupernode = -1;

    for ( set<int>::iterator rowSetIter = rowSet.begin();
          rowSetIter != rowSet.end(); rowSetIter++ )
    {
      if ( *rowSetIter < k2 ) continue;

      d = supermap[ *rowSetIter ];

      TRACE_ASSERT( d > s, "Ordering problem" );

      if ( d <= descendentSupernode )
      {
        // We've already done this one, so we can skip it.
        continue;
      }

      descendentSupernode = d;

      //kd2 = startColumns[ d + 1 ];
      kd2 = minFillInRows[ d ];

      set<int>            &dSet = rowSets[ d ];

      // Supernode s affects supernode d, so add the rest
      // of our row set to the row set for d.
      for ( set<int>::const_iterator subIter = rowSetIter;
            subIter != rowSet.end(); subIter++ )
      {
        if ( *subIter >= kd2 ) dSet.insert( *subIter );
      }
    }
  }
}
#endif
