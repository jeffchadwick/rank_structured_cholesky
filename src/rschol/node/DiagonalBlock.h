//////////////////////////////////////////////////////////////////////
// DiagonalBlock.h: Interface for the DiagonalBlock class
//
//////////////////////////////////////////////////////////////////////

#ifndef DIAGONAL_BLOCK_H
#define DIAGONAL_BLOCK_H

#include <rschol/linearalgebra/DenseBlock.h>
#include <rschol/linearalgebra/MATRIX.h>
#include <rschol/linearalgebra/MATRIX3.h>
#include <rschol/linearalgebra/VECTOR.h>

#include <rschol/solver/SolverError.h>

#include <rschol/SETTINGS.h>
#include <rschol/TYPES.h>

#include <boost/function.hpp>

class Supernode;

//////////////////////////////////////////////////////////////////////
// DiagonalBlock class
//
// Models a hierarchically compressed diagonal block in a supernode
//////////////////////////////////////////////////////////////////////
class DiagonalBlock {
	public:
		DiagonalBlock();

    DiagonalBlock( const DiagonalBlock &block );

		// Destructor
		virtual ~DiagonalBlock();

    DiagonalBlock &operator=( const DiagonalBlock &block );

    // Build types for a diagonal block.  Either build a block
    // recursively, or explicitly construct all of its contents
    // and then do compression afterwards.
    enum BuildType {
      RECURSIVE_BUILD = 0,
      EXPLICIT_BUILD
    };

    // Initialize the block
    void init( const vector<DenseBlock> &diagonalBlocks,
               const vector<long int> &diagonalOffsets,
               int explicitBlockThreshold );

    // Initialize diagonal block data entries
    void initData( Supernode &node,
                   const vector<long int> &diagonalOffsets );

    // Update the given diagonal block using previously computed
    // off-diagonal blocks
    void updateDiagonal( int block_idx );

    // Updates an explicit schur complement formed for the given
    // block.
    //
    // NOTE: block_entry_idx is an index in to _blockList
    void updateDiagonalSchurComplement( int block_entry_idx,
                                        Real *schurComplement );

    SolverErrorPtr factorDiagonal( int block_idx );

    typedef boost::function<int (int, int)> RankEstimator;

    // Given a workspace storing a Schur complement for the
    // full block rooted at the given block index, factor this
    // Schur complement and compress it's low-rank off-diagonal
    // blocks
    SolverErrorPtr compressExplicitBlock( int block_idx, Real *schurComplement,
                                         // Whether to decompose blocks directly
                                         // or decompose their transposes
                                         bool transpose,
                                         // To determine ranks for off-diagonal
                                         // blocks
                                         RankEstimator rankEstimator,
                                         // For storing compressed off-diagonal
                                         // data
                                         Real *extraData,
                                         long int &offset,
                                         long int &remainingSize,
                                         int powerIterations );

    // Multiply an off-diagonal block with a low-rank basis
    //
    // Note; this only accounts for contributions due to other off-diagonal
    // blocks in this diagonal block.  Any other contributions are assumed
    // to have been computed previously
    void offDiagonalMultiply( int block_idx,
                              const Real *G,
                              int rowsG, int colsG,
                              Real *output,
                              bool left, bool transpose );

    // Performs the lower-triangular solve associated with this
    // off-diagonal block (determined by the tree structure)
    void offDiagonalSolve( int block_idx,
                           Real *G, int nRHS,
                           bool transpose );

    // Assigns low-rank decomposition data to an off-diagonal block
    void assignOffDiagonalBlock( int block_idx,
                                 const Real *V, const Real *Utrans,
                                 int rank,
                                 Real *extraData,
                                 long int &offset,
                                 long int &remainingSize );

    void writeBlocks( const char *prefix );

    void fullDiagonalSolve( Real *rhs, int nRHS,
                            bool transpose, bool left ) const
    {
      TRACE_ASSERT( _root );

      return diagonalSolve( _root, rhs, nRHS, left, transpose );
    }

    // Same as the above, but for left solves on vector data
    void fullDiagonalSolve( Real *rhs, bool transpose ) const
    {
      TRACE_ASSERT( _root );

      return diagonalSolve( _root, rhs, transpose );
    }

    int numDiagonalBlocks() const
    {
      return _diagonalBlocks.size();
    }

    int numOffDiagonalBlocks() const
    {
      return _offDiagonalBlocks.size();
    }

    int numBlocks() const
    {
      return _nodeTraversal.size();
    }

    const DenseBlock &diagonalBlock( int i ) const
    {
      return _diagonalBlocks[ i ]->_block;
    }

    const DenseBlock &offDiagonalBlock( int i ) const
    {
      return _offDiagonalBlocks[ i ]->_block;
    }

    const DenseBlock &fullBlock( int i ) const
    {
      const BlockEntry      &entry = _blockList[ i ];

      if ( entry._isDiagonal ) {
        return _diagonalBlocks[ entry._index ]->_fullBlock;
      } else {
        return _offDiagonalBlocks[ entry._index ]->_fullBlock;
      }
    }

    struct BlockEntry {
      bool                   _isDiagonal;
      int                    _index;
      BuildType              _buildType;
    };

    const BlockEntry blockEntry( int i ) const
    {
      return _blockList[ i ];
    }

    bool hasExplicitBlocks() const
    {
      return _explicitBlocks;
    }

	protected:

  private:
    // For representing the hierarchy of blocks in a diagonal block
    class BlockNode {
      public:
        BlockNode();
        BlockNode( const IndexRange &rowRange, const IndexRange &columnRange );

        virtual ~BlockNode();

        // Set children for this block
        void setChildren( BlockNode *left, BlockNode *right );

        bool isDiagonal() const
        {
          return ( _left == NULL && _right == NULL );
        }

        // Diagonal solve centered at this node
        void diagonalSolve( Real *rhs, int nRHS, int ldRHS,
                            bool left, bool transpose ) const;

        // Diagonal solve (vector RHS) centered at this node.
        // Only left-solves are allowed - this is to be used when
        // solving systems with the factor matrix.
        void diagonalSolve( Real *rhs, bool transpose ) const;

        // Compress an explicitly provided matrix block in to this
        // node's block structure
        void compressExplicitBlock( Real *factor,
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
                                    int powerIterations );

        inline BuildType buildType() const {
          return _buildType;
        }

      private:
        friend class DiagonalBlock;

        // The actual matrix block represented by this node
        DenseBlock           _block;

        // The matrix block represented by this node and all of
        // its children
        DenseBlock           _fullBlock;

        BlockNode           *_parent;
        BlockNode           *_left;
        BlockNode           *_right;

        // Block storage
        //
        // A single block in the case of a diagonal
        Real                *_data;

        // ... or two blocks for an off-diagonal low-rank
        // decomposition
        Real                *_V;
        Real                *_U;
        int                  _rank;

        // The build method for this block.  If the block is
        // large enough, we will build it recursively.  However,
        // if the block is small enough then even if it is a leaf
        // we can build it by explicitly constructing the full
        // factor for this block and compressing blocks of this
        // afterwards.
        BuildType            _buildType;

    };

  private:
    // If this block has size smaller than explicitBlockThreshold,
    // then we stop the traversal here and change it to have
    // _buildType == EXPLICIT_BUILD
    //
    // This represents the fact that this block will now be constructed
    // explicitly, prior to compression.
    void buildTraversal( BlockNode *root,
                         int explicitBlockThreshold = -1 );

    // Assumes that buildTraversal has already been called
    // and that _diagonalBlocks and _offDiagonalBlocks have been
    // populated
    void buildBlockList();

    // Diagonal solve centered at this node
    void diagonalSolve( BlockNode *root,
                        Real *rhs, int nRHS,
                        bool left, bool transpose ) const;

    // Diagonal solve with a vector right-hand side, centered at
    // the given node.  This only allows solves from the left, and
    // is meant to be used when applying the fully constructed
    // factor matrix.
    void diagonalSolve( BlockNode *root, Real *rhs, bool transpose ) const;

    // Deep copy of a diagonal block tree
    void copy( const DiagonalBlock &block );

    void clear();

    // Called by update diagonal
    void updateDiagonal( const BlockNode *node, const DenseBlock &block,
                         Real *updateMatrix );

	private:
    Supernode               *_node;

    BlockNode               *_root;

    // Lists of diagonal and off-diagonal blocks
    vector<BlockNode *>      _diagonalBlocks;
    vector<BlockNode *>      _offDiagonalBlocks;

    // In-order traversal
    vector<BlockNode *>      _nodeTraversal;

    // List of nodes mapping to indices in the block arrays
    vector<BlockEntry>       _blockList;

    bool                     _explicitBlocks;

};

#endif
