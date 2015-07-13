//////////////////////////////////////////////////////////////////////
// Ordering.h: Interface for the Ordering class
//
//////////////////////////////////////////////////////////////////////

#ifndef ORDERING_H
#define ORDERING_H

#include <rschol/linearalgebra/MATRIX.h>
#include <rschol/linearalgebra/MATRIX3.h>
#include <rschol/linearalgebra/SPARSE_MATRIX.h>
#include <rschol/linearalgebra/VECTOR.h>

#include <rschol/SETTINGS.h>
#include <rschol/TYPES.h>

#include <rschol/library_interface/CHOLMOD_Environment.h>

#include <vector>

using namespace std;

class SparseCholeskyFactor;

//////////////////////////////////////////////////////////////////////
// Ordering class
//
// General matrix reordering routines
//////////////////////////////////////////////////////////////////////
class Ordering {
	public:
		Ordering();

		// Destructor
		virtual ~Ordering();

    //////////////////////////////////////////////////////////////////////
    // Simple specification structure for supernodes
    //////////////////////////////////////////////////////////////////////
    struct SupernodeSpecification {
      IndexRange             _columnRange;

      bool                   _compressOffDiagonal;
    };

    // Reorders a dense vector according to the given permutation
    static void reorderVector( VECTOR &input, VECTOR &output,
                               const IntArray &permutation,
                               bool invert = false );

    // Reorders a dense square matrix according to the given permutation
    static void reorderMatrix( MATRIX &input, MATRIX &output,
                               const IntArray &permutation );

    // Reorders a sparse square matrix according to the given permutation
    static void reorderMatrix( SPARSE_MATRIX &input, SPARSE_MATRIX &output,
                               const IntArray &permutation );

    static void reorderMatrix( SPARSE_MATRIX::SparseColumnMatrix &input,
                               SPARSE_MATRIX::SparseColumnMatrix &output,
                               const IntArray &permutation,
                               bool lowerOnly = false );

    // Used to extract a set of submatrices from a sparse matrix
    // using the provided range pairs.
    // This will be pretty slow if we are trying to extract
    // a lot of submatrices.  O(E*B*log(E)) to be specific, where
    // E is the number of matrix entries and B is the number of
    // submatrices to extract.
    //
    // if fullExtract is true, then this function will check
    // to make sure that everything in the matrix has actually
    // been extracted.
    static void extractSubMatrices( SPARSE_MATRIX &input,
                                    std::vector<SPARSE_MATRIX> &subMatrices,
                                    const std::vector<RangePair> &blocks,
                                    bool fullExtract = false );

#if 0
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
    // involve their immediate neighbours (eg. as in the 7 point Laplacian
    // stencil)
    //
    // divisions stores the number of grid divisions along each axis
    //
    // maxBlockSize is the smallest sparse block that we will handle
    // with conventional supernodal methods.  Anything larger we will
    // handle ourselves.
    //
    // maxDiagonalBlock represents the largest diagonal block we
    // will allow to appear in a Schur complement
    static SparseCholeskyFactor *supernodalFiniteDifference(
                                                     SPARSE_MATRIX &A,
                                                     Tuple3i divisions,
                                                     int maxBlockSize,
                                                     int maxDiagonalBlock );
#endif

    struct SparseFactorInfo {
      SparseFactorInfo()
        : _env( NULL ),
          _envFull( NULL ),
          _leafNNZ( 0 ),
          _offDiagonalNNZ( 0 ),
          _separatorNNZ( 0 ),
          _sparseSeparatorNNZ( 0 )
      {
      }

      CHOLMOD_Environment            *_env;
      CHOLMOD_Environment            *_envFull;

      long long int                   _leafNNZ;
      long long int                   _offDiagonalNNZ;
      long long int                   _separatorNNZ;
      long long int                   _sparseSeparatorNNZ;

      // Last separator column range
      IndexRange                      _finalSeparatorRange;

      vector<IndexRange>              _finalBlockRanges;

      vector<IndexRange>              _separatorRanges;

    };

#if 0
    // Similar to the above function, but it just spits out a
    // CHOLMOD_Environment object with a symbolic factor in which
    // parts of Schur complements have been explicitly removed.
    // This doesn't do anything to accomodate the sparsification
    // of low diagonal blocks, etc.  It just removes them entirely.
    static SparseFactorInfo supernodalFiniteDifferenceSparseCHOLMOD(
                                                     SPARSE_MATRIX &A,
                                                     Tuple3i divisions,
                                                     int maxBlockSize,
                                                     int maxDiagonalBlock );
#endif

    // Computes fill-in for a block matrix and returns the number of
    // non-zeros needed to represend the factor
    static long int BlockFillIn(
                      std::vector<std::set<int> > &blockInteractions,
                      const IntArray &blockSizes );

    // Builds a symbolic block schur complement given a set of interactions
    // between previously computed blocks and the desired variable space
    static void BuildBlockSchurComplement(
                      const std::vector<std::set<int> > &blockInteractions,
                      std::vector<std::set<int> > &schurInteractions,
                      int start_block, int end_block );

    // Builds a symbolic block schur complement given a set of interactions
    // between previously computed blocks and the desired variable space
    //
    // This version also accepts column ranges for the given block
    // interactions.  The purpose of this is to allow the user to specify
    // the exact range of column influence for an interaction and hence
    // possibly avoid introducing some interactions with numerically zero
    // value.
    static void BuildBlockSchurComplement(
        const std::vector<IntArray> &blockInteractions,
        std::vector<std::set<int> > &schurInteractions,
        int start_block, int end_block,
        const std::vector<std::vector<IndexRange> > *blockColumnRanges = NULL );

	protected:

	private:
#if 0
    // Helper function for adding a leaf supernode set to
    // our factor structure given some current indices, etc.
    static void addLeafSuperNodeSet(
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
                            int &x_size );

    // Helper function for adding a separator supernode to our
    // factor structure given some current indices, etc.
    static void addSeparatorSuperNode( 
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
                            int &x_size );

    // Generates a CHOLMOD factor, in which parts of Schur complements
    // have been explicitly sparsified.
    static SparseFactorInfo buildSparseCHOLMODFactor(
                    SPARSE_MATRIX &A,
                    IntArray &permutation,
                    vector<const NestedDissection::DissectionNode *> &nodeOrder,
                    vector<const NestedDissection::DissectionNode *> &leafNodes,
                    IntArray &numSeparatorLevels,
                    vector<IntArray> &separatorNodeSizes,
                    vector<IntArray> &separatorColumnStarts,
                    vector<vector<set<int> > > &extraLeafRows,
                    vector<set<int> > &extraSeparatorRows,
                    vector<CHOLMOD_Environment::FactorWrapper> &leafFactors );

    // Corrects the set of row indices associated with each super
    // node to account for fill-in.
    //
    // This is quadratic in the number of supernodes, which is bad
    static void addFillInRows( vector<int> &startColumns,
                               vector<int> &minFillInRows,
                               vector<set<int> > &rowSets );
#endif

};

#endif
