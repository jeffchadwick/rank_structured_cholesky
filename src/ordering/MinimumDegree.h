//////////////////////////////////////////////////////////////////////
// MinimumDegree.h: Interface for MinimumDegree class
//
//////////////////////////////////////////////////////////////////////

#ifndef MINIMUM_DEGREE_H
#define MINIMUM_DEGREE_H

#include <linearalgebra/MATRIX.h>
#include <linearalgebra/MATRIX3.h>
#include <linearalgebra/VECTOR.h>

#include <SETTINGS.h>
#include <TYPES.h>

#include <set>
#include <vector>

#include <map>
//#include <boost/unordered_map.hpp>

//////////////////////////////////////////////////////////////////////
// MinimumDegree class
//
// Static routines for performing minimum degree ordering, given
// some information about block interactions and sizes
//////////////////////////////////////////////////////////////////////
class MinimumDegree {
	public:
    typedef std::multimap<long int, int> DegreeHeap;
    //typedef boost::unordered_multimap<long int, int> DegreeHeap;

    // Computes a minimum degree ordering for a block matrix given the
    // set of interactions between blocks, and block sizes.
    //
    // blockInteractions will be modified
    static void BlockMinimumDegree(
                      const std::vector<std::set<int> > &blockInteractions,
                      const IntArray &blockSizes,
                      IntArray &permutation, IntArray &inversePermutation );

    // Prints out usage statistics comparing the number of non-zeros
    // in unpermuted and permuted versions of a block system
    static void PrintBlockOrderingStats(
                      std::vector<std::set<int> > &blockInteractions,
                      const IntArray &blockSizes,
                      const IntArray &permutation,
                      const IntArray &inversePermutation );

  private:
    // Given a lower triangular interaction set, compute the symmetric version
    // (effectively a graph adjacency matrix)
    static void SymmetrizeInteractionSets(
                      std::vector<std::set<int> > &blockInteractions );

    // Initializes the degree heap
    static void InitializeHeap(
                      std::vector<std::set<int> > &blockInteractions,
                      DegreeHeap &degreeHeap,
                      std::vector<long int> &blockDegrees,
                      const IntArray &blockSizes );

    // Eliminates a node from the given system and generates the required
    // fill-in.  Modifies the priority queue accordingly.
    static void EliminateBlock(
                      std::vector<std::set<int> > &blockInteractions,
                      DegreeHeap &degreeHeap,
                      std::vector<long int> &blockDegrees,
                      const IntArray &blockSizes,
                      int block_idx );

    // Computes a degree for a given block, given its interactions
    static long int BlockDegree(
                      std::vector<std::set<int> > &blockInteractions,
                      const IntArray &blockSizes,
                      int block_idx );

    // Computes fill-in for a block matrix and returns the number of
    // non-zeros needed to represend the factor
    static long int BlockFillIn(
                      std::vector<std::set<int> > &blockInteractions,
                      const IntArray &blockSizes );

    // Heap-like data structure operations implemented using a multimap

    // Insert
    static void HeapInsert(
                      DegreeHeap &degreeMap,
                      std::pair<long int, int> &entry );

    static void HeapInsert(
                      DegreeHeap &degreeMap,
                      long int degree, int node_idx )
    {
      std::pair<long int, int> toInsert( degree, node_idx );

      HeapInsert( degreeMap, toInsert );
    }

    // Delete
    static void HeapDelete(
                      DegreeHeap &degreeMap,
                      std::pair<long int, int> &entry );

    static void HeapDelete(
                      DegreeHeap &degreeMap,
                      long int degree, int node_idx )
    {
      std::pair<long int, int> toDelete( degree, node_idx );

      HeapDelete( degreeMap, toDelete );
    }

	private:
		MinimumDegree();

		// Destructor
		virtual ~MinimumDegree();

};

#endif
