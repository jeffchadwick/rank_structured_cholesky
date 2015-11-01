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
// MinimumDegree.h: Interface for MinimumDegree class
//
//////////////////////////////////////////////////////////////////////

#ifndef MINIMUM_DEGREE_H
#define MINIMUM_DEGREE_H

#include <rschol/linearalgebra/MATRIX.h>
#include <rschol/linearalgebra/MATRIX3.h>
#include <rschol/linearalgebra/VECTOR.h>

#include <rschol/SETTINGS.h>
#include <rschol/TYPES.h>

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
