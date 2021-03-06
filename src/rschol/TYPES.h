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
// Types.h: Project-wide type definitions
//
//////////////////////////////////////////////////////////////////////

#ifndef TYPES_H
#define TYPES_H

#include <rschol/SETTINGS.h>

#include <set>
#include <vector>

#include <rschol/linearalgebra/VEC3.h>

#include <boost/multi_array.hpp>
#include <boost/unordered_set.hpp>
#include <boost/unordered_map.hpp>

typedef std::vector<int>		       IntArray;
typedef std::vector<unsigned int>  UIntArray;
typedef std::vector<long int>      LongIntArray;
typedef std::vector<size_t>        SizeArray;
typedef std::vector<bool>		       BoolArray;
typedef std::vector<VEC3F>	       Vector3Array;
typedef std::vector<Real>		       FloatArray;

typedef TVEC3<int>					       Tuple3i;
typedef TVEC3<int>					       Tuple3i;

typedef std::vector<Tuple3i>	Tuple3Array;

typedef std::pair<int,int>    IndexRange;

inline int range_size( const IndexRange &range )
{
  return ( range.second - range.first + 1 );
}

inline bool in_range( const IndexRange &range, int index )
{
  return index >= range.first && index <= range.second;
}

inline bool range_overlap( const IndexRange &range1, const IndexRange &range2 )
{
  return ( range1.first <= range2.second && range2.first <= range1.second );
}

inline IndexRange range_intersection( const IndexRange &range1,
                                      const IndexRange &range2 )
{
  return IndexRange( std::max( range1.first, range2.first ),
                     std::min( range1.second, range2.second ) );
}

typedef std::pair<int,int>    IndexPair;

typedef std::pair<IndexRange, IndexRange> RangePair;

typedef std::vector<IndexPair>        PairArray;
typedef std::vector<IndexRange>       RangeArray;

// 2D and 3D arrays for various data types
typedef boost::multi_array<Real, 2>		ScalarArray2D;
typedef boost::multi_array<Real, 3>		ScalarArray3D;
typedef boost::multi_array<VEC3F, 2>	Vec3Array2D;
typedef boost::multi_array<VEC3F, 3>	Vec3Array3D;
typedef boost::multi_array<bool, 2>		BoolArray2D;
typedef boost::multi_array<bool, 3>		BoolArray3D;
typedef boost::multi_array<int, 2>		IntArray2D;
typedef boost::multi_array<int, 3>		IntArray3D;

template <class T>
struct UnorderedSet {
  typedef boost::unordered_set<T> type;
};

template <class T, class U>
struct UnorderedMap {
  typedef boost::unordered_map<T, U> type;
};

template <class T>
struct OrderedSet {
  typedef std::set<T> type;
};

enum AXIS {
	X_AXIS = 0,
	Y_AXIS,
	Z_AXIS
};

enum BoundaryCondition {
	BC_NONE = 0,
	BC_DIRICHLET,
	BC_NEUMANN
};

// Different factor types
enum ExtendedFactorType {
  CHOLESKY = 0,
  LDL
};

// Different decomposition strategies (in place, or using extended
// variables)
enum CompressionType {
  EXTENDED_VARIABLE = 0,
  IN_PLACE
};

// Whether to form an orthonormal basis for lower triangular blocks
// directly, or to form a basis for their transposes
enum DecompositionBlock {
  DECOMPOSE_STANDARD = 0,
  DECOMPOSE_TRANSPOSE
};

// Different quantities to be represented with low-rank decompositions
enum DecompositionType {
  DECOMPOSE_SCHUR = 0,
  DECOMPOSE_FACTOR
};

// Off-diagonal interactions can either be represented individually
// (i.e. on a per-node basis) or decomposed as a single block
enum OffDiagonalCompressionType {
  COMPRESS_INDIVIDUAL_INTERACTIONS = 0,
  COMPRESS_FULL_OFF_DIAGONAL
};

#endif
