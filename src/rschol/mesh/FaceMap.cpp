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

// FaceMap.cpp - Definition for the FaceMap class

#include "FaceMap.h"
#include "LinearElement.h"

#include <algorithm>
#include <vector>

using namespace std;

// Class constructor
FaceMap::FaceMap(int numVertices, int verticesPerFace)
	: _numVertices(numVertices), _verticesPerFace(verticesPerFace)
{
	_numVertices2 = _numVertices * _numVertices;
	_numVertices3 = _numVertices2 * _numVertices;
}

// Insert a face in to the map
void FaceMap::insertFace(int v1, int v2, int v3, const LinearElement *e,
												 int faceNum, int elementIndex)
{
	long long hash = hashValue( v1, v2, v3 );

	// Check to see if there is alread an entry in the hash
	// table for this face, and put one there if not
	if ( _table.find(hash) != _table.end() )
	{
		// Insert this edge in to the existing pair
		FaceMap::ElementPair &pair = _table[hash];
		FaceMap::ElementPair *pairPtr = &pair;

		if (pair.e1 != e)
		{
			pair.e2 = e;
			pair.faceNum2 = faceNum;

			// No longer on the boundary, so take it
			// out of the boundary set
			_boundaryFaces.erase( pairPtr );
		}
	}
	else
	{
		// Insert a new pair with this element in it
		_table[hash] = FaceMap::ElementPair(e, faceNum, elementIndex);

		FaceMap::ElementPair &pair = _table[hash];
		FaceMap::ElementPair *pairPtr = &pair;

		// Since this face currently has only one adjacent element,
		// it should be in the boundary set.
		_boundaryFaces.insert( pairPtr );
	}
}

void FaceMap::insertFace(int v1, int v2, int v3, int v4,
												 const LinearElement *e, int faceNum,
												 int elementIndex)
{
	long long hash = hashValue( v1, v2, v3, v4 );

	// Check to see if there is alread an entry in the hash
	// table for this face, and put one there if not
	if ( _table.find(hash) != _table.end() )
	{
		// Insert this edge in to the existing pair
		FaceMap::ElementPair &pair = _table[hash];
		FaceMap::ElementPair *pairPtr = &pair;

		if (pair.e1 != e)
		{
			pair.e2 = e;
			pair.faceNum2 = faceNum;

			// No longer on the boundary, so take it
			// out of the boundary set
			_boundaryFaces.erase( pairPtr );
		}
	}
	else
	{
		// Insert a new pair with this element in it
		_table[hash] = FaceMap::ElementPair(e, faceNum, elementIndex);

		FaceMap::ElementPair &pair = _table[hash];
		FaceMap::ElementPair *pairPtr = &pair;

		// Since this face currently has only one adjacent element,
		// it should be in the boundary set.
		_boundaryFaces.insert( pairPtr );
	}
}

// Lookup a face
FaceMap::ElementPair FaceMap::lookupFace(int v1, int v2, int v3)
{
	// Initialize a null pair
	FaceMap::ElementPair pair;

	long long hash = hashValue(v1, v2, v3);

	// If we can't find this in the table, return a
	// null pair
	if (_table.find(hash) == _table.end())
		return pair;

	// Otherwise, return the pair associated with this face
	pair = _table[hash];

	return pair;
}

FaceMap::ElementPair FaceMap::lookupFace(int v1, int v2, int v3, int v4)
{
	// Initialize a null pair
	FaceMap::ElementPair pair;

	long long hash = hashValue(v1, v2, v3, v4);

	// If we can't find this in the table, return a
	// null pair
	if (_table.find(hash) == _table.end())
		return pair;

	// Otherwise, return the pair associated with this face
	pair = _table[hash];

	return pair;
}

// Generate a hash value
long long FaceMap::hashValue(int v1, int v2, int v3) const
{
	// We will use the convention that the low-order "bits"
	// of the hash value correspond to the lowest indeces.
	long long hashValue;

	vector<int> indices;
	indices.push_back(v1);
	indices.push_back(v2);
	indices.push_back(v3);

	sort(indices.begin(), indices.end());

	hashValue = indices.at(0) + _numVertices * indices.at(1)
							+ _numVertices2 * indices.at(2);

	return hashValue;
}

long long FaceMap::hashValue(int v1, int v2, int v3, int v4) const
{
	// We will use the convention that the low-order "bits"
	// of the hash value correspond to the lowest indeces.
	long long hashValue;

	vector<int> indices;
	indices.push_back(v1);
	indices.push_back(v2);
	indices.push_back(v3);
	indices.push_back(v4);

	sort(indices.begin(), indices.end());

	hashValue = indices.at(0) + _numVertices * indices.at(1)
							+ _numVertices2 * indices.at(2)
							+ _numVertices3 * indices.at(3);

	return hashValue;
}
