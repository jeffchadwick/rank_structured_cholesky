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

// FaceMap.h - Definition for the FaceMap class

#ifndef FACEMAP_H
#define FACEMAP_H

#include <rschol/SETTINGS.h>
#include <rschol/TYPES.h>

#include <map>
#include <set>

class LinearElement;

// FaceMap class
// Provides a data structure used to determine connectivity
// between elements.  Also helps with identifying boundary
// faces
class FaceMap
{
	public:
		// Helper class to represent a pair of adjacent
		// elements
		class ElementPair {
			public:
				ElementPair()
					: e1(NULL), e2(NULL), faceNum1(0), faceNum2(0),
						index1(-1), index2(-1)
				{
				}

				ElementPair(const LinearElement *e1, int faceNum1, int index1)
					: e1(e1), e2(NULL), faceNum1(faceNum1), faceNum2(0),
						index1(index1), index2(-1)
				{
				}

				const LinearElement *e1;
				const LinearElement *e2;

				// Face numbers from each element
				int faceNum1;
				int faceNum2;

				// Element numbers for each element in some
				// external array
				int index1;
				int index2;
		};

		// We can support faces with 3 or 4 vertices
		FaceMap(int numVertices, int verticesPerFace);

		// Insert a face in to the face map.
		void insertFace(int v1, int v2, int v3, const LinearElement *e,
										int faceNum, int elementIndex);
		void insertFace(int v1, int v2, int v3, int v4, const LinearElement *e,
										int faceNum, int elementIndex);

		// Get the Element pair associated with an edge
		ElementPair lookupFace(int v1, int v2, int v3);
		ElementPair lookupFace(int v1, int v2, int v3, int v4);

		// A list of boundary faces
		std::set<ElementPair *> &boundaryFaces() { return _boundaryFaces; }

	private:
		// Generate a hash value for a vertex index pair
		long long hashValue(int v1, int v2, int v3) const;
		long long hashValue(int v1, int v2, int v3, int v4) const;

		// Store a map of long ints to element pairs
		std::map<long long, ElementPair> _table;

		// Need to know the number of vertices and
		// various powers of this to generate hash values
		int _numVertices;
		int _numVertices2;
		int _numVertices3;

		int _verticesPerFace;

		// A set of boundary pairs (faces with only one
		// adjacent element)
		std::set<ElementPair *> _boundaryFaces;
};

#endif
