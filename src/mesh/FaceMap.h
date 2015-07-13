// FaceMap.h - Definition for the FaceMap class
// Jeffrey Chadwick (jnc52)

#ifndef FACEMAP_H
#define FACEMAP_H

#include <SETTINGS.h>
#include <TYPES.h>

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
