//////////////////////////////////////////////////////////////////////
// Mesh.cpp - Definitions for the mesh class
// Jeffrey Chadwick (jnc52)
//
//////////////////////////////////////////////////////////////////////

#include "Mesh.h"

#include "FaceMap.h"
#include "LinearElement.h"
#include "IsoHexElement.h"

#include <set>

using namespace std;

//////////////////////////////////////////////////////////////////////
// Construct a new mesh
//////////////////////////////////////////////////////////////////////
Mesh::Mesh( const Vector3Array &vertices, Mesh::ElementType type,
						const IntArray &connectivity,
						const vector<bool> &constrained,
						Real E, Real v, Real density,
						int numGaussPoints, const Tuple3i &divisions )
	: _restPose( vertices ),
    _type( type ),
    _faceMap( NULL ),
		_E( E ),
    _v( v ),
    _density( density ),
    _numGaussPoints( numGaussPoints ),
    _divisions( divisions )
{
	IntArray newConnectivity( connectivity );

  cout << "   Mesh constructor: reindex constraints" << endl;
	// Re-index to respect constraints
	reindexConstraints( constrained, newConnectivity );

  cout << "   Mesh constructor: build elements" << endl;
	buildElements( newConnectivity );

  cout << "   Mesh constructor: build boundary" << endl;
	buildBoundary();

  cout << "   Mesh constructor: done" << endl;
}

//////////////////////////////////////////////////////////////////////
// Assembles mass and stiffness matrices for this mesh
//////////////////////////////////////////////////////////////////////
void Mesh::assembleMatrices()
{
	assembleMass();
	assembleStiffness();

	_internalForces.resizeAndWipe( 3 * _numUnconstrained );
}

//////////////////////////////////////////////////////////////////////
// Internal force computation
//////////////////////////////////////////////////////////////////////
void Mesh::computeInternalForces( VECTOR &meshState )
{
	_internalForces.clear();

	int localIndex, globalIndex;

	// Iterate over each element and assemble the internal
	// force vector
	for (int i = 0; i < _elements.size(); i++)
	{
		VECTOR &localForces
			= _elements.at(i)->computeLocalInternalForces( meshState );
		IntArray &localMap
			= _elements.at(i)->localMap();

		for (int k = 0; k < localMap.size(); k++)
		{
			localIndex = 3 * k;
			
			if ( localMap.at(k) < 0 || localMap.at(k) >= _numUnconstrained )
				continue;

			globalIndex = 3 * localMap.at(k);

			// Add the 3-vector corresponding to this node
			_internalForces( globalIndex     ) += localForces( localIndex     );
			_internalForces( globalIndex + 1 ) += localForces( localIndex + 1 );
			_internalForces( globalIndex + 2 ) += localForces( localIndex + 2 );
		}
	}
}

//////////////////////////////////////////////////////////////////////
// Determine strain at a triangle on the boundary
//////////////////////////////////////////////////////////////////////
VECTOR Mesh::boundaryStrain( VECTOR &meshState, int idx )
{
	return _boundaryElements.at(idx)->getFaceStrain( meshState,
																									 _boundaryFaces.at(idx) );
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
VECTOR Mesh::elementStrain( VECTOR &meshState, int idx )
{
	return _elements.at(idx)->getStrain( meshState );
}

//////////////////////////////////////////////////////////////////////
// Evaluates a vector field at all points of this mesh
//////////////////////////////////////////////////////////////////////
void Mesh::evaluateVectorField( Vec3Evaluator &evaluator,
                                Vector3Array &results ) const
{
  if ( results.size() < _numUnconstrained )
  {
    results.resize( _numUnconstrained );
  }

  for ( int point_idx = 0; point_idx < _numUnconstrained; point_idx++ )
  {
    const VEC3F             &x = _restPose[ point_idx ];

    results[ point_idx ] = evaluator( x );
  }
}

//////////////////////////////////////////////////////////////////////
// Evalutes a scalar field at all points of this mesh
//////////////////////////////////////////////////////////////////////
void Mesh::evaluateScalarField( ScalarEvaluator &evaluator,
                                FloatArray &results ) const
{
  if ( results.size() < _numUnconstrained )
  {
    results.resize( _numUnconstrained );
  }

  for ( int point_idx = 0; point_idx < _numUnconstrained; point_idx++ )
  {
    const VEC3F             &x = _restPose[ point_idx ];

    results[ point_idx ] = evaluator( x );
  }
}

//////////////////////////////////////////////////////////////////////
// Stiffness assembly
//////////////////////////////////////////////////////////////////////
void Mesh::assembleStiffness()
{
	int localIndex1, localIndex2, globalIndex1, globalIndex2;

	_stiffnessMatrix = SPARSE_MATRIX( 3 * _numUnconstrained,
																		3 * _numUnconstrained );

	_stiffnessDiagonal.resizeAndWipe( 3 * _numUnconstrained );

	for (int i = 0; i < _elements.size(); i++)
	{
		// The local to global mapping for this element
		const IntArray &localMap = _elements.at(i)->localMap();

		MATRIX &localStiffness = _elements.at(i)->localStiffness();

		// Iterate over each entry in the stiffness matrix
		for (int k1 = 0; k1 < localMap.size(); k1++)
		{
			// No entries in the matrix for constrained nodes
			if (localMap.at(k1) < 0 || localMap.at(k1) >= _numUnconstrained)
				continue;

			localIndex1 = 3 * k1;
			globalIndex1 = 3 * localMap.at(k1);

			for (int k2 = 0; k2 < localMap.size(); k2++)
			{
				if (localMap.at(k2) < 0 || localMap.at(k2) >= _numUnconstrained)
					continue;

				localIndex2 = 3 * k2;
				globalIndex2 = 3 * localMap.at(k2);

				// Fill in the 3x3 submatrix for this vertex
				_stiffnessMatrix(globalIndex1, globalIndex2)
					+= localStiffness(localIndex1, localIndex2);
				_stiffnessMatrix(globalIndex1, globalIndex2 + 1)
					+= localStiffness(localIndex1, localIndex2 + 1);
				_stiffnessMatrix(globalIndex1, globalIndex2 + 2)
					+= localStiffness(localIndex1, localIndex2 + 2);

				_stiffnessMatrix(globalIndex1 + 1, globalIndex2)
					+= localStiffness(localIndex1 + 1, localIndex2);
				_stiffnessMatrix(globalIndex1 + 1, globalIndex2 + 1)
					+= localStiffness(localIndex1 + 1, localIndex2 + 1);
				_stiffnessMatrix(globalIndex1 + 1, globalIndex2 + 2)
					+= localStiffness(localIndex1 + 1, localIndex2 + 2);

				_stiffnessMatrix(globalIndex1 + 2, globalIndex2)
					+= localStiffness(localIndex1 + 2, localIndex2);
				_stiffnessMatrix(globalIndex1 + 2, globalIndex2 + 1)
					+= localStiffness(localIndex1 + 2, localIndex2 + 1);
				_stiffnessMatrix(globalIndex1 + 2, globalIndex2 + 2)
					+= localStiffness(localIndex1 + 2, localIndex2 + 2);

				// Add to the row sum for the lumped stiffness
				_stiffnessDiagonal(globalIndex1)
					+= localStiffness(localIndex1, localIndex2);
				_stiffnessDiagonal(globalIndex1)
					+= localStiffness(localIndex1, localIndex2 + 1);
				_stiffnessDiagonal(globalIndex1)
					+= localStiffness(localIndex1, localIndex2 + 2);

				_stiffnessDiagonal(globalIndex1 + 1)
					+= localStiffness(localIndex1 + 1, localIndex2);
				_stiffnessDiagonal(globalIndex1 + 1)
					+= localStiffness(localIndex1 + 1, localIndex2 + 1);
				_stiffnessDiagonal(globalIndex1 + 1)
					+= localStiffness(localIndex1 + 1, localIndex2 + 2);

				_stiffnessDiagonal(globalIndex1 + 2)
					+= localStiffness(localIndex1 + 2, localIndex2);
				_stiffnessDiagonal(globalIndex1 + 2)
					+= localStiffness(localIndex1 + 2, localIndex2 + 1);
				_stiffnessDiagonal(globalIndex1 + 2)
					+= localStiffness(localIndex1 + 2, localIndex2 + 2);
			}
		}
	}
}

//////////////////////////////////////////////////////////////////////
// Mass matrix assembly
//////////////////////////////////////////////////////////////////////
void Mesh::assembleMass()
{
	int localIndex1, localIndex2, globalIndex1, globalIndex2;

	_massMatrix = SPARSE_MATRIX( 3 * _numUnconstrained,
															 3 * _numUnconstrained );

	_massDiagonal.resizeAndWipe( 3 * _numUnconstrained );

	for (int i = 0; i < _elements.size(); i++)
	{
		// The local to global mapping for this element
		const IntArray &localMap = _elements.at(i)->localMap();

		MATRIX &localMass = _elements.at(i)->localMass();

		// Iterate over each entry in the mass matrix
		for (int k1 = 0; k1 < localMap.size(); k1++)
		{
			// No entries in the matrix for constrained nodes
			if (localMap.at(k1) < 0 || localMap.at(k1) >= _numUnconstrained)
				continue;

			localIndex1 = k1;
			globalIndex1 = 3 * localMap.at(k1);

			for (int k2 = 0; k2 < localMap.size(); k2++)
			{
				if (localMap.at(k2) < 0 || localMap.at(k2) >= _numUnconstrained)
					continue;

				localIndex2 = k2;
				globalIndex2 = 3 * localMap.at(k2);

				// Add a diagonal component for each node
				_massMatrix(globalIndex1, globalIndex2)
					+= localMass(localIndex1, localIndex2);
				_massMatrix(globalIndex1 + 1, globalIndex2 + 1)
					+= localMass(localIndex1, localIndex2);
				_massMatrix(globalIndex1 + 2, globalIndex2 + 2)
					+= localMass(localIndex1, localIndex2);

				// Add to the row sum for the lumped stiffness
				_massDiagonal(globalIndex1)
					+= localMass(localIndex1, localIndex2);
				_massDiagonal(globalIndex1 + 1)
					+= localMass(localIndex1, localIndex2);
				_massDiagonal(globalIndex1 + 2)
					+= localMass(localIndex1, localIndex2);
			}
		}
	}
}

//////////////////////////////////////////////////////////////////////
// Builds all elements for this mesh, depending on the type
// of element used
//////////////////////////////////////////////////////////////////////
void Mesh::buildElements( const IntArray &connectivity )
{
	switch ( _type )
	{
		case Mesh::ISO_8_NODE:
		{
			// Initialize a face map with four elements per face
			_faceMap = new FaceMap( _restPose.size(), 8 );

			if ( connectivity.size() % 8 != 0 )
				cerr << "Connectivity information is wrong for 8 node element" << endl;

			int elementIndex = 0;

			// Start building elements
			for (int i = 0; i < connectivity.size(); i += 8)
			{
				_elements.push_back( new IsoHexElement( this, elementIndex,
																							  _E, _v, _density,
																							  connectivity.at(i),
																							  connectivity.at(i+1),
																							  connectivity.at(i+2),
																							  connectivity.at(i+3),
																							  connectivity.at(i+4),
																							  connectivity.at(i+5),
																							  connectivity.at(i+6),
																							  connectivity.at(i+7),
																							  _numGaussPoints ) );

				// Precompute local stiffness, mass, etc. for the element
				_elements.back()->precompute();

				elementIndex += 1;
			}

			break;
		}
		default:
		{
			break;
		}
	}

	// Have each element put all of its faces in to the face
	// map so that we can get connectivity information.
	for (int i = 0; i < _elements.size(); i++)
	{
		_elements.at(i)->addFaces( _faceMap );
	}
}

//////////////////////////////////////////////////////////////////////
// Find all boundary triangles
//////////////////////////////////////////////////////////////////////
void Mesh::buildBoundary()
{
	set<FaceMap::ElementPair *> &boundaryFaces = _faceMap->boundaryFaces();

	set<FaceMap::ElementPair *>::iterator i = boundaryFaces.begin();

	_boundary.clear();
	_boundaryElements.clear();
	_boundaryFaces.clear();

	for ( ; i != boundaryFaces.end(); i++)
	{
		// Get visible faces from each boundary
		FaceMap::ElementPair *pair = *i;

		if ( pair->e2 != NULL )
		{
			cout << "ERROR: Invalid boundary face" << endl;
		}

		IntArray faceIndices;

		// Get triangles from this face
		pair->e1->getFaceTriangles( pair->faceNum1, faceIndices );

		if ( faceIndices.size() % 3 != 0 )
		{
			cout << "ERROR: Invalid face triangles" << endl;
		}

		for (int j = 0; j < faceIndices.size(); j++)
		{
			_boundary.push_back( faceIndices.at(j) );
		}

		for (int j = 0; j < faceIndices.size() / 3; j++)
		{
			_boundaryElements.push_back( const_cast<LinearElement *>(pair->e1) );
			_boundaryFaces.push_back( pair->faceNum1 );
		}
	}
}

//////////////////////////////////////////////////////////////////////
// Reindex vertices to deal with constraints
//////////////////////////////////////////////////////////////////////
void Mesh::reindexConstraints( const std::vector<bool> &constrained,
															 IntArray &connectivity )
{
	// Start by counting the number of unconstrained nodes
	_numUnconstrained = _restPose.size();
	for (int i = 0; i < constrained.size(); i++)
	{
		if (constrained.at(i))
		{
			_numUnconstrained -= 1;
		}
	}

	int numConstrained = _restPose.size() - _numUnconstrained;

	// Initialize a vector of new indeces for each of our vertices
	IntArray newIndeces(_restPose.size());

	// We'll also need a new rest pose
	Vector3Array newRestPose(_restPose.size());

	// Start reindexing constrained and uncontrained nodes.
	int unconstrainedIndex = 0;
	int constrainedIndex = _restPose.size() - numConstrained;

	for (int i = 0; i < _restPose.size(); i++)
	{
		if (constrained.at(i) == true)
		{
			// Put this node near the end of the vertex array
			newIndeces.at(i) = constrainedIndex;

			newRestPose.at(constrainedIndex) = _restPose.at(i);

			constrainedIndex++;
		}
		else
		{
			// Put this node near the beginning of the vertex array
			newIndeces.at(i) = unconstrainedIndex;

			newRestPose.at(unconstrainedIndex) = _restPose.at(i);

			unconstrainedIndex++;
		}
	}

	// Rearrange the triangle index array using new indeces
	for (int i = 0; i < connectivity.size(); i++)
	{
		int newIndex = newIndeces.at(connectivity.at(i));

		connectivity.at(i) = newIndex;
	}

	// Assign the new rest pose
	_restPose = newRestPose;
}
