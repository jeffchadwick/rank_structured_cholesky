//////////////////////////////////////////////////////////////////////
// MeshReader.cpp - Definition of the MeshReader class
// Jeffrey Chadwick (jnc52)
//
//////////////////////////////////////////////////////////////////////

#include "Mesh.h"
#include "MeshReader.h"

#include <rschol/TYPES.h>
#include <rschol/SETTINGS.h>

#include <rschol/linearalgebra/VEC3.h>

#include <vector>

#include <iostream>
#include <fstream>

using namespace std;

//////////////////////////////////////////////////////////////////////
// Read a mesh
//////////////////////////////////////////////////////////////////////
Mesh *MeshReader::BuildMesh( const string &inputFile )
{
	// Open the input file
	ifstream f( inputFile.c_str() );

	Real E, v, density;
	int numGaussPoints;
	int elementType;
	int numVertices, numElements, numConstrained;

  Tuple3i divisions;

	f >> E;
	f >> v;
	f >> density;
	f >> numGaussPoints;
	f >> elementType;
	f >> numVertices;
	f >> numElements;
	f >> numConstrained;
  f >> divisions[ 0 ];
  f >> divisions[ 1 ];
  f >> divisions[ 2 ];

	Vector3Array vertices;
	for ( int i = 0; i < numVertices; i++ )
	{
		VEC3F p;

		f >> p[0];
		f >> p[1];
		f >> p[2];

		vertices.push_back( p );
	}

	// Read in elements
	if ( elementType < 0 || elementType >= Mesh::NUM_ELEMENT_TYPES )
	{
		cout << "Invalid element type!" << endl;
		return NULL;
	}

	IntArray connectivity;

	if ( elementType == Mesh::ISO_8_NODE )
	{
		for ( int i = 0; i < numElements; i++ )
		{
			for ( int j = 0; j < 8; j++ )
			{
				int nodeNum;
				f >> nodeNum;
				connectivity.push_back( nodeNum );
			}
		}
	}

	vector<bool> constrained(vertices.size());

	for (int i = 0; i < constrained.size(); i++)
	{
		constrained.at(i) = false;
	}

	// Read in constrained nodes
	for ( int i = 0; i < numConstrained; i++ )
	{
		int nodeNum;
		f >> nodeNum;
		constrained.at(nodeNum) = true;
	}

	f.close();

	// Finish building the mesh
	return new Mesh( vertices, (Mesh::ElementType)elementType, connectivity,
									 constrained, E, v, density, numGaussPoints, divisions );
}
