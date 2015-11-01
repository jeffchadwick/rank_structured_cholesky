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
// MeshReader.cpp - Definition of the MeshReader class
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
