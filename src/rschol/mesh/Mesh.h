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
// Mesh.h - Definition for the Mesh class
//
//////////////////////////////////////////////////////////////////////

#ifndef MESH_H
#define MESH_H

#include <rschol/SETTINGS.h>
#include <rschol/TYPES.h>

#include <rschol/geometry/BoundingBox.h>

#include <rschol/linearalgebra/VECTOR.h>
#include <rschol/linearalgebra/SPARSE_MATRIX.h>

#include <rschol/util/Evaluator.h>

#include <vector>

class FaceMap;
class LinearElement;

//////////////////////////////////////////////////////////////////////
// Mesh class
// Defines a finite element mesh
//////////////////////////////////////////////////////////////////////
class Mesh
{
	public:
    class MeshPoint {
      public:
        MeshPoint( int index, const VEC3F &x )
          : _index( index )
        {
          _bbox.setMinBound( x );
          _bbox.setMaxBound( x );
        }

        int index() const { return _index; }

        const BoundingBox &bbox() const { return _bbox; }

      private:
        int                  _index;

        // For KD tree compatibility
        BoundingBox          _bbox;

    };

		enum ElementType {
			ISO_8_NODE = 0,
			NUM_ELEMENT_TYPES
		};

		// Constructs a new mesh, supplying a list of vertices,
		// the type of element to be used, and connectivity information
		// for the elements
		Mesh( const Vector3Array &vertices, ElementType type,
					const IntArray &connectivity,
					const std::vector<bool> &constrained,
					Real E, Real v, Real density,
					int numGaussPoints, const Tuple3i &divisions );

		// Assembles the global mass and stiffness matrices
		void assembleMatrices();

		// The mesh rest pose
		Vector3Array &restPose() { return _restPose; }
		const Vector3Array &restPose() const { return _restPose; }

		VECTOR &internalForces() { return _internalForces; }

		// Accessors
		SPARSE_MATRIX &stiffnessMatrix() { return _stiffnessMatrix; }
		SPARSE_MATRIX &massMatrix() { return _massMatrix; }
		VECTOR &stiffnessDiagonal() { return _stiffnessDiagonal; }
		VECTOR &massDiagonal() { return _massDiagonal; }

		// Compute internal forces for the mesh, given
		// displacements and velocities.  This includes
		// elastic and damping forces.
		// alpha and beta are Raleigh damping parameters.
		void computeInternalForces( VECTOR &meshState, VECTOR &meshVelocities,
																Real alpha, Real beta );

		// Computes only elasticity forces
		void computeInternalForces( VECTOR &meshState );

		int numUnconstrained() const { return _numUnconstrained; }

		int numElements() const { return _elements.size(); }

		// Boundary triangles
		IntArray &boundary() { return _boundary; }

		// Returns the strain at the given boundary triangle
		VECTOR boundaryStrain( VECTOR &meshState, int idx );

		// Evaluates strain somewhere in the element
		VECTOR elementStrain( VECTOR &meshState, int idx );

    // Evaluates a vector field at all points of this mesh
    void evaluateVectorField( Vec3Evaluator &evaluator,
                              Vector3Array &results ) const;

    // Evalutes a scalar field at all points of this mesh
    void evaluateScalarField( ScalarEvaluator &evaluator,
                              FloatArray &results ) const;

    const Tuple3i &divisions() const { return _divisions; }

	private:
		// Take a boolean array indicating which points are constrained
		// and perform reindexing to place all unconstrained nodes at
		// the bottom of the vertex list and all constrained nodes
		// at the top.  This also requires that we rearrange the
		// connectivity index array.
		void reindexConstraints( const std::vector<bool> &constrained,
														 IntArray &connectivity );

		// Assembles the global stiffness matrix, and its
		// diagonally lumped version.
		void assembleStiffness();

		// Assembles the global mass matrix, and its diagonally
		// lumped version.
		void assembleMass();

		void buildElements( const IntArray &connectivity );

		// Get boundary triangles
		void buildBoundary();

	private:
		// Mesh vertices
		Vector3Array _restPose;

		// Mesh elements
		std::vector<LinearElement *> _elements;

		ElementType _type;

		// Diagonal of the mass matrix, assuming that
		// mass lumping is in effect
		VECTOR _massDiagonal;

		// Diagonally lumped stiffness matrix, which
		// is useful for damping purposes
		VECTOR _stiffnessDiagonal;

		// A vector used to store internal forces
		VECTOR _internalForces;

		// The global stiffness matrix
		SPARSE_MATRIX _stiffnessMatrix;

		SPARSE_MATRIX _massMatrix;

		// The number of unconstrained nodes in the mesh
		int _numUnconstrained;

		FaceMap *_faceMap;

		// Physical parameters
		Real _E, _v, _density;

		// Number of integration points, if needed
		int _numGaussPoints;

		// Boundary triangle indices
		IntArray _boundary;
		std::vector<LinearElement *> _boundaryElements;
		IntArray _boundaryFaces;

    Tuple3i _divisions;
		
};

#endif
