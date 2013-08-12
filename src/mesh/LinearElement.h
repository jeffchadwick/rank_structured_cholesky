// LinearElement.h - Header for the LinearElement class
// Jeffrey Chadwick (jnc52)

#ifndef LINEARELEMENT_H
#define LINEARELEMENT_H

#include <SETTINGS.h>
#include <TYPES.h>

#include <linearalgebra/MATRIX.h>
#include <linearalgebra/VECTOR.h>
#include <linearalgebra/VEC3.h>

class FaceMap;

// Element class
// Abstract class used to define basic element routines for a
// linear element.
class LinearElement
{
	public:
		LinearElement( int index )
			: _index(index), _strain(6), _stress(6)
		{
		}

		virtual ~LinearElement()
		{
		}

		// Accessors
		MATRIX &localStiffness() { return _localStiffness; }
		MATRIX &localMass() { return _localMass; }
		IntArray &localMap() { return _localMap; }

		// Precompute necessary data (ie. local stiffness and mass matrices)
		// for this element, given the initial mesh state.
		void precompute();

		// Computes local internal forces based on the given
		// mesh state
		VECTOR& computeLocalInternalForces( VECTOR &meshState );

		// Computes stiffness proportional damping forces for this
		// element, given the current mesh velocity.  Beta is the
		// stiffness damping coefficient.
		VECTOR& computeStiffnessProportionalDamping( VECTOR &meshVelocity,
																								 Real beta );

		int index() { return _index; }

		// Adds this element's faces to a face map data structure
		virtual void addFaces( FaceMap *faceMap ) = 0;

		// Given a face number, provide an array of indices that can
		// be used to draw triangles representing the given face
		virtual void getFaceTriangles( int faceNum, IntArray &indices ) const = 0;

		// Returns the strain experienced at the center of the given
		// element face
		virtual VECTOR &getFaceStrain( VECTOR &meshState, int faceNum ) = 0;

		// Return strain somewhere (ie. probably the center) of the
		// element
		virtual VECTOR &getStrain( VECTOR &meshState ) = 0;

	protected:
		// Virtual pre-computation functions
		virtual void precomputeStiffness() = 0;
		virtual void precomputeMass() = 0;

		// Sets the local state, given an overall mesh state
		void setLocalState( VECTOR &meshState );

	protected:
		// The local stiffness matrix, which can be precomputed
		// due to linearity
		MATRIX _localStiffness;

		// Local mass matrix
		MATRIX _localMass;

		// Local map, which maps element node numbers to
		// global node numbers
		IntArray _localMap;

		// Storage for internal forces computed via this
		// element
		VECTOR _internalForces;

		// Used as temporary storage for the degrees of freedom
		// associated with this particular element
		VECTOR _localState;

		// The global element index
		int _index;

		// Temporary storage for stress and strain
		VECTOR _strain, _stress;
};

#endif
