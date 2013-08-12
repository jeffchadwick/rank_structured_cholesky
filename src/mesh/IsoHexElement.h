// IsoHexElement.h - Definition for the IsoHexElement class
// Jeffrey Chadwick (jnc52)

#ifndef ISOHEXELEMENT_H
#define ISOHEXELEMENT_H

#include <SETTINGS.h>
#include <TYPES.h>

#include "LinearElement.h"

class Mesh;

// IsoHexElement class
// Represents and Isoparametric, 8 node hexahedral element
class IsoHexElement : public LinearElement
{
	public:
		// Constructs a new element given a mesh and
		// the indices of the vertices included in this
		// element.  Also needed are the elastic modulus,
		// Poisson ratio and density.
		IsoHexElement( Mesh *mesh, int index,
									 Real E, Real v, Real density,
									 int v0, int v1, int v2, int v3,
									 int v4, int v5, int v6, int v7,
									 int numGaussPoints );

		virtual ~IsoHexElement()
		{
		}

		// Adds this element's faces to a face map data structure
		virtual void addFaces( FaceMap *faceMap );

		// Given a face number, provide an array of indices that can
		// be used to draw triangles representing the given face
		virtual void getFaceTriangles( int faceNum, IntArray &indices ) const;

		// Returns the strain experienced at the center of the given
		// element face
		virtual VECTOR &getFaceStrain( VECTOR &meshState, int faceNum );

		// Return strain at the center of the element
		virtual VECTOR &getStrain( VECTOR &meshState );

	protected:
		// Virtual pre-computation functions
		virtual void precomputeStiffness();
		virtual void precomputeMass();

	private:
		// Adds the stiffness integrand, evaluated at the given
		// coordinates and scaled by the given weight, to the
		// stiffness matrix
		void addStiffnessIntegrandEvaluation( VEC3F &p, Real weight );

		// Adds the mass integrand, evaluated at the given
		// coordinates and scaled by the given weight, to the
		// mass matrix
		void addMassIntegrandEvaluation( VEC3F &p, Real weight );

		// Evaluates strain at the given xi, eta, zeta coordinates
		void evaluateStrain( VEC3F &p, VECTOR &meshState );

	private:
		Real E, v, density;

		int vertices[8];

		// Store x y and z values explicitly to sync up
		// more easily with generated code.
		Real x0, x1, x2, x3, x4, x5, x6, x7;
		Real y0, y1, y2, y3, y4, y5, y6, y7;
		Real z0, z1, z2, z3, z4, z5, z6, z7;

		// The number of quadrature points to
		// use when evaluating integrals
		int _numGaussPoints;
};

#endif
