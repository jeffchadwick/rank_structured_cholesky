// IsoHexElement.cpp - Definition for the IsoHexElement class
// Jeffrey Chadwick (jnc52)

#include "IsoHexElement.h"

#include <math/GaussPoints.h>

#include "FaceMap.h"
#include "Mesh.h"

// Construct a new element
IsoHexElement::IsoHexElement( Mesh *mesh, int index,
															Real E, Real v, Real density,
														  int v0, int v1, int v2, int v3,
														  int v4, int v5, int v6, int v7,
															int numGaussPoints )
	: E(E), v(v), density(density), _numGaussPoints(numGaussPoints),
		LinearElement(index)
{
	// Add to the local map
	_localMap.push_back( v0 );
	_localMap.push_back( v1 );
	_localMap.push_back( v2 );
	_localMap.push_back( v3 );
	_localMap.push_back( v4 );
	_localMap.push_back( v5 );
	_localMap.push_back( v6 );
	_localMap.push_back( v7 );

	vertices[0] = v0;
	vertices[1] = v1;
	vertices[2] = v2;
	vertices[3] = v3;
	vertices[4] = v4;
	vertices[5] = v5;
	vertices[6] = v6;
	vertices[7] = v7;

	// Find initial x,y,z coordinates of all vertices
	const VEC3F &p0 = mesh->restPose().at( v0 );
	x0 = p0[0];
	y0 = p0[1];
	z0 = p0[2];

	const VEC3F &p1 = mesh->restPose().at( v1 );
	x1 = p1[0];
	y1 = p1[1];
	z1 = p1[2];

	const VEC3F &p2 = mesh->restPose().at( v2 );
	x2 = p2[0];
	y2 = p2[1];
	z2 = p2[2];

	const VEC3F &p3 = mesh->restPose().at( v3 );
	x3 = p3[0];
	y3 = p3[1];
	z3 = p3[2];

	const VEC3F &p4 = mesh->restPose().at( v4 );
	x4 = p4[0];
	y4 = p4[1];
	z4 = p4[2];

	const VEC3F &p5 = mesh->restPose().at( v5 );
	x5 = p5[0];
	y5 = p5[1];
	z5 = p5[2];

	const VEC3F &p6 = mesh->restPose().at( v6 );
	x6 = p6[0];
	y6 = p6[1];
	z6 = p6[2];

	const VEC3F &p7 = mesh->restPose().at( v7 );
	x7 = p7[0];
	y7 = p7[1];
	z7 = p7[2];

	_localStiffness.resizeAndWipe( 24, 24 );
	_localMass.resizeAndWipe( 8, 8 );
}

// Adds this element's faces to the given face map
void IsoHexElement::addFaces( FaceMap *faceMap )
{
	// Add the six faces of this element to the map
	faceMap->insertFace( vertices[0], vertices[1], vertices[2], vertices[3],
											 this, 0, _index );
	faceMap->insertFace( vertices[4], vertices[5], vertices[6], vertices[7],
											 this, 1, _index );
	faceMap->insertFace( vertices[0], vertices[1], vertices[5], vertices[4],
											 this, 2, _index );
	faceMap->insertFace( vertices[2], vertices[3], vertices[7], vertices[6],
											 this, 3, _index );
	faceMap->insertFace( vertices[0], vertices[3], vertices[7], vertices[4],
											 this, 4, _index );
	faceMap->insertFace( vertices[1], vertices[2], vertices[6], vertices[5],
											 this, 5, _index );
}

void IsoHexElement::getFaceTriangles( int faceNum, IntArray &indices ) const
{
	indices.clear();

	switch ( faceNum )
	{
		case 0:
		{
			indices.push_back( vertices[0] );
			indices.push_back( vertices[2] );
			indices.push_back( vertices[1] );

			indices.push_back( vertices[2] );
			indices.push_back( vertices[0] );
			indices.push_back( vertices[3] );

			break;
		}
		case 1:
		{
			indices.push_back( vertices[4] );
			indices.push_back( vertices[5] );
			indices.push_back( vertices[6] );

			indices.push_back( vertices[4] );
			indices.push_back( vertices[6] );
			indices.push_back( vertices[7] );

			break;
		}
		case 2:
		{
			indices.push_back( vertices[5] );
			indices.push_back( vertices[0] );
			indices.push_back( vertices[1] );

			indices.push_back( vertices[4] );
			indices.push_back( vertices[0] );
			indices.push_back( vertices[5] );

			break;
		}
		case 3:
		{
			indices.push_back( vertices[3] );
			indices.push_back( vertices[6] );
			indices.push_back( vertices[2] );

			indices.push_back( vertices[6] );
			indices.push_back( vertices[3] );
			indices.push_back( vertices[7] );

			break;
		}
		case 4:
		{
			indices.push_back( vertices[4] );
			indices.push_back( vertices[3] );
			indices.push_back( vertices[0] );

			indices.push_back( vertices[7] );
			indices.push_back( vertices[3] );
			indices.push_back( vertices[4] );

			break;
		}
		case 5:
		{
			indices.push_back( vertices[2] );
			indices.push_back( vertices[5] );
			indices.push_back( vertices[1] );

			indices.push_back( vertices[5] );
			indices.push_back( vertices[2] );
			indices.push_back( vertices[6] );

			break;
		}
		default:
		{
			break;
		}
	}
}

// Evaluate strain at the center of the element
VECTOR &IsoHexElement::getStrain( VECTOR &meshState )
{
	VEC3F p( 0, 0, 0 );
	evaluateStrain( p, meshState );

	return _strain;
}

// Get the strain exprerienced at the center of the
// given face
VECTOR &IsoHexElement::getFaceStrain( VECTOR &meshState, int faceNum )
{
	switch ( faceNum )
	{
		case 0:
		{
			VEC3F p( 0, 0, -1 );

			evaluateStrain( p, meshState );
			break;
		}
		case 1:
		{
			VEC3F p( 0, 0, 1 );

			evaluateStrain( p, meshState );
			break;
		}
		case 2:
		{
			VEC3F p( 0, -1, 0 );

			evaluateStrain( p, meshState );
			break;
		}
		case 3:
		{
			VEC3F p( 0, 1, 0 );

			evaluateStrain( p, meshState );
			break;
		}
		case 4:
		{
			VEC3F p( -1, 0, 0 );

			evaluateStrain( p, meshState );
			break;
		}
		case 5:
		{
			VEC3F p( 1, 0, 0 );

			evaluateStrain( p, meshState );
			break;
		}
	}

	return _strain;
}

// Precomputation functions for the stiffness and mass
// matrices
void IsoHexElement::precomputeStiffness()
{
	VEC3F evalPoint;
	VEC3F evalWeights;
	Real weight;

	// Perform multidimensional quadrature to evaluate
	// the integral.
	for (int i = 1; i <= _numGaussPoints; i++)
	{
		// Set the xi coordinate
		evalPoint[0] = GaussPoints::gaussPoint( _numGaussPoints, i );
		evalWeights[0] = GaussPoints::gaussWeight( _numGaussPoints, i );

		for (int j = 1; j <= _numGaussPoints; j++)
		{
			// Set the eta coordinate
			evalPoint[1] = GaussPoints::gaussPoint( _numGaussPoints, j );
			evalWeights[1] = GaussPoints::gaussWeight( _numGaussPoints, j );

			for (int k = 1; k <= _numGaussPoints; k++)
			{
				// Set the zeta coordinate
				evalPoint[2] = GaussPoints::gaussPoint( _numGaussPoints, k );
				evalWeights[2] = GaussPoints::gaussWeight( _numGaussPoints, k );

				weight = evalWeights[0] * evalWeights[1] * evalWeights[2];

				addStiffnessIntegrandEvaluation( evalPoint, weight );
			}
		}
	}
}

void IsoHexElement::precomputeMass()
{
	VEC3F evalPoint;
	VEC3F evalWeights;
	Real weight;

	// Perform multidimensional quadrature to evaluate
	// the integral.
	for (int i = 1; i <= _numGaussPoints; i++)
	{
		// Set the xi coordinate
		evalPoint[0] = GaussPoints::gaussPoint( _numGaussPoints, i );
		evalWeights[0] = GaussPoints::gaussWeight( _numGaussPoints, i );

		for (int j = 1; j <= _numGaussPoints; j++)
		{
			// Set the eta coordinate
			evalPoint[1] = GaussPoints::gaussPoint( _numGaussPoints, j );
			evalWeights[1] = GaussPoints::gaussWeight( _numGaussPoints, j );

			for (int k = 1; k <= _numGaussPoints; k++)
			{
				// Set the zeta coordinate
				evalPoint[2] = GaussPoints::gaussPoint( _numGaussPoints, k );
				evalWeights[2] = GaussPoints::gaussWeight( _numGaussPoints, k );

				weight = evalWeights[0] * evalWeights[1] * evalWeights[2];

				addMassIntegrandEvaluation( evalPoint, weight );
			}
		}
	}
}

// Compute the stiffness integrand at the given xi, eta, zeta coordinates
// and add the evaluation to the total stiffness, scaled by the given
// weight.
// This code is generated by matlab/GenerateStiffnessCode.m
void IsoHexElement::addStiffnessIntegrandEvaluation( VEC3F &p, Real weight )
{
	Real xi = p[0];
	Real eta = p[1];
	Real zeta_ = p[2];

	Real t1 = 1.0-eta;
	Real t2 = 1.0-zeta_;
	Real t3 = t1*t2;
	Real t6 = 1.0+eta;
	Real t7 = t6*t2;
	Real t10 = 1.0+zeta_;
	Real t11 = t1*t10;
	Real t14 = t6*t10;
	Real t17 = -t3*x0+t3*x1+t7*x2-t7*x3-t11*x4+t11*x5+t14*x6-t14*x7;
	Real t18 = 1.0-xi;
	Real t19 = t18*t2;
	Real t21 = 1.0+xi;
	Real t22 = t21*t2;
	Real t26 = t18*t10;
	Real t28 = t21*t10;
	Real t32 = -t19*y0-t22*y1+t22*y2+t19*y3-t26*y4-t28*y5+t28*y6+t26*y7;
	Real t33 = t18*t1;
	Real t35 = t21*t1;
	Real t37 = t21*t6;
	Real t39 = t18*t6;
	Real t45 = -t33*z0-t35*z1-t37*z2-t39*z3+t33*z4+t35*z5+t37*z6+t39*z7;
	Real t55 = -t33*y0-t35*y1-t37*y2-t39*y3+t33*y4+t35*y5+t37*y6+t39*y7;
	Real t64 = -t19*z0-t22*z1+t22*z2+t19*z3-t26*z4-t28*z5+t28*z6+t26*z7;
	Real t66 = t32*t45/64.0-t55*t64/64.0;
	Real t76 = -t19*x0-t22*x1+t22*x2+t19*x3-t26*x4-t28*x5+t28*x6+t26*x7;
	Real t85 = -t3*z0+t3*z1+t7*z2-t7*z3-t11*z4+t11*z5+t14*z6-t14*z7;
	Real t95 = -t3*y0+t3*y1+t7*y2-t7*y3-t11*y4+t11*y5+t14*y6-t14*y7;
	Real t97 = t55*t85/64.0-t95*t45/64.0;
	Real t107 = -t33*x0-t35*x1-t37*x2-t39*x3+t33*x4+t35*x5+t37*x6+t39*x7;
	Real t110 = t95*t64/64.0-t32*t85/64.0;
	Real t112 = t17*t66/8.0+t76*t97/8.0+t107*t110/8.0;
	Real t113 = t3*t66;
	Real t114 = t19*t97;
	Real t115 = t33*t110;
	Real t116 = -t113-t114-t115;
	Real t118 = t112*t112;
	Real t119 = 1/t118;
	Real t121 = t116*t116*t119*E/64.0;
	Real t127 = 1/(1.0+v)/(1.0-2.0*v);
	Real t128 = 1.0-v;
	Real t129 = t127*t128;
	Real t133 = t76*t55/64.0-t107*t32/64.0;
	Real t134 = t3*t133;
	Real t137 = t107*t95/64.0-t17*t55/64.0;
	Real t138 = t19*t137;
	Real t141 = t17*t32/64.0-t76*t95/64.0;
	Real t142 = t33*t141;
	Real t143 = -t134-t138-t142;
	Real t146 = t143*t143*t119*E/64.0;
	Real t147 = 1.0/2.0-v;
	Real t148 = t127*t147;
	Real t149 = t146*t148;
	Real t152 = t64*t107/64.0-t45*t76/64.0;
	Real t153 = t3*t152;
	Real t156 = t45*t17/64.0-t85*t107/64.0;
	Real t157 = t19*t156;
	Real t160 = t85*t76/64.0-t64*t17/64.0;
	Real t161 = t33*t160;
	Real t162 = -t153-t157-t161;
	Real t165 = t162*t162*t119*E/64.0;
	Real t166 = t165*t148;
	Real t170 = t116*t119*E/8.0;
	Real t175 = t162*t119*E/8.0;
	Real t177 = t127*t147*t116/8.0;
	Real t180 = t112*(t170*t127*v*t162/8.0+t175*t177);
	Real t182 = t127*v*t143/8.0;
	Real t185 = t143*t119*E/8.0;
	Real t188 = t112*(t170*t182+t185*t177);
	Real t189 = t22*t97;
	Real t190 = t35*t110;
	Real t191 = t113-t189-t190;
	Real t195 = t22*t137;
	Real t196 = t35*t141;
	Real t197 = t134-t195-t196;
	Real t199 = t127*t147*t197/8.0;
	Real t200 = t185*t199;
	Real t201 = t22*t156;
	Real t202 = t35*t160;
	Real t203 = t153-t201-t202;
	Real t205 = t127*t147*t203/8.0;
	Real t206 = t175*t205;
	Real t208 = t112*(t170*t127*t128*t191/8.0+t200+t206);
	Real t210 = t127*v*t203/8.0;
	Real t213 = t127*t147*t191/8.0;
	Real t216 = t112*(t170*t210+t175*t213);
	Real t218 = t127*v*t197/8.0;
	Real t222 = t112*(t170*t218+t185*t213);
	Real t223 = t7*t66;
	Real t224 = t37*t110;
	Real t225 = t223+t189-t224;
	Real t227 = t127*t128*t225/8.0;
	Real t229 = t7*t133;
	Real t230 = t37*t141;
	Real t231 = t229+t195-t230;
	Real t233 = t127*t147*t231/8.0;
	Real t234 = t185*t233;
	Real t235 = t7*t152;
	Real t236 = t37*t160;
	Real t237 = t235+t201-t236;
	Real t239 = t127*t147*t237/8.0;
	Real t240 = t175*t239;
	Real t242 = t112*(t170*t227+t234+t240);
	Real t244 = t127*v*t237/8.0;
	Real t247 = t127*t147*t225/8.0;
	Real t250 = t112*(t170*t244+t175*t247);
	Real t252 = t127*v*t231/8.0;
	Real t256 = t112*(t170*t252+t185*t247);
	Real t257 = t39*t110;
	Real t258 = -t223+t114-t257;
	Real t260 = t127*t128*t258/8.0;
	Real t262 = t39*t141;
	Real t263 = -t229+t138-t262;
	Real t265 = t127*t147*t263/8.0;
	Real t266 = t185*t265;
	Real t267 = t39*t160;
	Real t268 = -t235+t157-t267;
	Real t270 = t127*t147*t268/8.0;
	Real t271 = t175*t270;
	Real t273 = t112*(t170*t260+t266+t271);
	Real t275 = t127*v*t268/8.0;
	Real t278 = t127*t147*t258/8.0;
	Real t281 = t112*(t170*t275+t175*t278);
	Real t283 = t127*v*t263/8.0;
	Real t287 = t112*(t170*t283+t185*t278);
	Real t288 = t11*t66;
	Real t289 = t26*t97;
	Real t290 = -t288-t289+t115;
	Real t292 = t127*t128*t290/8.0;
	Real t294 = t11*t133;
	Real t295 = t26*t137;
	Real t296 = -t294-t295+t142;
	Real t298 = t127*t147*t296/8.0;
	Real t299 = t185*t298;
	Real t300 = t11*t152;
	Real t301 = t26*t156;
	Real t302 = -t300-t301+t161;
	Real t304 = t127*t147*t302/8.0;
	Real t305 = t175*t304;
	Real t307 = t112*(t170*t292+t299+t305);
	Real t309 = t127*v*t302/8.0;
	Real t312 = t127*t147*t290/8.0;
	Real t315 = t112*(t170*t309+t175*t312);
	Real t317 = t127*v*t296/8.0;
	Real t321 = t112*(t170*t317+t185*t312);
	Real t322 = t28*t97;
	Real t323 = t288-t322+t190;
	Real t325 = t127*t128*t323/8.0;
	Real t327 = t28*t137;
	Real t328 = t294-t327+t196;
	Real t330 = t127*t147*t328/8.0;
	Real t331 = t185*t330;
	Real t332 = t28*t156;
	Real t333 = t300-t332+t202;
	Real t335 = t127*t147*t333/8.0;
	Real t336 = t175*t335;
	Real t338 = t112*(t170*t325+t331+t336);
	Real t340 = t127*v*t333/8.0;
	Real t343 = t127*t147*t323/8.0;
	Real t346 = t112*(t170*t340+t175*t343);
	Real t348 = t127*v*t328/8.0;
	Real t352 = t112*(t170*t348+t185*t343);
	Real t353 = t14*t66;
	Real t354 = t353+t322+t224;
	Real t356 = t127*t128*t354/8.0;
	Real t358 = t14*t133;
	Real t359 = t358+t327+t230;
	Real t361 = t127*t147*t359/8.0;
	Real t362 = t185*t361;
	Real t363 = t14*t152;
	Real t364 = t363+t332+t236;
	Real t366 = t127*t147*t364/8.0;
	Real t367 = t175*t366;
	Real t369 = t112*(t170*t356+t362+t367);
	Real t371 = t127*v*t364/8.0;
	Real t374 = t127*t147*t354/8.0;
	Real t377 = t112*(t170*t371+t175*t374);
	Real t379 = t127*v*t359/8.0;
	Real t383 = t112*(t170*t379+t185*t374);
	Real t384 = -t353+t289+t257;
	Real t386 = t127*t128*t384/8.0;
	Real t388 = -t358+t295+t262;
	Real t390 = t127*t147*t388/8.0;
	Real t391 = t185*t390;
	Real t392 = -t363+t301+t267;
	Real t394 = t127*t147*t392/8.0;
	Real t395 = t175*t394;
	Real t397 = t112*(t170*t386+t391+t395);
	Real t399 = t127*v*t392/8.0;
	Real t402 = t127*t147*t384/8.0;
	Real t405 = t112*(t170*t399+t175*t402);
	Real t407 = t127*v*t388/8.0;
	Real t411 = t112*(t170*t407+t185*t402);
	Real t413 = t121*t148;
	Real t421 = t112*(t175*t182+t185*t127*t147*t162/8.0);
	Real t423 = t127*v*t191/8.0;
	Real t427 = t112*(t175*t423+t170*t205);
	Real t431 = t170*t213;
	Real t433 = t112*(t175*t127*t128*t203/8.0+t200+t431);
	Real t437 = t112*(t175*t218+t185*t205);
	Real t439 = t127*v*t225/8.0;
	Real t443 = t112*(t175*t439+t170*t239);
	Real t445 = t127*t128*t237/8.0;
	Real t447 = t170*t247;
	Real t449 = t112*(t175*t445+t234+t447);
	Real t453 = t112*(t175*t252+t185*t239);
	Real t455 = t127*v*t258/8.0;
	Real t459 = t112*(t175*t455+t170*t270);
	Real t461 = t127*t128*t268/8.0;
	Real t463 = t170*t278;
	Real t465 = t112*(t175*t461+t266+t463);
	Real t469 = t112*(t175*t283+t185*t270);
	Real t471 = t127*v*t290/8.0;
	Real t475 = t112*(t175*t471+t170*t304);
	Real t477 = t127*t128*t302/8.0;
	Real t479 = t170*t312;
	Real t481 = t112*(t175*t477+t299+t479);
	Real t485 = t112*(t175*t317+t185*t304);
	Real t487 = t127*v*t323/8.0;
	Real t491 = t112*(t175*t487+t170*t335);
	Real t493 = t127*t128*t333/8.0;
	Real t495 = t170*t343;
	Real t497 = t112*(t175*t493+t331+t495);
	Real t501 = t112*(t175*t348+t185*t335);
	Real t503 = t127*v*t354/8.0;
	Real t507 = t112*(t175*t503+t170*t366);
	Real t509 = t127*t128*t364/8.0;
	Real t511 = t170*t374;
	Real t513 = t112*(t175*t509+t362+t511);
	Real t517 = t112*(t175*t379+t185*t366);
	Real t519 = t127*v*t384/8.0;
	Real t523 = t112*(t175*t519+t170*t394);
	Real t525 = t127*t128*t392/8.0;
	Real t527 = t170*t402;
	Real t529 = t112*(t175*t525+t391+t527);
	Real t533 = t112*(t175*t407+t185*t394);
	Real t540 = t112*(t185*t423+t170*t199);
	Real t544 = t112*(t185*t210+t175*t199);
	Real t549 = t112*(t185*t127*t128*t197/8.0+t206+t431);
	Real t553 = t112*(t185*t439+t170*t233);
	Real t557 = t112*(t185*t244+t175*t233);
	Real t559 = t127*t128*t231/8.0;
	Real t562 = t112*(t185*t559+t240+t447);
	Real t566 = t112*(t185*t455+t170*t265);
	Real t570 = t112*(t185*t275+t175*t265);
	Real t572 = t127*t128*t263/8.0;
	Real t575 = t112*(t185*t572+t271+t463);
	Real t579 = t112*(t185*t471+t170*t298);
	Real t583 = t112*(t185*t309+t175*t298);
	Real t585 = t127*t128*t296/8.0;
	Real t588 = t112*(t185*t585+t305+t479);
	Real t592 = t112*(t185*t487+t170*t330);
	Real t596 = t112*(t185*t340+t175*t330);
	Real t598 = t127*t128*t328/8.0;
	Real t601 = t112*(t185*t598+t336+t495);
	Real t605 = t112*(t185*t503+t170*t361);
	Real t609 = t112*(t185*t371+t175*t361);
	Real t611 = t127*t128*t359/8.0;
	Real t614 = t112*(t185*t611+t367+t511);
	Real t618 = t112*(t185*t519+t170*t390);
	Real t622 = t112*(t185*t399+t175*t390);
	Real t624 = t127*t128*t388/8.0;
	Real t627 = t112*(t185*t624+t395+t527);
	Real t630 = t191*t191*t119*E/64.0;
	Real t634 = t197*t197*t119*E/64.0;
	Real t635 = t634*t148;
	Real t638 = t203*t203*t119*E/64.0;
	Real t639 = t638*t148;
	Real t643 = t191*t119*E/8.0;
	Real t646 = t203*t119*E/8.0;
	Real t649 = t112*(t643*t210+t646*t213);
	Real t652 = t197*t119*E/8.0;
	Real t655 = t112*(t643*t218+t652*t213);
	Real t657 = t652*t233;
	Real t658 = t646*t239;
	Real t660 = t112*(t643*t227+t657+t658);
	Real t664 = t112*(t643*t244+t646*t247);
	Real t668 = t112*(t643*t252+t652*t247);
	Real t670 = t652*t265;
	Real t671 = t646*t270;
	Real t673 = t112*(t643*t260+t670+t671);
	Real t677 = t112*(t643*t275+t646*t278);
	Real t681 = t112*(t643*t283+t652*t278);
	Real t683 = t652*t298;
	Real t684 = t646*t304;
	Real t686 = t112*(t643*t292+t683+t684);
	Real t690 = t112*(t643*t309+t646*t312);
	Real t694 = t112*(t643*t317+t652*t312);
	Real t696 = t652*t330;
	Real t697 = t646*t335;
	Real t699 = t112*(t643*t325+t696+t697);
	Real t703 = t112*(t643*t340+t646*t343);
	Real t707 = t112*(t643*t348+t652*t343);
	Real t709 = t652*t361;
	Real t710 = t646*t366;
	Real t712 = t112*(t643*t356+t709+t710);
	Real t716 = t112*(t643*t371+t646*t374);
	Real t720 = t112*(t643*t379+t652*t374);
	Real t722 = t652*t390;
	Real t723 = t646*t394;
	Real t725 = t112*(t643*t386+t722+t723);
	Real t729 = t112*(t643*t399+t646*t402);
	Real t733 = t112*(t643*t407+t652*t402);
	Real t735 = t630*t148;
	Real t741 = t112*(t646*t218+t652*t205);
	Real t745 = t112*(t646*t439+t643*t239);
	Real t747 = t643*t247;
	Real t749 = t112*(t646*t445+t657+t747);
	Real t753 = t112*(t646*t252+t652*t239);
	Real t757 = t112*(t646*t455+t643*t270);
	Real t759 = t643*t278;
	Real t761 = t112*(t646*t461+t670+t759);
	Real t765 = t112*(t646*t283+t652*t270);
	Real t769 = t112*(t646*t471+t643*t304);
	Real t771 = t643*t312;
	Real t773 = t112*(t646*t477+t683+t771);
	Real t777 = t112*(t646*t317+t652*t304);
	Real t781 = t112*(t646*t487+t643*t335);
	Real t783 = t643*t343;
	Real t785 = t112*(t646*t493+t696+t783);
	Real t789 = t112*(t646*t348+t652*t335);
	Real t793 = t112*(t646*t503+t643*t366);
	Real t795 = t643*t374;
	Real t797 = t112*(t646*t509+t709+t795);
	Real t801 = t112*(t646*t379+t652*t366);
	Real t805 = t112*(t646*t519+t643*t394);
	Real t807 = t643*t402;
	Real t809 = t112*(t646*t525+t722+t807);
	Real t813 = t112*(t646*t407+t652*t394);
	Real t820 = t112*(t652*t439+t643*t233);
	Real t824 = t112*(t652*t244+t646*t233);
	Real t827 = t112*(t652*t559+t658+t747);
	Real t831 = t112*(t652*t455+t643*t265);
	Real t835 = t112*(t652*t275+t646*t265);
	Real t838 = t112*(t652*t572+t671+t759);
	Real t842 = t112*(t652*t471+t643*t298);
	Real t846 = t112*(t652*t309+t646*t298);
	Real t849 = t112*(t652*t585+t684+t771);
	Real t853 = t112*(t652*t487+t643*t330);
	Real t857 = t112*(t652*t340+t646*t330);
	Real t860 = t112*(t652*t598+t697+t783);
	Real t864 = t112*(t652*t503+t643*t361);
	Real t868 = t112*(t652*t371+t646*t361);
	Real t871 = t112*(t652*t611+t710+t795);
	Real t875 = t112*(t652*t519+t643*t390);
	Real t879 = t112*(t652*t399+t646*t390);
	Real t882 = t112*(t652*t624+t723+t807);
	Real t885 = t225*t225*t119*E/64.0;
	Real t889 = t231*t231*t119*E/64.0;
	Real t890 = t889*t148;
	Real t893 = t237*t237*t119*E/64.0;
	Real t894 = t893*t148;
	Real t898 = t225*t119*E/8.0;
	Real t901 = t237*t119*E/8.0;
	Real t904 = t112*(t898*t244+t901*t247);
	Real t907 = t231*t119*E/8.0;
	Real t910 = t112*(t898*t252+t907*t247);
	Real t912 = t907*t265;
	Real t913 = t901*t270;
	Real t915 = t112*(t898*t260+t912+t913);
	Real t919 = t112*(t898*t275+t901*t278);
	Real t923 = t112*(t898*t283+t907*t278);
	Real t925 = t907*t298;
	Real t926 = t901*t304;
	Real t928 = t112*(t898*t292+t925+t926);
	Real t932 = t112*(t898*t309+t901*t312);
	Real t936 = t112*(t898*t317+t907*t312);
	Real t938 = t907*t330;
	Real t939 = t901*t335;
	Real t941 = t112*(t898*t325+t938+t939);
	Real t945 = t112*(t898*t340+t901*t343);
	Real t949 = t112*(t898*t348+t907*t343);
	Real t951 = t907*t361;
	Real t952 = t901*t366;
	Real t954 = t112*(t898*t356+t951+t952);
	Real t958 = t112*(t898*t371+t901*t374);
	Real t962 = t112*(t898*t379+t907*t374);
	Real t964 = t907*t390;
	Real t965 = t901*t394;
	Real t967 = t112*(t898*t386+t964+t965);
	Real t971 = t112*(t898*t399+t901*t402);
	Real t975 = t112*(t898*t407+t907*t402);
	Real t977 = t885*t148;
	Real t983 = t112*(t901*t252+t907*t239);
	Real t987 = t112*(t901*t455+t898*t270);
	Real t989 = t898*t278;
	Real t991 = t112*(t901*t461+t912+t989);
	Real t995 = t112*(t901*t283+t907*t270);
	Real t999 = t112*(t901*t471+t898*t304);
	Real t1001 = t898*t312;
	Real t1003 = t112*(t901*t477+t925+t1001);
	Real t1007 = t112*(t901*t317+t907*t304);
	Real t1011 = t112*(t901*t487+t898*t335);
	Real t1013 = t898*t343;
	Real t1015 = t112*(t901*t493+t938+t1013);
	Real t1019 = t112*(t901*t348+t907*t335);
	Real t1023 = t112*(t901*t503+t898*t366);
	Real t1025 = t898*t374;
	Real t1027 = t112*(t901*t509+t951+t1025);
	Real t1031 = t112*(t901*t379+t907*t366);
	Real t1035 = t112*(t901*t519+t898*t394);
	Real t1037 = t898*t402;
	Real t1039 = t112*(t901*t525+t964+t1037);
	Real t1043 = t112*(t901*t407+t907*t394);
	Real t1050 = t112*(t907*t455+t898*t265);
	Real t1054 = t112*(t907*t275+t901*t265);
	Real t1057 = t112*(t907*t572+t913+t989);
	Real t1061 = t112*(t907*t471+t898*t298);
	Real t1065 = t112*(t907*t309+t901*t298);
	Real t1068 = t112*(t907*t585+t926+t1001);
	Real t1072 = t112*(t907*t487+t898*t330);
	Real t1076 = t112*(t907*t340+t901*t330);
	Real t1079 = t112*(t907*t598+t939+t1013);
	Real t1083 = t112*(t907*t503+t898*t361);
	Real t1087 = t112*(t907*t371+t901*t361);
	Real t1090 = t112*(t907*t611+t952+t1025);
	Real t1094 = t112*(t907*t519+t898*t390);
	Real t1098 = t112*(t907*t399+t901*t390);
	Real t1101 = t112*(t907*t624+t965+t1037);
	Real t1104 = t258*t258*t119*E/64.0;
	Real t1108 = t263*t263*t119*E/64.0;
	Real t1109 = t1108*t148;
	Real t1112 = t268*t268*t119*E/64.0;
	Real t1113 = t1112*t148;
	Real t1117 = t258*t119*E/8.0;
	Real t1120 = t268*t119*E/8.0;
	Real t1123 = t112*(t1117*t275+t1120*t278);
	Real t1126 = t263*t119*E/8.0;
	Real t1129 = t112*(t1117*t283+t1126*t278);
	Real t1131 = t1126*t298;
	Real t1132 = t1120*t304;
	Real t1134 = t112*(t1117*t292+t1131+t1132);
	Real t1138 = t112*(t1117*t309+t1120*t312);
	Real t1142 = t112*(t1117*t317+t1126*t312);
	Real t1144 = t1126*t330;
	Real t1145 = t1120*t335;
	Real t1147 = t112*(t1117*t325+t1144+t1145);
	Real t1151 = t112*(t1117*t340+t1120*t343);
	Real t1155 = t112*(t1117*t348+t1126*t343);
	Real t1157 = t1126*t361;
	Real t1158 = t1120*t366;
	Real t1160 = t112*(t1117*t356+t1157+t1158);
	Real t1164 = t112*(t1117*t371+t1120*t374);
	Real t1168 = t112*(t1117*t379+t1126*t374);
	Real t1170 = t1126*t390;
	Real t1171 = t1120*t394;
	Real t1173 = t112*(t1117*t386+t1170+t1171);
	Real t1177 = t112*(t1117*t399+t1120*t402);
	Real t1181 = t112*(t1117*t407+t1126*t402);
	Real t1183 = t1104*t148;
	Real t1189 = t112*(t1120*t283+t1126*t270);
	Real t1193 = t112*(t1120*t471+t1117*t304);
	Real t1195 = t1117*t312;
	Real t1197 = t112*(t1120*t477+t1131+t1195);
	Real t1201 = t112*(t1120*t317+t1126*t304);
	Real t1205 = t112*(t1120*t487+t1117*t335);
	Real t1207 = t1117*t343;
	Real t1209 = t112*(t1120*t493+t1144+t1207);
	Real t1213 = t112*(t1120*t348+t1126*t335);
	Real t1217 = t112*(t1120*t503+t1117*t366);
	Real t1219 = t1117*t374;
	Real t1221 = t112*(t1120*t509+t1157+t1219);
	Real t1225 = t112*(t1120*t379+t1126*t366);
	Real t1229 = t112*(t1120*t519+t1117*t394);
	Real t1231 = t1117*t402;
	Real t1233 = t112*(t1120*t525+t1170+t1231);
	Real t1237 = t112*(t1120*t407+t1126*t394);
	Real t1244 = t112*(t1126*t471+t1117*t298);
	Real t1248 = t112*(t1126*t309+t1120*t298);
	Real t1251 = t112*(t1126*t585+t1132+t1195);
	Real t1255 = t112*(t1126*t487+t1117*t330);
	Real t1259 = t112*(t1126*t340+t1120*t330);
	Real t1262 = t112*(t1126*t598+t1145+t1207);
	Real t1266 = t112*(t1126*t503+t1117*t361);
	Real t1270 = t112*(t1126*t371+t1120*t361);
	Real t1273 = t112*(t1126*t611+t1158+t1219);
	Real t1277 = t112*(t1126*t519+t1117*t390);
	Real t1281 = t112*(t1126*t399+t1120*t390);
	Real t1284 = t112*(t1126*t624+t1171+t1231);
	Real t1287 = t290*t290*t119*E/64.0;
	Real t1291 = t296*t296*t119*E/64.0;
	Real t1292 = t1291*t148;
	Real t1295 = t302*t302*t119*E/64.0;
	Real t1296 = t1295*t148;
	Real t1300 = t290*t119*E/8.0;
	Real t1303 = t302*t119*E/8.0;
	Real t1306 = t112*(t1300*t309+t1303*t312);
	Real t1309 = t296*t119*E/8.0;
	Real t1312 = t112*(t1300*t317+t1309*t312);
	Real t1314 = t1309*t330;
	Real t1315 = t1303*t335;
	Real t1317 = t112*(t1300*t325+t1314+t1315);
	Real t1321 = t112*(t1300*t340+t1303*t343);
	Real t1325 = t112*(t1300*t348+t1309*t343);
	Real t1327 = t1309*t361;
	Real t1328 = t1303*t366;
	Real t1330 = t112*(t1300*t356+t1327+t1328);
	Real t1334 = t112*(t1300*t371+t1303*t374);
	Real t1338 = t112*(t1300*t379+t1309*t374);
	Real t1340 = t1309*t390;
	Real t1341 = t1303*t394;
	Real t1343 = t112*(t1300*t386+t1340+t1341);
	Real t1347 = t112*(t1300*t399+t1303*t402);
	Real t1351 = t112*(t1300*t407+t1309*t402);
	Real t1353 = t1287*t148;
	Real t1359 = t112*(t1303*t317+t1309*t304);
	Real t1363 = t112*(t1303*t487+t1300*t335);
	Real t1365 = t1300*t343;
	Real t1367 = t112*(t1303*t493+t1314+t1365);
	Real t1371 = t112*(t1303*t348+t1309*t335);
	Real t1375 = t112*(t1303*t503+t1300*t366);
	Real t1377 = t1300*t374;
	Real t1379 = t112*(t1303*t509+t1327+t1377);
	Real t1383 = t112*(t1303*t379+t1309*t366);
	Real t1387 = t112*(t1303*t519+t1300*t394);
	Real t1389 = t1300*t402;
	Real t1391 = t112*(t1303*t525+t1340+t1389);
	Real t1395 = t112*(t1303*t407+t1309*t394);
	Real t1402 = t112*(t1309*t487+t1300*t330);
	Real t1406 = t112*(t1309*t340+t1303*t330);
	Real t1409 = t112*(t1309*t598+t1315+t1365);
	Real t1413 = t112*(t1309*t503+t1300*t361);
	Real t1417 = t112*(t1309*t371+t1303*t361);
	Real t1420 = t112*(t1309*t611+t1328+t1377);
	Real t1424 = t112*(t1309*t519+t1300*t390);
	Real t1428 = t112*(t1309*t399+t1303*t390);
	Real t1431 = t112*(t1309*t624+t1341+t1389);
	Real t1434 = t323*t323*t119*E/64.0;
	Real t1438 = t328*t328*t119*E/64.0;
	Real t1439 = t1438*t148;
	Real t1442 = t333*t333*t119*E/64.0;
	Real t1443 = t1442*t148;
	Real t1447 = t323*t119*E/8.0;
	Real t1450 = t333*t119*E/8.0;
	Real t1453 = t112*(t1447*t340+t1450*t343);
	Real t1456 = t328*t119*E/8.0;
	Real t1459 = t112*(t1447*t348+t1456*t343);
	Real t1461 = t1456*t361;
	Real t1462 = t1450*t366;
	Real t1464 = t112*(t1447*t356+t1461+t1462);
	Real t1468 = t112*(t1447*t371+t1450*t374);
	Real t1472 = t112*(t1447*t379+t1456*t374);
	Real t1474 = t1456*t390;
	Real t1475 = t1450*t394;
	Real t1477 = t112*(t1447*t386+t1474+t1475);
	Real t1481 = t112*(t1447*t399+t1450*t402);
	Real t1485 = t112*(t1447*t407+t1456*t402);
	Real t1487 = t1434*t148;
	Real t1493 = t112*(t1450*t348+t1456*t335);
	Real t1497 = t112*(t1450*t503+t1447*t366);
	Real t1499 = t1447*t374;
	Real t1501 = t112*(t1450*t509+t1461+t1499);
	Real t1505 = t112*(t1450*t379+t1456*t366);
	Real t1509 = t112*(t1450*t519+t1447*t394);
	Real t1511 = t1447*t402;
	Real t1513 = t112*(t1450*t525+t1474+t1511);
	Real t1517 = t112*(t1450*t407+t1456*t394);
	Real t1524 = t112*(t1456*t503+t1447*t361);
	Real t1528 = t112*(t1456*t371+t1450*t361);
	Real t1531 = t112*(t1456*t611+t1462+t1499);
	Real t1535 = t112*(t1456*t519+t1447*t390);
	Real t1539 = t112*(t1456*t399+t1450*t390);
	Real t1542 = t112*(t1456*t624+t1475+t1511);
	Real t1545 = t354*t354*t119*E/64.0;
	Real t1549 = t359*t359*t119*E/64.0;
	Real t1550 = t1549*t148;
	Real t1553 = t364*t364*t119*E/64.0;
	Real t1554 = t1553*t148;
	Real t1558 = t354*t119*E/8.0;
	Real t1561 = t364*t119*E/8.0;
	Real t1564 = t112*(t1558*t371+t1561*t374);
	Real t1567 = t359*t119*E/8.0;
	Real t1570 = t112*(t1558*t379+t1567*t374);
	Real t1572 = t1567*t390;
	Real t1573 = t1561*t394;
	Real t1575 = t112*(t1558*t386+t1572+t1573);
	Real t1579 = t112*(t1558*t399+t1561*t402);
	Real t1583 = t112*(t1558*t407+t1567*t402);
	Real t1585 = t1545*t148;
	Real t1591 = t112*(t1561*t379+t1567*t366);
	Real t1595 = t112*(t1561*t519+t1558*t394);
	Real t1597 = t1558*t402;
	Real t1599 = t112*(t1561*t525+t1572+t1597);
	Real t1603 = t112*(t1561*t407+t1567*t394);
	Real t1610 = t112*(t1567*t519+t1558*t390);
	Real t1614 = t112*(t1567*t399+t1561*t390);
	Real t1617 = t112*(t1567*t624+t1573+t1597);
	Real t1620 = t384*t384*t119*E/64.0;
	Real t1624 = t388*t388*t119*E/64.0;
	Real t1625 = t1624*t148;
	Real t1628 = t392*t392*t119*E/64.0;
	Real t1629 = t1628*t148;
	Real t1633 = t384*t119*E/8.0;
	Real t1636 = t392*t119*E/8.0;
	Real t1639 = t112*(t1633*t399+t1636*t402);
	Real t1642 = t388*t119*E/8.0;
	Real t1645 = t112*(t1633*t407+t1642*t402);
	Real t1647 = t1620*t148;
	Real t1653 = t112*(t1636*t407+t1642*t394);

	// Add to the local stiffness matrix
	_localStiffness(0, 0) += weight * t112*(t121*t129+t149+t166);
	_localStiffness(0, 1) += weight * t180;
	_localStiffness(0, 2) += weight * t188;
	_localStiffness(0, 3) += weight * t208;
	_localStiffness(0, 4) += weight * t216;
	_localStiffness(0, 5) += weight * t222;
	_localStiffness(0, 6) += weight * t242;
	_localStiffness(0, 7) += weight * t250;
	_localStiffness(0, 8) += weight * t256;
	_localStiffness(0, 9) += weight * t273;
	_localStiffness(0, 10) += weight * t281;
	_localStiffness(0, 11) += weight * t287;
	_localStiffness(0, 12) += weight * t307;
	_localStiffness(0, 13) += weight * t315;
	_localStiffness(0, 14) += weight * t321;
	_localStiffness(0, 15) += weight * t338;
	_localStiffness(0, 16) += weight * t346;
	_localStiffness(0, 17) += weight * t352;
	_localStiffness(0, 18) += weight * t369;
	_localStiffness(0, 19) += weight * t377;
	_localStiffness(0, 20) += weight * t383;
	_localStiffness(0, 21) += weight * t397;
	_localStiffness(0, 22) += weight * t405;
	_localStiffness(0, 23) += weight * t411;
	_localStiffness(1, 0) += weight * t180;
	_localStiffness(1, 1) += weight * t112*(t165*t129+t149+t413);
	_localStiffness(1, 2) += weight * t421;
	_localStiffness(1, 3) += weight * t427;
	_localStiffness(1, 4) += weight * t433;
	_localStiffness(1, 5) += weight * t437;
	_localStiffness(1, 6) += weight * t443;
	_localStiffness(1, 7) += weight * t449;
	_localStiffness(1, 8) += weight * t453;
	_localStiffness(1, 9) += weight * t459;
	_localStiffness(1, 10) += weight * t465;
	_localStiffness(1, 11) += weight * t469;
	_localStiffness(1, 12) += weight * t475;
	_localStiffness(1, 13) += weight * t481;
	_localStiffness(1, 14) += weight * t485;
	_localStiffness(1, 15) += weight * t491;
	_localStiffness(1, 16) += weight * t497;
	_localStiffness(1, 17) += weight * t501;
	_localStiffness(1, 18) += weight * t507;
	_localStiffness(1, 19) += weight * t513;
	_localStiffness(1, 20) += weight * t517;
	_localStiffness(1, 21) += weight * t523;
	_localStiffness(1, 22) += weight * t529;
	_localStiffness(1, 23) += weight * t533;
	_localStiffness(2, 0) += weight * t188;
	_localStiffness(2, 1) += weight * t421;
	_localStiffness(2, 2) += weight * t112*(t146*t129+t166+t413);
	_localStiffness(2, 3) += weight * t540;
	_localStiffness(2, 4) += weight * t544;
	_localStiffness(2, 5) += weight * t549;
	_localStiffness(2, 6) += weight * t553;
	_localStiffness(2, 7) += weight * t557;
	_localStiffness(2, 8) += weight * t562;
	_localStiffness(2, 9) += weight * t566;
	_localStiffness(2, 10) += weight * t570;
	_localStiffness(2, 11) += weight * t575;
	_localStiffness(2, 12) += weight * t579;
	_localStiffness(2, 13) += weight * t583;
	_localStiffness(2, 14) += weight * t588;
	_localStiffness(2, 15) += weight * t592;
	_localStiffness(2, 16) += weight * t596;
	_localStiffness(2, 17) += weight * t601;
	_localStiffness(2, 18) += weight * t605;
	_localStiffness(2, 19) += weight * t609;
	_localStiffness(2, 20) += weight * t614;
	_localStiffness(2, 21) += weight * t618;
	_localStiffness(2, 22) += weight * t622;
	_localStiffness(2, 23) += weight * t627;
	_localStiffness(3, 0) += weight * t208;
	_localStiffness(3, 1) += weight * t427;
	_localStiffness(3, 2) += weight * t540;
	_localStiffness(3, 3) += weight * t112*(t630*t129+t635+t639);
	_localStiffness(3, 4) += weight * t649;
	_localStiffness(3, 5) += weight * t655;
	_localStiffness(3, 6) += weight * t660;
	_localStiffness(3, 7) += weight * t664;
	_localStiffness(3, 8) += weight * t668;
	_localStiffness(3, 9) += weight * t673;
	_localStiffness(3, 10) += weight * t677;
	_localStiffness(3, 11) += weight * t681;
	_localStiffness(3, 12) += weight * t686;
	_localStiffness(3, 13) += weight * t690;
	_localStiffness(3, 14) += weight * t694;
	_localStiffness(3, 15) += weight * t699;
	_localStiffness(3, 16) += weight * t703;
	_localStiffness(3, 17) += weight * t707;
	_localStiffness(3, 18) += weight * t712;
	_localStiffness(3, 19) += weight * t716;
	_localStiffness(3, 20) += weight * t720;
	_localStiffness(3, 21) += weight * t725;
	_localStiffness(3, 22) += weight * t729;
	_localStiffness(3, 23) += weight * t733;
	_localStiffness(4, 0) += weight * t216;
	_localStiffness(4, 1) += weight * t433;
	_localStiffness(4, 2) += weight * t544;
	_localStiffness(4, 3) += weight * t649;
	_localStiffness(4, 4) += weight * t112*(t638*t129+t635+t735);
	_localStiffness(4, 5) += weight * t741;
	_localStiffness(4, 6) += weight * t745;
	_localStiffness(4, 7) += weight * t749;
	_localStiffness(4, 8) += weight * t753;
	_localStiffness(4, 9) += weight * t757;
	_localStiffness(4, 10) += weight * t761;
	_localStiffness(4, 11) += weight * t765;
	_localStiffness(4, 12) += weight * t769;
	_localStiffness(4, 13) += weight * t773;
	_localStiffness(4, 14) += weight * t777;
	_localStiffness(4, 15) += weight * t781;
	_localStiffness(4, 16) += weight * t785;
	_localStiffness(4, 17) += weight * t789;
	_localStiffness(4, 18) += weight * t793;
	_localStiffness(4, 19) += weight * t797;
	_localStiffness(4, 20) += weight * t801;
	_localStiffness(4, 21) += weight * t805;
	_localStiffness(4, 22) += weight * t809;
	_localStiffness(4, 23) += weight * t813;
	_localStiffness(5, 0) += weight * t222;
	_localStiffness(5, 1) += weight * t437;
	_localStiffness(5, 2) += weight * t549;
	_localStiffness(5, 3) += weight * t655;
	_localStiffness(5, 4) += weight * t741;
	_localStiffness(5, 5) += weight * t112*(t634*t129+t639+t735);
	_localStiffness(5, 6) += weight * t820;
	_localStiffness(5, 7) += weight * t824;
	_localStiffness(5, 8) += weight * t827;
	_localStiffness(5, 9) += weight * t831;
	_localStiffness(5, 10) += weight * t835;
	_localStiffness(5, 11) += weight * t838;
	_localStiffness(5, 12) += weight * t842;
	_localStiffness(5, 13) += weight * t846;
	_localStiffness(5, 14) += weight * t849;
	_localStiffness(5, 15) += weight * t853;
	_localStiffness(5, 16) += weight * t857;
	_localStiffness(5, 17) += weight * t860;
	_localStiffness(5, 18) += weight * t864;
	_localStiffness(5, 19) += weight * t868;
	_localStiffness(5, 20) += weight * t871;
	_localStiffness(5, 21) += weight * t875;
	_localStiffness(5, 22) += weight * t879;
	_localStiffness(5, 23) += weight * t882;
	_localStiffness(6, 0) += weight * t242;
	_localStiffness(6, 1) += weight * t443;
	_localStiffness(6, 2) += weight * t553;
	_localStiffness(6, 3) += weight * t660;
	_localStiffness(6, 4) += weight * t745;
	_localStiffness(6, 5) += weight * t820;
	_localStiffness(6, 6) += weight * t112*(t885*t129+t890+t894);
	_localStiffness(6, 7) += weight * t904;
	_localStiffness(6, 8) += weight * t910;
	_localStiffness(6, 9) += weight * t915;
	_localStiffness(6, 10) += weight * t919;
	_localStiffness(6, 11) += weight * t923;
	_localStiffness(6, 12) += weight * t928;
	_localStiffness(6, 13) += weight * t932;
	_localStiffness(6, 14) += weight * t936;
	_localStiffness(6, 15) += weight * t941;
	_localStiffness(6, 16) += weight * t945;
	_localStiffness(6, 17) += weight * t949;
	_localStiffness(6, 18) += weight * t954;
	_localStiffness(6, 19) += weight * t958;
	_localStiffness(6, 20) += weight * t962;
	_localStiffness(6, 21) += weight * t967;
	_localStiffness(6, 22) += weight * t971;
	_localStiffness(6, 23) += weight * t975;
	_localStiffness(7, 0) += weight * t250;
	_localStiffness(7, 1) += weight * t449;
	_localStiffness(7, 2) += weight * t557;
	_localStiffness(7, 3) += weight * t664;
	_localStiffness(7, 4) += weight * t749;
	_localStiffness(7, 5) += weight * t824;
	_localStiffness(7, 6) += weight * t904;
	_localStiffness(7, 7) += weight * t112*(t893*t129+t890+t977);
	_localStiffness(7, 8) += weight * t983;
	_localStiffness(7, 9) += weight * t987;
	_localStiffness(7, 10) += weight * t991;
	_localStiffness(7, 11) += weight * t995;
	_localStiffness(7, 12) += weight * t999;
	_localStiffness(7, 13) += weight * t1003;
	_localStiffness(7, 14) += weight * t1007;
	_localStiffness(7, 15) += weight * t1011;
	_localStiffness(7, 16) += weight * t1015;
	_localStiffness(7, 17) += weight * t1019;
	_localStiffness(7, 18) += weight * t1023;
	_localStiffness(7, 19) += weight * t1027;
	_localStiffness(7, 20) += weight * t1031;
	_localStiffness(7, 21) += weight * t1035;
	_localStiffness(7, 22) += weight * t1039;
	_localStiffness(7, 23) += weight * t1043;
	_localStiffness(8, 0) += weight * t256;
	_localStiffness(8, 1) += weight * t453;
	_localStiffness(8, 2) += weight * t562;
	_localStiffness(8, 3) += weight * t668;
	_localStiffness(8, 4) += weight * t753;
	_localStiffness(8, 5) += weight * t827;
	_localStiffness(8, 6) += weight * t910;
	_localStiffness(8, 7) += weight * t983;
	_localStiffness(8, 8) += weight * t112*(t889*t129+t894+t977);
	_localStiffness(8, 9) += weight * t1050;
	_localStiffness(8, 10) += weight * t1054;
	_localStiffness(8, 11) += weight * t1057;
	_localStiffness(8, 12) += weight * t1061;
	_localStiffness(8, 13) += weight * t1065;
	_localStiffness(8, 14) += weight * t1068;
	_localStiffness(8, 15) += weight * t1072;
	_localStiffness(8, 16) += weight * t1076;
	_localStiffness(8, 17) += weight * t1079;
	_localStiffness(8, 18) += weight * t1083;
	_localStiffness(8, 19) += weight * t1087;
	_localStiffness(8, 20) += weight * t1090;
	_localStiffness(8, 21) += weight * t1094;
	_localStiffness(8, 22) += weight * t1098;
	_localStiffness(8, 23) += weight * t1101;
	_localStiffness(9, 0) += weight * t273;
	_localStiffness(9, 1) += weight * t459;
	_localStiffness(9, 2) += weight * t566;
	_localStiffness(9, 3) += weight * t673;
	_localStiffness(9, 4) += weight * t757;
	_localStiffness(9, 5) += weight * t831;
	_localStiffness(9, 6) += weight * t915;
	_localStiffness(9, 7) += weight * t987;
	_localStiffness(9, 8) += weight * t1050;
	_localStiffness(9, 9) += weight * t112*(t1104*t129+t1109+t1113);
	_localStiffness(9, 10) += weight * t1123;
	_localStiffness(9, 11) += weight * t1129;
	_localStiffness(9, 12) += weight * t1134;
	_localStiffness(9, 13) += weight * t1138;
	_localStiffness(9, 14) += weight * t1142;
	_localStiffness(9, 15) += weight * t1147;
	_localStiffness(9, 16) += weight * t1151;
	_localStiffness(9, 17) += weight * t1155;
	_localStiffness(9, 18) += weight * t1160;
	_localStiffness(9, 19) += weight * t1164;
	_localStiffness(9, 20) += weight * t1168;
	_localStiffness(9, 21) += weight * t1173;
	_localStiffness(9, 22) += weight * t1177;
	_localStiffness(9, 23) += weight * t1181;
	_localStiffness(10, 0) += weight * t281;
	_localStiffness(10, 1) += weight * t465;
	_localStiffness(10, 2) += weight * t570;
	_localStiffness(10, 3) += weight * t677;
	_localStiffness(10, 4) += weight * t761;
	_localStiffness(10, 5) += weight * t835;
	_localStiffness(10, 6) += weight * t919;
	_localStiffness(10, 7) += weight * t991;
	_localStiffness(10, 8) += weight * t1054;
	_localStiffness(10, 9) += weight * t1123;
	_localStiffness(10, 10) += weight * t112*(t1112*t129+t1109+t1183);
	_localStiffness(10, 11) += weight * t1189;
	_localStiffness(10, 12) += weight * t1193;
	_localStiffness(10, 13) += weight * t1197;
	_localStiffness(10, 14) += weight * t1201;
	_localStiffness(10, 15) += weight * t1205;
	_localStiffness(10, 16) += weight * t1209;
	_localStiffness(10, 17) += weight * t1213;
	_localStiffness(10, 18) += weight * t1217;
	_localStiffness(10, 19) += weight * t1221;
	_localStiffness(10, 20) += weight * t1225;
	_localStiffness(10, 21) += weight * t1229;
	_localStiffness(10, 22) += weight * t1233;
	_localStiffness(10, 23) += weight * t1237;
	_localStiffness(11, 0) += weight * t287;
	_localStiffness(11, 1) += weight * t469;
	_localStiffness(11, 2) += weight * t575;
	_localStiffness(11, 3) += weight * t681;
	_localStiffness(11, 4) += weight * t765;
	_localStiffness(11, 5) += weight * t838;
	_localStiffness(11, 6) += weight * t923;
	_localStiffness(11, 7) += weight * t995;
	_localStiffness(11, 8) += weight * t1057;
	_localStiffness(11, 9) += weight * t1129;
	_localStiffness(11, 10) += weight * t1189;
	_localStiffness(11, 11) += weight * t112*(t1108*t129+t1113+t1183);
	_localStiffness(11, 12) += weight * t1244;
	_localStiffness(11, 13) += weight * t1248;
	_localStiffness(11, 14) += weight * t1251;
	_localStiffness(11, 15) += weight * t1255;
	_localStiffness(11, 16) += weight * t1259;
	_localStiffness(11, 17) += weight * t1262;
	_localStiffness(11, 18) += weight * t1266;
	_localStiffness(11, 19) += weight * t1270;
	_localStiffness(11, 20) += weight * t1273;
	_localStiffness(11, 21) += weight * t1277;
	_localStiffness(11, 22) += weight * t1281;
	_localStiffness(11, 23) += weight * t1284;
	_localStiffness(12, 0) += weight * t307;
	_localStiffness(12, 1) += weight * t475;
	_localStiffness(12, 2) += weight * t579;
	_localStiffness(12, 3) += weight * t686;
	_localStiffness(12, 4) += weight * t769;
	_localStiffness(12, 5) += weight * t842;
	_localStiffness(12, 6) += weight * t928;
	_localStiffness(12, 7) += weight * t999;
	_localStiffness(12, 8) += weight * t1061;
	_localStiffness(12, 9) += weight * t1134;
	_localStiffness(12, 10) += weight * t1193;
	_localStiffness(12, 11) += weight * t1244;
	_localStiffness(12, 12) += weight * t112*(t1287*t129+t1292+t1296);
	_localStiffness(12, 13) += weight * t1306;
	_localStiffness(12, 14) += weight * t1312;
	_localStiffness(12, 15) += weight * t1317;
	_localStiffness(12, 16) += weight * t1321;
	_localStiffness(12, 17) += weight * t1325;
	_localStiffness(12, 18) += weight * t1330;
	_localStiffness(12, 19) += weight * t1334;
	_localStiffness(12, 20) += weight * t1338;
	_localStiffness(12, 21) += weight * t1343;
	_localStiffness(12, 22) += weight * t1347;
	_localStiffness(12, 23) += weight * t1351;
	_localStiffness(13, 0) += weight * t315;
	_localStiffness(13, 1) += weight * t481;
	_localStiffness(13, 2) += weight * t583;
	_localStiffness(13, 3) += weight * t690;
	_localStiffness(13, 4) += weight * t773;
	_localStiffness(13, 5) += weight * t846;
	_localStiffness(13, 6) += weight * t932;
	_localStiffness(13, 7) += weight * t1003;
	_localStiffness(13, 8) += weight * t1065;
	_localStiffness(13, 9) += weight * t1138;
	_localStiffness(13, 10) += weight * t1197;
	_localStiffness(13, 11) += weight * t1248;
	_localStiffness(13, 12) += weight * t1306;
	_localStiffness(13, 13) += weight * t112*(t1295*t129+t1292+t1353);
	_localStiffness(13, 14) += weight * t1359;
	_localStiffness(13, 15) += weight * t1363;
	_localStiffness(13, 16) += weight * t1367;
	_localStiffness(13, 17) += weight * t1371;
	_localStiffness(13, 18) += weight * t1375;
	_localStiffness(13, 19) += weight * t1379;
	_localStiffness(13, 20) += weight * t1383;
	_localStiffness(13, 21) += weight * t1387;
	_localStiffness(13, 22) += weight * t1391;
	_localStiffness(13, 23) += weight * t1395;
	_localStiffness(14, 0) += weight * t321;
	_localStiffness(14, 1) += weight * t485;
	_localStiffness(14, 2) += weight * t588;
	_localStiffness(14, 3) += weight * t694;
	_localStiffness(14, 4) += weight * t777;
	_localStiffness(14, 5) += weight * t849;
	_localStiffness(14, 6) += weight * t936;
	_localStiffness(14, 7) += weight * t1007;
	_localStiffness(14, 8) += weight * t1068;
	_localStiffness(14, 9) += weight * t1142;
	_localStiffness(14, 10) += weight * t1201;
	_localStiffness(14, 11) += weight * t1251;
	_localStiffness(14, 12) += weight * t1312;
	_localStiffness(14, 13) += weight * t1359;
	_localStiffness(14, 14) += weight * t112*(t1291*t129+t1296+t1353);
	_localStiffness(14, 15) += weight * t1402;
	_localStiffness(14, 16) += weight * t1406;
	_localStiffness(14, 17) += weight * t1409;
	_localStiffness(14, 18) += weight * t1413;
	_localStiffness(14, 19) += weight * t1417;
	_localStiffness(14, 20) += weight * t1420;
	_localStiffness(14, 21) += weight * t1424;
	_localStiffness(14, 22) += weight * t1428;
	_localStiffness(14, 23) += weight * t1431;
	_localStiffness(15, 0) += weight * t338;
	_localStiffness(15, 1) += weight * t491;
	_localStiffness(15, 2) += weight * t592;
	_localStiffness(15, 3) += weight * t699;
	_localStiffness(15, 4) += weight * t781;
	_localStiffness(15, 5) += weight * t853;
	_localStiffness(15, 6) += weight * t941;
	_localStiffness(15, 7) += weight * t1011;
	_localStiffness(15, 8) += weight * t1072;
	_localStiffness(15, 9) += weight * t1147;
	_localStiffness(15, 10) += weight * t1205;
	_localStiffness(15, 11) += weight * t1255;
	_localStiffness(15, 12) += weight * t1317;
	_localStiffness(15, 13) += weight * t1363;
	_localStiffness(15, 14) += weight * t1402;
	_localStiffness(15, 15) += weight * t112*(t1434*t129+t1439+t1443);
	_localStiffness(15, 16) += weight * t1453;
	_localStiffness(15, 17) += weight * t1459;
	_localStiffness(15, 18) += weight * t1464;
	_localStiffness(15, 19) += weight * t1468;
	_localStiffness(15, 20) += weight * t1472;
	_localStiffness(15, 21) += weight * t1477;
	_localStiffness(15, 22) += weight * t1481;
	_localStiffness(15, 23) += weight * t1485;
	_localStiffness(16, 0) += weight * t346;
	_localStiffness(16, 1) += weight * t497;
	_localStiffness(16, 2) += weight * t596;
	_localStiffness(16, 3) += weight * t703;
	_localStiffness(16, 4) += weight * t785;
	_localStiffness(16, 5) += weight * t857;
	_localStiffness(16, 6) += weight * t945;
	_localStiffness(16, 7) += weight * t1015;
	_localStiffness(16, 8) += weight * t1076;
	_localStiffness(16, 9) += weight * t1151;
	_localStiffness(16, 10) += weight * t1209;
	_localStiffness(16, 11) += weight * t1259;
	_localStiffness(16, 12) += weight * t1321;
	_localStiffness(16, 13) += weight * t1367;
	_localStiffness(16, 14) += weight * t1406;
	_localStiffness(16, 15) += weight * t1453;
	_localStiffness(16, 16) += weight * t112*(t1442*t129+t1439+t1487);
	_localStiffness(16, 17) += weight * t1493;
	_localStiffness(16, 18) += weight * t1497;
	_localStiffness(16, 19) += weight * t1501;
	_localStiffness(16, 20) += weight * t1505;
	_localStiffness(16, 21) += weight * t1509;
	_localStiffness(16, 22) += weight * t1513;
	_localStiffness(16, 23) += weight * t1517;
	_localStiffness(17, 0) += weight * t352;
	_localStiffness(17, 1) += weight * t501;
	_localStiffness(17, 2) += weight * t601;
	_localStiffness(17, 3) += weight * t707;
	_localStiffness(17, 4) += weight * t789;
	_localStiffness(17, 5) += weight * t860;
	_localStiffness(17, 6) += weight * t949;
	_localStiffness(17, 7) += weight * t1019;
	_localStiffness(17, 8) += weight * t1079;
	_localStiffness(17, 9) += weight * t1155;
	_localStiffness(17, 10) += weight * t1213;
	_localStiffness(17, 11) += weight * t1262;
	_localStiffness(17, 12) += weight * t1325;
	_localStiffness(17, 13) += weight * t1371;
	_localStiffness(17, 14) += weight * t1409;
	_localStiffness(17, 15) += weight * t1459;
	_localStiffness(17, 16) += weight * t1493;
	_localStiffness(17, 17) += weight * t112*(t1438*t129+t1443+t1487);
	_localStiffness(17, 18) += weight * t1524;
	_localStiffness(17, 19) += weight * t1528;
	_localStiffness(17, 20) += weight * t1531;
	_localStiffness(17, 21) += weight * t1535;
	_localStiffness(17, 22) += weight * t1539;
	_localStiffness(17, 23) += weight * t1542;
	_localStiffness(18, 0) += weight * t369;
	_localStiffness(18, 1) += weight * t507;
	_localStiffness(18, 2) += weight * t605;
	_localStiffness(18, 3) += weight * t712;
	_localStiffness(18, 4) += weight * t793;
	_localStiffness(18, 5) += weight * t864;
	_localStiffness(18, 6) += weight * t954;
	_localStiffness(18, 7) += weight * t1023;
	_localStiffness(18, 8) += weight * t1083;
	_localStiffness(18, 9) += weight * t1160;
	_localStiffness(18, 10) += weight * t1217;
	_localStiffness(18, 11) += weight * t1266;
	_localStiffness(18, 12) += weight * t1330;
	_localStiffness(18, 13) += weight * t1375;
	_localStiffness(18, 14) += weight * t1413;
	_localStiffness(18, 15) += weight * t1464;
	_localStiffness(18, 16) += weight * t1497;
	_localStiffness(18, 17) += weight * t1524;
	_localStiffness(18, 18) += weight * t112*(t1545*t129+t1550+t1554);
	_localStiffness(18, 19) += weight * t1564;
	_localStiffness(18, 20) += weight * t1570;
	_localStiffness(18, 21) += weight * t1575;
	_localStiffness(18, 22) += weight * t1579;
	_localStiffness(18, 23) += weight * t1583;
	_localStiffness(19, 0) += weight * t377;
	_localStiffness(19, 1) += weight * t513;
	_localStiffness(19, 2) += weight * t609;
	_localStiffness(19, 3) += weight * t716;
	_localStiffness(19, 4) += weight * t797;
	_localStiffness(19, 5) += weight * t868;
	_localStiffness(19, 6) += weight * t958;
	_localStiffness(19, 7) += weight * t1027;
	_localStiffness(19, 8) += weight * t1087;
	_localStiffness(19, 9) += weight * t1164;
	_localStiffness(19, 10) += weight * t1221;
	_localStiffness(19, 11) += weight * t1270;
	_localStiffness(19, 12) += weight * t1334;
	_localStiffness(19, 13) += weight * t1379;
	_localStiffness(19, 14) += weight * t1417;
	_localStiffness(19, 15) += weight * t1468;
	_localStiffness(19, 16) += weight * t1501;
	_localStiffness(19, 17) += weight * t1528;
	_localStiffness(19, 18) += weight * t1564;
	_localStiffness(19, 19) += weight * t112*(t1553*t129+t1550+t1585);
	_localStiffness(19, 20) += weight * t1591;
	_localStiffness(19, 21) += weight * t1595;
	_localStiffness(19, 22) += weight * t1599;
	_localStiffness(19, 23) += weight * t1603;
	_localStiffness(20, 0) += weight * t383;
	_localStiffness(20, 1) += weight * t517;
	_localStiffness(20, 2) += weight * t614;
	_localStiffness(20, 3) += weight * t720;
	_localStiffness(20, 4) += weight * t801;
	_localStiffness(20, 5) += weight * t871;
	_localStiffness(20, 6) += weight * t962;
	_localStiffness(20, 7) += weight * t1031;
	_localStiffness(20, 8) += weight * t1090;
	_localStiffness(20, 9) += weight * t1168;
	_localStiffness(20, 10) += weight * t1225;
	_localStiffness(20, 11) += weight * t1273;
	_localStiffness(20, 12) += weight * t1338;
	_localStiffness(20, 13) += weight * t1383;
	_localStiffness(20, 14) += weight * t1420;
	_localStiffness(20, 15) += weight * t1472;
	_localStiffness(20, 16) += weight * t1505;
	_localStiffness(20, 17) += weight * t1531;
	_localStiffness(20, 18) += weight * t1570;
	_localStiffness(20, 19) += weight * t1591;
	_localStiffness(20, 20) += weight * t112*(t1549*t129+t1554+t1585);
	_localStiffness(20, 21) += weight * t1610;
	_localStiffness(20, 22) += weight * t1614;
	_localStiffness(20, 23) += weight * t1617;
	_localStiffness(21, 0) += weight * t397;
	_localStiffness(21, 1) += weight * t523;
	_localStiffness(21, 2) += weight * t618;
	_localStiffness(21, 3) += weight * t725;
	_localStiffness(21, 4) += weight * t805;
	_localStiffness(21, 5) += weight * t875;
	_localStiffness(21, 6) += weight * t967;
	_localStiffness(21, 7) += weight * t1035;
	_localStiffness(21, 8) += weight * t1094;
	_localStiffness(21, 9) += weight * t1173;
	_localStiffness(21, 10) += weight * t1229;
	_localStiffness(21, 11) += weight * t1277;
	_localStiffness(21, 12) += weight * t1343;
	_localStiffness(21, 13) += weight * t1387;
	_localStiffness(21, 14) += weight * t1424;
	_localStiffness(21, 15) += weight * t1477;
	_localStiffness(21, 16) += weight * t1509;
	_localStiffness(21, 17) += weight * t1535;
	_localStiffness(21, 18) += weight * t1575;
	_localStiffness(21, 19) += weight * t1595;
	_localStiffness(21, 20) += weight * t1610;
	_localStiffness(21, 21) += weight * t112*(t1620*t129+t1625+t1629);
	_localStiffness(21, 22) += weight * t1639;
	_localStiffness(21, 23) += weight * t1645;
	_localStiffness(22, 0) += weight * t405;
	_localStiffness(22, 1) += weight * t529;
	_localStiffness(22, 2) += weight * t622;
	_localStiffness(22, 3) += weight * t729;
	_localStiffness(22, 4) += weight * t809;
	_localStiffness(22, 5) += weight * t879;
	_localStiffness(22, 6) += weight * t971;
	_localStiffness(22, 7) += weight * t1039;
	_localStiffness(22, 8) += weight * t1098;
	_localStiffness(22, 9) += weight * t1177;
	_localStiffness(22, 10) += weight * t1233;
	_localStiffness(22, 11) += weight * t1281;
	_localStiffness(22, 12) += weight * t1347;
	_localStiffness(22, 13) += weight * t1391;
	_localStiffness(22, 14) += weight * t1428;
	_localStiffness(22, 15) += weight * t1481;
	_localStiffness(22, 16) += weight * t1513;
	_localStiffness(22, 17) += weight * t1539;
	_localStiffness(22, 18) += weight * t1579;
	_localStiffness(22, 19) += weight * t1599;
	_localStiffness(22, 20) += weight * t1614;
	_localStiffness(22, 21) += weight * t1639;
	_localStiffness(22, 22) += weight * t112*(t1628*t129+t1625+t1647);
	_localStiffness(22, 23) += weight * t1653;
	_localStiffness(23, 0) += weight * t411;
	_localStiffness(23, 1) += weight * t533;
	_localStiffness(23, 2) += weight * t627;
	_localStiffness(23, 3) += weight * t733;
	_localStiffness(23, 4) += weight * t813;
	_localStiffness(23, 5) += weight * t882;
	_localStiffness(23, 6) += weight * t975;
	_localStiffness(23, 7) += weight * t1043;
	_localStiffness(23, 8) += weight * t1101;
	_localStiffness(23, 9) += weight * t1181;
	_localStiffness(23, 10) += weight * t1237;
	_localStiffness(23, 11) += weight * t1284;
	_localStiffness(23, 12) += weight * t1351;
	_localStiffness(23, 13) += weight * t1395;
	_localStiffness(23, 14) += weight * t1431;
	_localStiffness(23, 15) += weight * t1485;
	_localStiffness(23, 16) += weight * t1517;
	_localStiffness(23, 17) += weight * t1542;
	_localStiffness(23, 18) += weight * t1583;
	_localStiffness(23, 19) += weight * t1603;
	_localStiffness(23, 20) += weight * t1617;
	_localStiffness(23, 21) += weight * t1645;
	_localStiffness(23, 22) += weight * t1653;
	_localStiffness(23, 23) += weight * t112*(t1624*t129+t1629+t1647);
}

// Compute the mass integrand at the given xi, eta, zeta coordinates
// and add the evaluation to the total mass, scaled by the given
// weight.
// This code is generated by matlab/GenerateMassCode.m
void IsoHexElement::addMassIntegrandEvaluation( VEC3F &p, Real weight )
{
	Real xi = p[0];
	Real eta = p[1];
	Real zeta_ = p[2];

	Real t1 = 1.0-eta;
	Real t2 = 1.0-zeta_;
	Real t3 = t1*t2;
	Real t6 = 1.0+eta;
	Real t7 = t6*t2;
	Real t10 = 1.0+zeta_;
	Real t11 = t1*t10;
	Real t14 = t6*t10;
	Real t18 = 1.0-xi;
	Real t19 = t18*t2;
	Real t21 = 1.0+xi;
	Real t22 = t21*t2;
	Real t26 = t18*t10;
	Real t28 = t21*t10;
	Real t32 = -t19*y0-t22*y1+t22*y2+t19*y3-t26*y4-t28*y5+t28*y6+t26*y7;
	Real t33 = t18*t1;
	Real t35 = t21*t1;
	Real t37 = t21*t6;
	Real t39 = t18*t6;
	Real t45 = -t33*z0-t35*z1-t37*z2-t39*z3+t33*z4+t35*z5+t37*z6+t39*z7;
	Real t55 = -t33*y0-t35*y1-t37*y2-t39*y3+t33*y4+t35*y5+t37*y6+t39*y7;
	Real t64 = -t19*z0-t22*z1+t22*z2+t19*z3-t26*z4-t28*z5+t28*z6+t26*z7;
	Real t85 = -t3*z0+t3*z1+t7*z2-t7*z3-t11*z4+t11*z5+t14*z6-t14*z7;
	Real t95 = -t3*y0+t3*y1+t7*y2-t7*y3-t11*y4+t11*y5+t14*y6-t14*y7;
	Real t113 = ((-t3*x0+t3*x1+t7*x2-t7*x3-t11*x4+t11*x5+t14*x6-t14*x7)*(t32*t45/64.0-t55*t64/64.0)/8.0+(-t19*x0-t22*x1+t22*x2+t19*x3-t26*x4-t28*x5+t28*x6+t26*x7)*(t55*t85/64.0-t95*t45/64.0)/8.0+(-t33*x0-t35*x1-t37*x2-t39*x3+t33*x4+t35*x5+t37*x6+t39*x7)*(t95*t64/64.0-t32*t85/64.0)/8.0)*density;
	Real t114 = t18*t18;
	Real t115 = t1*t1;
	Real t116 = t114*t115;
	Real t117 = t2*t2;
	Real t121 = t113*t18;
	Real t125 = t121*t115*t117*t21/64.0;
	Real t126 = t1*t117;
	Real t129 = t121*t126*t37/64.0;
	Real t130 = t113*t114;
	Real t131 = t126*t6;
	Real t133 = t130*t131/64.0;
	Real t134 = t115*t2;
	Real t135 = t134*t10;
	Real t137 = t130*t135/64.0;
	Real t140 = t121*t134*t28/64.0;
	Real t144 = t113*t33*t22*t14/64.0;
	Real t145 = t3*t14;
	Real t147 = t130*t145/64.0;
	Real t148 = t21*t21;
	Real t149 = t148*t115;
	Real t153 = t113*t148;
	Real t155 = t153*t131/64.0;
	Real t157 = t153*t135/64.0;
	Real t159 = t153*t145/64.0;
	Real t160 = t6*t6;
	Real t161 = t148*t160;
	Real t165 = t113*t21;
	Real t169 = t165*t160*t117*t18/64.0;
	Real t170 = t160*t2;
	Real t171 = t170*t10;
	Real t173 = t153*t171/64.0;
	Real t176 = t165*t170*t26/64.0;
	Real t177 = t114*t160;
	Real t182 = t130*t171/64.0;
	Real t183 = t10*t10;
	Real t190 = t121*t115*t183*t21/64.0;
	Real t191 = t1*t183;
	Real t194 = t121*t191*t37/64.0;
	Real t195 = t191*t6;
	Real t197 = t130*t195/64.0;
	Real t202 = t153*t195/64.0;
	Real t209 = t165*t160*t183*t18/64.0;

	// Add this contribution to the local mass matrix
	_localMass(0, 0) += weight * t113*t116*t117/64.0;
	_localMass(0, 1) += weight * t125;
	_localMass(0, 2) += weight * t129;
	_localMass(0, 3) += weight * t133;
	_localMass(0, 4) += weight * t137;
	_localMass(0, 5) += weight * t140;
	_localMass(0, 6) += weight * t144;
	_localMass(0, 7) += weight * t147;
	_localMass(1, 0) += weight * t125;
	_localMass(1, 1) += weight * t113*t149*t117/64.0;
	_localMass(1, 2) += weight * t155;
	_localMass(1, 3) += weight * t129;
	_localMass(1, 4) += weight * t140;
	_localMass(1, 5) += weight * t157;
	_localMass(1, 6) += weight * t159;
	_localMass(1, 7) += weight * t144;
	_localMass(2, 0) += weight * t129;
	_localMass(2, 1) += weight * t155;
	_localMass(2, 2) += weight * t113*t161*t117/64.0;
	_localMass(2, 3) += weight * t169;
	_localMass(2, 4) += weight * t144;
	_localMass(2, 5) += weight * t159;
	_localMass(2, 6) += weight * t173;
	_localMass(2, 7) += weight * t176;
	_localMass(3, 0) += weight * t133;
	_localMass(3, 1) += weight * t129;
	_localMass(3, 2) += weight * t169;
	_localMass(3, 3) += weight * t113*t177*t117/64.0;
	_localMass(3, 4) += weight * t147;
	_localMass(3, 5) += weight * t144;
	_localMass(3, 6) += weight * t176;
	_localMass(3, 7) += weight * t182;
	_localMass(4, 0) += weight * t137;
	_localMass(4, 1) += weight * t140;
	_localMass(4, 2) += weight * t144;
	_localMass(4, 3) += weight * t147;
	_localMass(4, 4) += weight * t113*t116*t183/64.0;
	_localMass(4, 5) += weight * t190;
	_localMass(4, 6) += weight * t194;
	_localMass(4, 7) += weight * t197;
	_localMass(5, 0) += weight * t140;
	_localMass(5, 1) += weight * t157;
	_localMass(5, 2) += weight * t159;
	_localMass(5, 3) += weight * t144;
	_localMass(5, 4) += weight * t190;
	_localMass(5, 5) += weight * t113*t149*t183/64.0;
	_localMass(5, 6) += weight * t202;
	_localMass(5, 7) += weight * t194;
	_localMass(6, 0) += weight * t144;
	_localMass(6, 1) += weight * t159;
	_localMass(6, 2) += weight * t173;
	_localMass(6, 3) += weight * t176;
	_localMass(6, 4) += weight * t194;
	_localMass(6, 5) += weight * t202;
	_localMass(6, 6) += weight * t113*t161*t183/64.0;
	_localMass(6, 7) += weight * t209;
	_localMass(7, 0) += weight * t147;
	_localMass(7, 1) += weight * t144;
	_localMass(7, 2) += weight * t176;
	_localMass(7, 3) += weight * t182;
	_localMass(7, 4) += weight * t197;
	_localMass(7, 5) += weight * t194;
	_localMass(7, 6) += weight * t209;
	_localMass(7, 7) += weight * t113*t177*t183/64.0;
}

// Evaluates strain at the given xi, eta, zeta coordinates
void IsoHexElement::evaluateStrain( VEC3F &p, VECTOR &meshState )
{
	setLocalState( meshState );

	Real xi = p[0];
	Real eta = p[1];
	Real zeta_ = p[2];

	// Define the strain displacement matrix
	Real t1 = 1.0-eta;
	Real t2 = 1.0-zeta_;
	Real t3 = t1*t2;
	Real t4 = 1.0-xi;
	Real t5 = t4*t2;
	Real t7 = 1.0+xi;
	Real t8 = t7*t2;
	Real t12 = 1.0+zeta_;
	Real t13 = t4*t12;
	Real t15 = t7*t12;
	Real t19 = -t5*y0-t8*y1+t8*y2+t5*y3-t13*y4-t15*y5+t15*y6+t13*y7;
	Real t20 = t4*t1;
	Real t22 = t7*t1;
	Real t24 = 1.0+eta;
	Real t25 = t7*t24;
	Real t27 = t4*t24;
	Real t33 = -t20*z0-t22*z1-t25*z2-t27*z3+t20*z4+t22*z5+t25*z6+t27*z7;
	Real t43 = -t20*y0-t22*y1-t25*y2-t27*y3+t20*y4+t22*y5+t25*y6+t27*y7;
	Real t52 = -t5*z0-t8*z1+t8*z2+t5*z3-t13*z4-t15*z5+t15*z6+t13*z7;
	Real t54 = t19*t33/64.0-t43*t52/64.0;
	Real t55 = t3*t54;
	Real t58 = t24*t2;
	Real t61 = t1*t12;
	Real t64 = t24*t12;
	Real t67 = -t3*z0+t3*z1+t58*z2-t58*z3-t61*z4+t61*z5+t64*z6-t64*z7;
	Real t77 = -t3*y0+t3*y1+t58*y2-t58*y3-t61*y4+t61*y5+t64*y6-t64*y7;
	Real t79 = t43*t67/64.0-t77*t33/64.0;
	Real t80 = t5*t79;
	Real t83 = t77*t52/64.0-t19*t67/64.0;
	Real t84 = t20*t83;
	Real t94 = -t3*x0+t3*x1+t58*x2-t58*x3-t61*x4+t61*x5+t64*x6-t64*x7;
	Real t104 = -t5*x0-t8*x1+t8*x2+t5*x3-t13*x4-t15*x5+t15*x6+t13*x7;
	Real t114 = -t20*x0-t22*x1-t25*x2-t27*x3+t20*x4+t22*x5+t25*x6+t27*x7;
	Real t117 = 1/(t94*t54/8.0+t104*t79/8.0+t114*t83/8.0);
	Real t118 = (-t55-t80-t84)*t117/8.0;
	Real t119 = t8*t79;
	Real t120 = t22*t83;
	Real t122 = (t55-t119-t120)*t117/8.0;
	Real t123 = t58*t54;
	Real t124 = t25*t83;
	Real t126 = (t123+t119-t124)*t117/8.0;
	Real t127 = t27*t83;
	Real t129 = (-t123+t80-t127)*t117/8.0;
	Real t130 = t61*t54;
	Real t131 = t13*t79;
	Real t133 = (-t130-t131+t84)*t117/8.0;
	Real t134 = t15*t79;
	Real t136 = (t130-t134+t120)*t117/8.0;
	Real t137 = t64*t54;
	Real t139 = (t137+t134+t124)*t117/8.0;
	Real t141 = (-t137+t131+t127)*t117/8.0;
	Real t144 = t52*t114/64.0-t33*t104/64.0;
	Real t145 = t3*t144;
	Real t148 = t33*t94/64.0-t67*t114/64.0;
	Real t149 = t5*t148;
	Real t152 = t67*t104/64.0-t52*t94/64.0;
	Real t153 = t20*t152;
	Real t155 = (-t145-t149-t153)*t117/8.0;
	Real t156 = t8*t148;
	Real t157 = t22*t152;
	Real t159 = (t145-t156-t157)*t117/8.0;
	Real t160 = t58*t144;
	Real t161 = t25*t152;
	Real t163 = (t160+t156-t161)*t117/8.0;
	Real t164 = t27*t152;
	Real t166 = (-t160+t149-t164)*t117/8.0;
	Real t167 = t61*t144;
	Real t168 = t13*t148;
	Real t170 = (-t167-t168+t153)*t117/8.0;
	Real t171 = t15*t148;
	Real t173 = (t167-t171+t157)*t117/8.0;
	Real t174 = t64*t144;
	Real t176 = (t174+t171+t161)*t117/8.0;
	Real t178 = (-t174+t168+t164)*t117/8.0;
	Real t181 = t104*t43/64.0-t114*t19/64.0;
	Real t182 = t3*t181;
	Real t185 = t114*t77/64.0-t94*t43/64.0;
	Real t186 = t5*t185;
	Real t189 = t94*t19/64.0-t104*t77/64.0;
	Real t190 = t20*t189;
	Real t192 = (-t182-t186-t190)*t117/8.0;
	Real t193 = t8*t185;
	Real t194 = t22*t189;
	Real t196 = (t182-t193-t194)*t117/8.0;
	Real t197 = t58*t181;
	Real t198 = t25*t189;
	Real t200 = (t197+t193-t198)*t117/8.0;
	Real t201 = t27*t189;
	Real t203 = (-t197+t186-t201)*t117/8.0;
	Real t204 = t61*t181;
	Real t205 = t13*t185;
	Real t207 = (-t204-t205+t190)*t117/8.0;
	Real t208 = t15*t185;
	Real t210 = (t204-t208+t194)*t117/8.0;
	Real t211 = t64*t181;
	Real t213 = (t211+t208+t198)*t117/8.0;
	Real t215 = (-t211+t205+t201)*t117/8.0;

	MATRIX strainDisplacement(6, 24);

	strainDisplacement(0, 0) = t118;
	strainDisplacement(0, 1) = 0.0;
	strainDisplacement(0, 2) = 0.0;
	strainDisplacement(0, 3) = t122;
	strainDisplacement(0, 4) = 0.0;
	strainDisplacement(0, 5) = 0.0;
	strainDisplacement(0, 6) = t126;
	strainDisplacement(0, 7) = 0.0;
	strainDisplacement(0, 8) = 0.0;
	strainDisplacement(0, 9) = t129;
	strainDisplacement(0, 10) = 0.0;
	strainDisplacement(0, 11) = 0.0;
	strainDisplacement(0, 12) = t133;
	strainDisplacement(0, 13) = 0.0;
	strainDisplacement(0, 14) = 0.0;
	strainDisplacement(0, 15) = t136;
	strainDisplacement(0, 16) = 0.0;
	strainDisplacement(0, 17) = 0.0;
	strainDisplacement(0, 18) = t139;
	strainDisplacement(0, 19) = 0.0;
	strainDisplacement(0, 20) = 0.0;
	strainDisplacement(0, 21) = t141;
	strainDisplacement(0, 22) = 0.0;
	strainDisplacement(0, 23) = 0.0;
	strainDisplacement(1, 0) = 0.0;
	strainDisplacement(1, 1) = t155;
	strainDisplacement(1, 2) = 0.0;
	strainDisplacement(1, 3) = 0.0;
	strainDisplacement(1, 4) = t159;
	strainDisplacement(1, 5) = 0.0;
	strainDisplacement(1, 6) = 0.0;
	strainDisplacement(1, 7) = t163;
	strainDisplacement(1, 8) = 0.0;
	strainDisplacement(1, 9) = 0.0;
	strainDisplacement(1, 10) = t166;
	strainDisplacement(1, 11) = 0.0;
	strainDisplacement(1, 12) = 0.0;
	strainDisplacement(1, 13) = t170;
	strainDisplacement(1, 14) = 0.0;
	strainDisplacement(1, 15) = 0.0;
	strainDisplacement(1, 16) = t173;
	strainDisplacement(1, 17) = 0.0;
	strainDisplacement(1, 18) = 0.0;
	strainDisplacement(1, 19) = t176;
	strainDisplacement(1, 20) = 0.0;
	strainDisplacement(1, 21) = 0.0;
	strainDisplacement(1, 22) = t178;
	strainDisplacement(1, 23) = 0.0;
	strainDisplacement(2, 0) = 0.0;
	strainDisplacement(2, 1) = 0.0;
	strainDisplacement(2, 2) = t192;
	strainDisplacement(2, 3) = 0.0;
	strainDisplacement(2, 4) = 0.0;
	strainDisplacement(2, 5) = t196;
	strainDisplacement(2, 6) = 0.0;
	strainDisplacement(2, 7) = 0.0;
	strainDisplacement(2, 8) = t200;
	strainDisplacement(2, 9) = 0.0;
	strainDisplacement(2, 10) = 0.0;
	strainDisplacement(2, 11) = t203;
	strainDisplacement(2, 12) = 0.0;
	strainDisplacement(2, 13) = 0.0;
	strainDisplacement(2, 14) = t207;
	strainDisplacement(2, 15) = 0.0;
	strainDisplacement(2, 16) = 0.0;
	strainDisplacement(2, 17) = t210;
	strainDisplacement(2, 18) = 0.0;
	strainDisplacement(2, 19) = 0.0;
	strainDisplacement(2, 20) = t213;
	strainDisplacement(2, 21) = 0.0;
	strainDisplacement(2, 22) = 0.0;
	strainDisplacement(2, 23) = t215;
	strainDisplacement(3, 0) = 0.0;
	strainDisplacement(3, 1) = t192;
	strainDisplacement(3, 2) = t155;
	strainDisplacement(3, 3) = 0.0;
	strainDisplacement(3, 4) = t196;
	strainDisplacement(3, 5) = t159;
	strainDisplacement(3, 6) = 0.0;
	strainDisplacement(3, 7) = t200;
	strainDisplacement(3, 8) = t163;
	strainDisplacement(3, 9) = 0.0;
	strainDisplacement(3, 10) = t203;
	strainDisplacement(3, 11) = t166;
	strainDisplacement(3, 12) = 0.0;
	strainDisplacement(3, 13) = t207;
	strainDisplacement(3, 14) = t170;
	strainDisplacement(3, 15) = 0.0;
	strainDisplacement(3, 16) = t210;
	strainDisplacement(3, 17) = t173;
	strainDisplacement(3, 18) = 0.0;
	strainDisplacement(3, 19) = t213;
	strainDisplacement(3, 20) = t176;
	strainDisplacement(3, 21) = 0.0;
	strainDisplacement(3, 22) = t215;
	strainDisplacement(3, 23) = t178;
	strainDisplacement(4, 0) = t192;
	strainDisplacement(4, 1) = 0.0;
	strainDisplacement(4, 2) = t118;
	strainDisplacement(4, 3) = t196;
	strainDisplacement(4, 4) = 0.0;
	strainDisplacement(4, 5) = t122;
	strainDisplacement(4, 6) = t200;
	strainDisplacement(4, 7) = 0.0;
	strainDisplacement(4, 8) = t126;
	strainDisplacement(4, 9) = t203;
	strainDisplacement(4, 10) = 0.0;
	strainDisplacement(4, 11) = t129;
	strainDisplacement(4, 12) = t207;
	strainDisplacement(4, 13) = 0.0;
	strainDisplacement(4, 14) = t133;
	strainDisplacement(4, 15) = t210;
	strainDisplacement(4, 16) = 0.0;
	strainDisplacement(4, 17) = t136;
	strainDisplacement(4, 18) = t213;
	strainDisplacement(4, 19) = 0.0;
	strainDisplacement(4, 20) = t139;
	strainDisplacement(4, 21) = t215;
	strainDisplacement(4, 22) = 0.0;
	strainDisplacement(4, 23) = t141;
	strainDisplacement(5, 0) = t155;
	strainDisplacement(5, 1) = t118;
	strainDisplacement(5, 2) = 0.0;
	strainDisplacement(5, 3) = t159;
	strainDisplacement(5, 4) = t122;
	strainDisplacement(5, 5) = 0.0;
	strainDisplacement(5, 6) = t163;
	strainDisplacement(5, 7) = t126;
	strainDisplacement(5, 8) = 0.0;
	strainDisplacement(5, 9) = t166;
	strainDisplacement(5, 10) = t129;
	strainDisplacement(5, 11) = 0.0;
	strainDisplacement(5, 12) = t170;
	strainDisplacement(5, 13) = t133;
	strainDisplacement(5, 14) = 0.0;
	strainDisplacement(5, 15) = t173;
	strainDisplacement(5, 16) = t136;
	strainDisplacement(5, 17) = 0.0;
	strainDisplacement(5, 18) = t176;
	strainDisplacement(5, 19) = t139;
	strainDisplacement(5, 20) = 0.0;
	strainDisplacement(5, 21) = t178;
	strainDisplacement(5, 22) = t141;
	strainDisplacement(5, 23) = 0.0;

	// Strain is given by the product of this matrix with the
	// displacements for this element
	strainDisplacement.multiplyInplace( _localState, _strain );
}
