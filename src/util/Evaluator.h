//////////////////////////////////////////////////////////////////////
// Evaluator.h: Some useful function evaluator definitions
//
//////////////////////////////////////////////////////////////////////

#ifndef EVALUATOR_H
#define EVALUATOR_H

#include <linearalgebra/MATRIX3.h>
#include <linearalgebra/VECTOR.h>
#include <linearalgebra/MATRIX.h>

#include <SETTINGS.h>
#include <TYPES.h>

#include <boost/bind.hpp>
#include <boost/function.hpp>

// Functions for evaluating a function from R^3 --> ?
// Types with the suffix 'E' are evaluators which return
// error information.  They return true in the case of a
// successful evaluation and false otherwise.
// Ref versions take the return value as a reference
typedef boost::function<VEC3F (const VEC3F &x)> Vec3Evaluator;
typedef boost::function<void (const VEC3F &x, VEC3F &result)> Vec3Evaluator_Ref;

typedef boost::function<VECTOR & (const VEC3F &x)> VectorEvaluator;
typedef boost::function<void (const VEC3F &x, VECTOR &result)> VectorEvaluator_Ref;

typedef boost::function<MATRIX3 & (const VEC3F &x)> Matrix3Evaluator;
typedef boost::function<void (const VEC3F &x, MATRIX3 &result)> Matrix3Evalutor_Ref;

typedef boost::function<MATRIX & (const VEC3F &x)> MatrixEvaluator;
typedef boost::function<void (const VEC3F &x, MATRIX &result)> MatrixEvaluator_Ref;

typedef boost::function<Real (const VEC3F &x)> ScalarEvaluator;
typedef boost::function<void (const VEC3F &x, Real &result)> ScalarEvaluator_Ref;
typedef boost::function<bool (const VEC3F &x, Real &value)> ScalarEvaluator_E;

typedef boost::function<Real (Real)> RealFunction;

#endif
