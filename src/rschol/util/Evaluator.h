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
// Evaluator.h: Some useful function evaluator definitions
//
//////////////////////////////////////////////////////////////////////

#ifndef EVALUATOR_H
#define EVALUATOR_H

#include <rschol/linearalgebra/MATRIX3.h>
#include <rschol/linearalgebra/VECTOR.h>
#include <rschol/linearalgebra/MATRIX.h>

#include <rschol/SETTINGS.h>
#include <rschol/TYPES.h>

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
