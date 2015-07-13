// SETTINGS.h: Project-wide options set in one place
//
//////////////////////////////////////////////////////////////////////

#ifndef SETTINGS_H
#define SETTINGS_H

#include <cassert>

#include <boost/tuple/tuple.hpp>

#include <cstddef>

//#ifdef _WIN32
//#define USING_OPENMP 
//#define USING_RENDERMAN
//#define USING_MKL
//#endif

#define EPSILON 1e-16

// At least Intel MKL or ATLAS must be installed in order
// to compile and run "CubatureViewer". Select one and only one to use

// select single or double precision
#ifdef SINGLE_PRECISION
typedef float Real;
#else
typedef double Real;
#endif

// turn off asserts?
//#define NDEBUG

#endif
