# Copyright (c) 2015, Jeff Chadwick and David Bindel
# All rights reserved.
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
# 
# 1. Redistributions of source code must retain the above copyright notice,
# this list of conditions and the following disclaimer.
# 
# 2. Redistributions in binary form must reproduce the above copyright notice,
# this list of conditions and the following disclaimer in the documentation
# and/or other materials provided with the distribution.
# 
# 3. Neither the name of the copyright holder nor the names of its
# contributors may be used to endorse or promote products derived from this
# software without specific prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
# THE POSSIBILITY OF SUCH DAMAGE.

# CMake project file for the rank-structured Cholesky preconditioner
# Jeffrey Chadwick (chadwick.jeff@gmail.com)
#
# DSB: Switched OpenMP support

project (rank_structured_cholesky)

cmake_minimum_required(VERSION 2.8)

# General Configuration: find all the required libraries.
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/cmake-modules-master")

# compile option
option(USE_DEBUG "Turn on the debug mode" OFF)
option(USE_OPENMP "Turn on the OpenMP feature when compiling" ON)
#option(USE_64BIT_ARCH "Compile the 64bit executable" OFF)

#===================================================================
# RPATH management

SET(MACOSX_RPATH TRUE)
SET(CMAKE_MACOSX_RPATH TRUE)
SET(CMAKE_SKIP_BUILD_RPATH  FALSE)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE) 
SET(CMAKE_INSTALL_NAME_DIR "${CMAKE_INSTALL_PREFIX}/lib")
SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
# the RPATH to be used when installing, but only if it's not a system directory
LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
IF("${isSystemDir}" STREQUAL "-1")
   SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
ENDIF("${isSystemDir}" STREQUAL "-1")
 

#===================================================================
# Compilation options

set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_CXX_STANDARD 11)

#===================================================================
## Compiler
# set compiler flags for debug/release

if ( USE_DEBUG )
  MESSAGE(STATUS "Using debug mode")
  add_definitions(-DDEBUG)
  set(CMAKE_BUILD_TYPE Debug)
else ( USE_DEBUG )
  MESSAGE(STATUS "Using release mode")
  set(CMAKE_BUILD_TYPE Release)
endif ( USE_DEBUG )

if ( USE_OPENMP )
  find_package(OpenMP)
  if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  else (OPENMP_FOUND)
    MESSAGE(STATUS "Could not determine OpenMP support, disabling")
    set( USE_OPENMP OFF )
  endif()
endif ( USE_OPENMP )

# add_definitions("-std=c++11")

# check if we're using Intel's compiler
if ( CMAKE_CXX_COMPILER MATCHES ".*icpc$" )
  add_definitions(-wd981 -wd383 -wd444 -wd1224 -wd1572)
  if ( NOT USE_DEBUG )
    set(CMAKE_CXX_FLAGS_RELEASE "-ipo -O3 -no-prec-div -xHost -DNDEBUG")
  endif ( NOT USE_DEBUG )

  # Get around builtin fpclassify bug
  #add_definitions("-gcc-name=gcc-4.4")
  #add_definitions("-gxx-name=g++-4.4")
else ( CMAKE_CXX_COMPILER MATCHES ".*icpc$" )
  if ( USE_OPENMP )
    add_definitions(-fopenmp)
  endif ( USE_OPENMP )

  if ( NOT USE_DEBUG )
    set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
  endif ( NOT USE_DEBUG )
endif ( CMAKE_CXX_COMPILER MATCHES ".*icpc$" )

#if ( USE_64BIT_ARCH )
#    add_definitions(-m64)
#    set(LINK_FLAGS -m64)
#endif ( USE_64BIT_ARCH )

# Check C++11 features
include(CheckCXX11Features)

#===================================================================
## Libraries
#===================================================================

# check boost libraries
find_package(Boost 1.33 REQUIRED)

# BLAS/LAPACK implementation
option(USE_MKL  "Use Intel MKL BLAS/LAPACK"   ON)

if ( USE_MKL )

  # check for Intel MKL
  find_package(MKL REQUIRED)
  set(BLAS_LAPACK_LIBS mkl_intel_lp64 mkl_intel_thread mkl_core iomp5)
  set(BLAS_INCLUDE_DIR ${MKL_INCLUDE_DIR})
  set(LAPACK_INCLUDE_DIR ${MKL_INCLUDE_DIR})
  add_definitions(-DUSE_MKL)

else ( USE_MKL )

  # Look for standard BLAS/LAPACK implementations.
  find_package(CBLAS REQUIRED)
  find_package(LAPACK REQUIRED)
  find_package(LAPACKE REQUIRED)

  set(BLAS_LAPACK_LIBS ${LAPACKE_LIBRARY}
      ${LAPACK_LIBRARIES} ${CBLAS_LIBRARIES})
  set(BLAS_INCLUDE_DIR ${CBLAS_INCLUDE_DIR})
  set(LAPACK_INCLUDE_DIR ${LAPACKE_INCLUDE_DIR})

endif ( USE_MKL )

# Look for CHOLMOD
find_package(Cholmod REQUIRED)

#===================================================================
# Enable/disable compile targets
#===================================================================

option(BUILD_PCG_SOLVER
       "Compile PCG solver"                               ON)

option(BUILD_PCG_PYTHON
       "Compile PCG Python"                               OFF)

option(USE_ZLIB
       "Build with zlib"                                  ON)

if ( USE_ZLIB )
  FIND_PACKAGE(ZLIB REQUIRED)
  SET(HAS_ZLIB 1)
endif ( USE_ZLIB )

#===================================================================
# Source lists
#===================================================================

#===================================================================
# Linear algebra
#===================================================================

set(LINEAR_SRC
    src/rschol/linearalgebra/MATRIX3.cpp
    src/rschol/linearalgebra/MATRIX_FAST.cpp
    src/rschol/linearalgebra/PCG_SOLVER.cpp
    src/rschol/linearalgebra/SPARSE_MATRIX.cpp
    src/rschol/linearalgebra/VECTOR_FAST.cpp)

#===================================================================
# General utilities
#===================================================================

set(UTIL_SRC
    src/rschol/util/MathUtil.cpp                     src/rschol/util/StatsCounter.cpp)


#===================================================================
# Solver sources
#===================================================================

set(SOLVER_SRC
    src/rschol/solver/FactorManager.cpp
    src/rschol/solver/FactorPreconditioner.cpp       src/rschol/solver/SolverError.cpp
    src/rschol/solver/SparseSolver.cpp               src/rschol/solver/FactorPreconditioner.cpp)


#===================================================================
# Interface sources for CHOLMOD and CSparse
#===================================================================

set(CHOLMOD_INTERFACE_SRC
    src/rschol/library_interface/CHOLMOD_Environment.cpp
    src/rschol/library_interface/CHOLMOD_Supernodal.cpp)

set(CSPARSE_INTERFACE_SRC
    src/rschol/library_interface/CSparse_Interface.cpp)


#===================================================================
# Symbolic factorization/reordering sources (depends on CHOLMOD
# and CSparse)
#===================================================================

set(ORDERING_SRC
    src/rschol/ordering/MinimumDegree.cpp
    src/rschol/ordering/NestedDissection.cpp
    src/rschol/ordering/Ordering.cpp)


#===================================================================
# Supernode sources
#===================================================================

set(NODE_SRC
    src/rschol/node/DiagonalBlock.cpp                src/rschol/node/Supernode.cpp)


#===================================================================
# Geometry sources
#===================================================================

set(GEOMETRY_SRC
    src/rschol/geometry/FiniteDifference.cpp)


#===================================================================
# Set up general includes
#===================================================================

include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}/src/
  ${CMAKE_CURRENT_BINARY_DIR}/src/
  ${Boost_INCLUDE_DIR}
  ${ZLIB_INCLUDE_DIRS}
  ${BLAS_INCLUDE_DIR}
  ${LAPACK_INCLUDE_DIR}
  ${CHOLMOD_INCLUDE_DIR})

link_directories(${Boost_LIBRARY_DIRS}
                 ${CHOLMOD_LIBDIR}
                 ${MKL_LIBRARY_DIR})

#===================================================================
# Configure and install header files
#===================================================================

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/src/rschol/config.h.in
               ${CMAKE_CURRENT_BINARY_DIR}/src/rschol/config.h)

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/src/rschol/config.h
        DESTINATION include/rschol)
install(DIRECTORY src/rschol/ 
        DESTINATION include/rschol
        FILES_MATCHING PATTERN "*.h")
install(DIRECTORY src/rschol/datastructure/ 
        DESTINATION include/rschol/datastructure
        FILES_MATCHING PATTERN "*.cpp")

#===================================================================
# Build libraries
#===================================================================

add_library(rsclib STATIC
  ${UTIL_SRC}              ${LINEAR_SRC}
  ${SOLVER_SRC}            ${ORDERING_SRC}
  ${NODE_SRC}              ${GEOMETRY_SRC}
  ${CHOLMOD_INTERFACE_SRC} ${CSPARSE_INTERFACE_SRC})

target_link_libraries(rsclib
  ${BLAS_LAPACK_LIBS}
  ${CHOLMOD_LIBRARIES}
  ${ZLIB_LIBRARIES})

install(TARGETS rsclib ARCHIVE DESTINATION lib)

#===================================================================
# Build executables
#===================================================================

if ( BUILD_PCG_SOLVER )

  add_executable(pcg_solver src/pcg_solver.cpp)

  target_link_libraries(pcg_solver rsclib
    ${BLAS_LAPACK_LIBS}
    ${CHOLMOD_LIBRARIES}
    ${ZLIB_LIBRARIES})

  install(TARGETS pcg_solver RUNTIME DESTINATION bin)

endif ( BUILD_PCG_SOLVER )

#===================================================================
# Build Cython interface
#===================================================================

if ( BUILD_PCG_PYTHON )
  include(UseCython)
  find_package(NumPy)
  include_directories(${PYTHON_NUMPY_INCLUDE_DIR})
  add_subdirectory(src/python)
endif ( BUILD_PCG_PYTHON )

