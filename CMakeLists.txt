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
    MESSAGE(FATAL_ERROR "Could not determine OpenMP support")
  endif()
endif ( USE_OPENMP )

add_definitions("-std=c++0x")

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

#===================================================================
# Source lists
#===================================================================

#===================================================================
# Linear algebra
#===================================================================
set(LINEAR_SRC
    src/linearalgebra/MATRIX3.cpp
    src/linearalgebra/MATRIX_FAST.cpp
    src/linearalgebra/PCG_SOLVER.cpp
    src/linearalgebra/SPARSE_MATRIX.cpp
    src/linearalgebra/VECTOR_FAST.cpp)

#===================================================================
# General utilities
#===================================================================
set(UTIL_SRC
    src/util/MathUtil.cpp                     src/util/StatsCounter.cpp)


#===================================================================
# Solver sources
#===================================================================
set(SOLVER_SRC
    src/solver/FactorManager.cpp
    src/solver/FactorPreconditioner.cpp       src/solver/SolverError.cpp
    src/solver/SparseSolver.cpp               src/solver/FactorPreconditioner.cpp)


#===================================================================
# Interface sources for CHOLMOD and CSparse
#===================================================================
set(CHOLMOD_INTERFACE_SRC
    src/library_interface/CHOLMOD_Environment.cpp
    src/library_interface/CHOLMOD_Supernodal.cpp)

set(CSPARSE_INTERFACE_SRC
    src/library_interface/CSparse_Interface.cpp)


#===================================================================
# Symbolic factorization/reordering sources (depends on CHOLMOD
# and CSparse)
#===================================================================
set(ORDERING_SRC
    src/ordering/MinimumDegree.cpp
    src/ordering/NestedDissection.cpp
    src/ordering/Ordering.cpp)


#===================================================================
# Supernode sources
#===================================================================
set(NODE_SRC
    src/node/DiagonalBlock.cpp                src/node/Supernode.cpp)


#===================================================================
# Geometry sources
#===================================================================
set(GEOMETRY_SRC
    src/geometry/FiniteDifference.cpp)


#===================================================================
# Set up general includes
#===================================================================
include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}/src
  ${CMAKE_CURRENT_BINARY_DIR}/src
  ${Boost_INCLUDE_DIR}
  ${BLAS_INCLUDE_DIR}
  ${LAPACK_INCLUDE_DIR}
  ${CHOLMOD_INCLUDE_DIR})
link_directories(${Boost_LIBRARY_DIRS}
                 ${CHOLMOD_LIBDIR}
                 ${MKL_LIBRARY_DIR})

#===================================================================
# Configure header files
#===================================================================

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/src/config.h.in 
               ${CMAKE_CURRENT_BINARY_DIR}/src/config.h)

#===================================================================
# Build executables
#===================================================================
if ( BUILD_PCG_SOLVER )
  add_executable(pcg_solver
                 src/pcg_solver.cpp
                 ${UTIL_SRC}              ${LINEAR_SRC}
                 ${SOLVER_SRC}            ${ORDERING_SRC}
                 ${NODE_SRC}              ${GEOMETRY_SRC}
                 ${CHOLMOD_INTERFACE_SRC} ${CSPARSE_INTERFACE_SRC})

  target_link_libraries(pcg_solver
                        ${BLAS_LAPACK_LIBS}       ${CHOLMOD_LIBRARIES})

endif ( BUILD_PCG_SOLVER )

