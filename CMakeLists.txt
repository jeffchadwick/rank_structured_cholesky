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

# Look for Deal.II
option(BUILD_WITH_DEALII "Build Deal.II examples" OFF)

if ( BUILD_WITH_DEALII )
  find_package(Deal.II 7.2 REQUIRED)
endif ( BUILD_WITH_DEALII )

# Look for CHOLMOD
find_package(Cholmod REQUIRED)

# DSB: Use librt?
option(USE_LIBRT "Build with librt" OFF)
if ( USE_LIBRT )
  set(RT_LIBRARY rt)
endif ( USE_LIBRT )

# DSB: Use zlib?
option(USE_ZLIB "Built with zlib" ON)
if (USE_ZLIB)
  find_package(ZLIB REQUIRED)
  set(HAS_ZLIB 1)
endif (USE_ZLIB)

#===================================================================
# Enable/disable compile targets
#===================================================================
option(BUILD_PCG_SOLVER
       "Compile PCG solver"                               ON)

# Deal.II's nonlinear elasticity tutorial example, modified
# to use our solver
option(BUILD_NONLINEAR_ELASTICITY
       "Compile nonlinear_elasticity"                     OFF)

# Program for testing our Trilinos solver interface against other
# preconditioners provided by Trilinos
option(BUILD_TRILINOS_SOLVER
       "Compile trilinos_solver"                          OFF)


#===================================================================
# If we are building certain executables, we need to make sure that
# we have include and library directories for 
#===================================================================
if ( BUILD_NONLINEAR_ELASTICITY AND ( NOT TRILINOS_INCLUDE_DIR
                                      OR NOT TRILINOS_LIB_DIR ) )

  MESSAGE(FATAL_ERROR "Need to specify TRILINOS_INCLUDE_DIR and TRILINOS_LIB_DIR to build nonlinear_elasticity")

endif ()

#===================================================================
# Set up the TRILINOS library list.  Note that this build system
# does not explicitly check for the existence of all of these
# packages.
#
# TODO: do we need python2.7 here?
#===================================================================
if ( TRILINOS_INCLUDE_DIR AND TRILINOS_LIB_DIR )
  # FIXME: Need to figure out how to set the library list correctly
  #set(TRILINOS_LIBS epetra aztecoo teuchoscore teuchosparameterlist ifpack ml)

  if ( NOT TEUCHOS_LIBS )
    MESSAGE(FATAL_ERROR "No teuchos libraries specified")
  endif ( NOT TEUCHOS_LIBS )

  set(TRILINOS_LIBS epetra aztecoo ${TEUCHOS_LIBS} ifpack ml)
endif ( TRILINOS_INCLUDE_DIR AND TRILINOS_LIB_DIR )


#===================================================================
# Depending on how you built Trilinos and/or Deal.II, you may need
# to also include/link against an MPI implementation
#
# Enable or disable this here
#===================================================================
option(BUILD_WITH_MPI  "Include and link against MPI"   OFF)

if ( BUILD_WITH_MPI )
  find_package(MPI REQUIRED)
endif ( BUILD_WITH_MPI )


#===================================================================
# Similarly, we may need to include and link against PetSC,
# depending on how Deal.II was built
#===================================================================
option(BUILD_WITH_PETSC   "Include and link against PetSC"  OFF)

if ( BUILD_WITH_PETSC )
  find_package(PETSc REQUIRED)
endif ( BUILD_WITH_PETSC )


#===================================================================
# Set up the list of libraries necessary for deal.ii
#
# TODO: We might also need tbb, netcdf, netcdf_c++ and libz (ugh)
#===================================================================
if ( Deal.II_LIBRARY )
  # We also need to link against ZLIB for Deal.II
  find_package(ZLIB REQUIRED)

  set(Deal.II_LIBS deal_II tbb)
endif ( Deal.II_LIBRARY )

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
  ${ZLIB_INCLUDE_DIRS}
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
                        ${BLAS_LAPACK_LIBS}       ${CHOLMOD_LIBRARIES}
                        ${ZLIB_LIBRARIES}         ${RT_LIBRARY})
endif ( BUILD_PCG_SOLVER )

if ( BUILD_NONLINEAR_ELASTICITY )
  # Requires CHOLMOD (for our solver) as well as Trilinos and Deal.II
  add_executable(nonlinear_elasticity
                 deal.ii/nonlinear_elasticity/step-44.cpp
                 ${UTIL_SRC}              ${LINEAR_SRC}
                 ${SOLVER_SRC}            ${ORDERING_SRC}
                 ${NODE_SRC}              ${GEOMETRY_SRC}
                 ${CHOLMOD_INTERFACE_SRC} ${CSPARSE_INTERFACE_SRC}
                 ${Deal.II_UTIL_SRC}
                 ${TRILINOS_UTIL_SRC}     ${TRILINOS_SOLVER_SRC})

  target_link_libraries(nonlinear_elasticity
                        ${BLAS_LAPACK_LIBS}   ${CHOLMOD_LIBRARIES}
                        ${Deal.II_LIBS}       ${TRILINOS_LIBS}
                        ${MPI_CXX_LIBRARIES}  ${PETSC_LIBRARIES}
                        ${ZLIB_LIBRARIES}     
                        ${RT_LIBRARY})

endif ( BUILD_NONLINEAR_ELASTICITY )

if ( BUILD_TRILINOS_SOLVER )
  add_executable(trilinos_solver
                 trilinos/testML.cpp
                 ${UTIL_SRC}              ${LINEAR_SRC}
                 ${SOLVER_SRC}            ${ORDERING_SRC}
                 ${NODE_SRC}              ${GEOMETRY_SRC}
                 ${CHOLMOD_INTERFACE_SRC} ${CSPARSE_INTERFACE_SRC}
                 ${TRILINOS_UTIL_SRC}     ${TRILINOS_SOLVER_SRC})

  target_link_libraries(trilinos_solver
                        ${BLAS_LAPACK_LIBS}   ${CHOLMOD_LIBRARIES}
                        ${TRILINOS_LIBS}      ${MPI_CXX_LIBRARIES}
                        ${ZLIB_LIBRARIES}     
                        ${RT_LIBRARY})

endif ( BUILD_TRILINOS_SOLVER )
