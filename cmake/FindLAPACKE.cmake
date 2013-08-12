# - Try to find LAPACKE
#
# Once done this will define
#  LAPACKE_FOUND - System has LAPACKE
#  LAPACKE_INCLUDE_DIRS - The LAPACKE include directories
#  LAPACKE_LIBRARIES - The libraries needed to use LAPACKE
#  LAPACKE_DEFINITIONS - Compiler switches required for using LAPACKE
#
# Usually, LAPACKE requires LAPACK and the BLAS.  This module does
# not enforce anything about that.

find_path(LAPACKE_INCLUDE_DIR
          NAMES lapacke.h
          PATHS $ENV{LAPACKEDIR} ${INCLUDE_INSTALL_DIR}
          PATHS ENV INCLUDE)

find_library(LAPACKE_LIBRARY lapacke
             PATHS $ENV{LAPACKEDIR} ${LIB_INSTALL_DIR}
             PATHS ENV LIBRARY_PATH
             PATHS ENV LD_LIBRARY_PATH)

set(LAPACKE_LIBRARIES ${LAPACKE_LIBRARY})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(LAPACKE DEFAULT_MSG
                                  LAPACKE_INCLUDE_DIR 
                                  LAPACKE_LIBRARIES)
mark_as_advanced(LAPACKE_INCLUDE_DIR LAPACKE_LIBRARIES)
