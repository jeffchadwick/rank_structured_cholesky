# Find the Deal.II Finite Element Analysis library
FIND_PATH(Deal.II_INCLUDE_DIR deal.II
  PATHS ENV INCLUDE
  PATHS ${SYSTEM_INC_PATH}
  PATHS ${Deal.II_DIR}/include)

FIND_LIBRARY(Deal.II_LIBRARY deal_II
  PATHS ENV LD_LIBRARY_PATH
  PATHS ${SYSTEM_LIB_PATH}
  PATHS ENV LIBRARY_PATH
  PATHS ${Deal.II_DIR}/lib)

GET_FILENAME_COMPONENT(Deal.II_LIBRARY_DIR ${Deal.II_LIBRARY} PATH)

IF ( Deal.II_INCLUDE_DIR AND Deal.II_LIBRARY AND Deal.II_LIBRARY_DIR )
  MESSAGE(STATUS "Found Deal.II: ${Deal.II_LIBRARY_DIR}")
ELSE ( Deal.II_INCLUDE_DIR AND Deal.II_LIBRARY AND Deal.II_LIBRARY_DIR )
  IF ( Deal.II_FIND_REQUIRED )
    MESSAGE(FATAL_ERROR "Could not find Deal.II: $ENV{INCLUDE} $ENV{LD_LIBRARY_PATH}")
  ELSE ( Deal.II_FIND_REQUIRED )
    IF ( NOT Deal.II_LIBRARY )
      MESSAGE(STATUS "WARNING: Could not find Deal.II library: $ENV{LD_LIBRARY_PATH}")
    ELSE ( NOT Deal.II_LIBRARY )
      MESSAGE(STATUS "WARNING: Could not find MKL include: $ENV{INCLUDE}")
    ENDIF ( NOT Deal.II_LIBRARY )
  ENDIF ( Deal.II_FIND_REQUIRED )
ENDIF ( Deal.II_INCLUDE_DIR AND Deal.II_LIBRARY AND Deal.II_LIBRARY_DIR )
