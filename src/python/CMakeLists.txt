include_directories(.)

set_source_files_properties(rschol.pyx
  PROPERTIES CYTHON_IS_CXX TRUE)

cython_add_module(rschol rschol.pyx pcg_pysolver.cpp)

target_link_libraries(rschol rsclib
  ${BLAS_LAPACK_LIBS}
  ${CHOLMOD_LIBRARIES}
  ${ZLIB_LIBRARIES})

install(TARGETS rschol DESTINATION lib/python)
