SET (CBLAS_INCLUDE cblas.h cblas_f77.h cblas_test.h cblas_mangling.h)
message(STATUS "Include CBLAS_INCLUDE = ${CBLAS_INCLUDE}")

file(COPY ${CBLAS_INCLUDE} DESTINATION ${LAPACK_BINARY_DIR}/include)