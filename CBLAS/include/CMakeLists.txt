set(CBLAS_INCLUDE cblas.h cblas_f77.h)

file(COPY ${CBLAS_INCLUDE} DESTINATION ${LAPACK_BINARY_DIR}/include)

file(COPY cblas_test.h DESTINATION ${LAPACK_BINARY_DIR}/test_include)
