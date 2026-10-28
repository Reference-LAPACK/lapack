if(BUILD_DEFAULT_API)
  add_executable(xexample1_CBLAS cblas_example1.c)
  lapack_configure_target(xexample1_CBLAS)
  add_executable(xexample2_CBLAS cblas_example2.c)
  lapack_configure_target(xexample2_CBLAS)

  target_link_libraries(xexample1_CBLAS PRIVATE ${CBLASLIB})
  target_link_libraries(xexample2_CBLAS PRIVATE ${CBLASLIB} ${lapack_blas_target})

  add_test(NAME example1_CBLAS COMMAND xexample1_CBLAS)
  add_test(NAME example2_CBLAS COMMAND xexample2_CBLAS)
endif()

if(BUILD_INDEX64_EXT_API)
  add_executable(xexample1_64_CBLAS cblas_example1_64.c)
  lapack_configure_target(xexample1_64_CBLAS)
  add_executable(xexample2_64_CBLAS cblas_example2_64.c)
  lapack_configure_target(xexample2_64_CBLAS)

  target_link_libraries(xexample1_64_CBLAS PRIVATE ${CBLASLIB})
  target_link_libraries(xexample2_64_CBLAS PRIVATE ${CBLASLIB} ${lapack_blas_target})

  add_test(NAME example1_64_CBLAS COMMAND xexample1_64_CBLAS)
  add_test(NAME example2_64_CBLAS COMMAND xexample2_64_CBLAS)
endif()
