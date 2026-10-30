if(BUILD_DEFAULT_API)
  add_executable(xexample1_CBLAS cblas_example1.c)
  add_executable(xexample2_CBLAS cblas_example2.c)

  # Single-source examples gain nothing from unity builds.
  set_target_properties(xexample1_CBLAS xexample2_CBLAS PROPERTIES UNITY_BUILD OFF)

  target_link_libraries(xexample1_CBLAS ${CBLASLIB})
  target_link_libraries(xexample2_CBLAS ${CBLASLIB} ${BLAS_LIBRARIES})

  add_test(example1_CBLAS ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/xexample1_CBLAS)
  add_test(example2_CBLAS ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/xexample2_CBLAS)
endif()

if(BUILD_INDEX64_EXT_API)
  add_executable(xexample1_64_CBLAS cblas_example1_64.c)
  add_executable(xexample2_64_CBLAS cblas_example2_64.c)

  # Single-source examples gain nothing from unity builds.
  set_target_properties(xexample1_64_CBLAS xexample2_64_CBLAS PROPERTIES UNITY_BUILD OFF)

  target_link_libraries(xexample1_64_CBLAS ${CBLASLIB})
  target_link_libraries(xexample2_64_CBLAS ${CBLASLIB} ${BLAS_LIBRARIES})

  add_test(example1_64_CBLAS ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/xexample1_64_CBLAS)
  add_test(example2_64_CBLAS ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/xexample2_64_CBLAS)
endif()
