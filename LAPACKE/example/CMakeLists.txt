add_executable(xexample_DGESV_rowmajor example_DGESV_rowmajor.c lapacke_example_aux.c lapacke_example_aux.h)
add_executable(xexample_DGESV_colmajor example_DGESV_colmajor.c lapacke_example_aux.c lapacke_example_aux.h)
add_executable(xexample_DGELS_rowmajor example_DGELS_rowmajor.c lapacke_example_aux.c lapacke_example_aux.h)
add_executable(xexample_DGELS_colmajor example_DGELS_colmajor.c lapacke_example_aux.c lapacke_example_aux.h)

target_link_libraries(xexample_DGESV_rowmajor ${LAPACKELIB} ${BLAS_LIBRARIES})
target_link_libraries(xexample_DGESV_colmajor ${LAPACKELIB} ${BLAS_LIBRARIES})
target_link_libraries(xexample_DGELS_rowmajor ${LAPACKELIB} ${BLAS_LIBRARIES})
target_link_libraries(xexample_DGELS_colmajor ${LAPACKELIB} ${BLAS_LIBRARIES})

add_test(example_DGESV_rowmajor ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/xexample_DGESV_rowmajor)
add_test(example_DGESV_colmajor ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/xexample_DGESV_colmajor)
add_test(example_DGELS_rowmajor ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/xexample_DGELS_rowmajor)
add_test(example_DGELS_colmajor ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/xexample_DGELS_colmajor)

if(BUILD_INDEX64_EXT_API)
  add_executable(xexample_DGESV_rowmajor_64 example_DGESV_rowmajor_64.c lapacke_example_aux.c lapacke_example_aux.h)
  add_executable(xexample_DGESV_colmajor_64 example_DGESV_colmajor_64.c lapacke_example_aux.c lapacke_example_aux.h)
  add_executable(xexample_DGELS_rowmajor_64 example_DGELS_rowmajor_64.c lapacke_example_aux.c lapacke_example_aux.h)
  add_executable(xexample_DGELS_colmajor_64 example_DGELS_colmajor_64.c lapacke_example_aux.c lapacke_example_aux.h)

  target_link_libraries(xexample_DGESV_rowmajor_64 ${LAPACKELIB} ${BLAS_LIBRARIES})
  target_link_libraries(xexample_DGESV_colmajor_64 ${LAPACKELIB} ${BLAS_LIBRARIES})
  target_link_libraries(xexample_DGELS_rowmajor_64 ${LAPACKELIB} ${BLAS_LIBRARIES})
  target_link_libraries(xexample_DGELS_colmajor_64 ${LAPACKELIB} ${BLAS_LIBRARIES})

  add_test(example_DGESV_rowmajor_64 ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/xexample_DGESV_rowmajor_64)
  add_test(example_DGESV_colmajor_64 ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/xexample_DGESV_colmajor_64)
  add_test(example_DGELS_rowmajor_64 ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/xexample_DGELS_rowmajor_64)
  add_test(example_DGELS_colmajor_64 ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/xexample_DGELS_colmajor_64)
endif()
