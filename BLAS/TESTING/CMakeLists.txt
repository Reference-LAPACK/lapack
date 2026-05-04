function(_add_blas_test name src test_input)
  add_executable(${name} ${src})
  target_link_libraries(${name} ${BLASLIB})
  if(EXISTS "${test_input}")
    add_test(NAME BLAS-${name} COMMAND "${CMAKE_COMMAND}"
      -DTEST=$<TARGET_FILE:${name}>
      -DINPUT=${test_input}
      -DINTDIR=${CMAKE_CFG_INTDIR}
      -P "${LAPACK_SOURCE_DIR}/TESTING/runtest.cmake")
  else()
    add_test(NAME BLAS-${name} COMMAND "${CMAKE_COMMAND}"
      -DTEST=$<TARGET_FILE:${name}>
      -DINTDIR=${CMAKE_CFG_INTDIR}
      -P "${LAPACK_SOURCE_DIR}/TESTING/runtest.cmake")
  endif()

  # Disable constant propagation for NAG compiler to avoid issues with
  # special values (Inf, NaN) returned by SXVALS and DXVALS.
  if(CMAKE_Fortran_COMPILER_ID STREQUAL "NAG")
    target_compile_options(${name} PRIVATE "-Onopropagate")
  endif()
endfunction()

function(add_blas_test name source)
  get_filename_component(baseNAME ${source} NAME_WE)
  set(test_input "${CMAKE_CURRENT_SOURCE_DIR}/${baseNAME}.in")

  if(BUILD_DEFAULT_API)
    _add_blas_test(${name} ${source} "${test_input}")
  endif()

  if(BUILD_INDEX64_EXT_API)
    include(ExtendedAPIHelpers)
    generate_64bit_suffixed_sources(${name} source source_64 NO_STRING_REPLACEMENTS)

    # Create 64-bit version of test input if it exists, replacing output
    # file names with their 64-bit counterparts
    if(EXISTS "${test_input}")
      file(READ "${test_input}" test_input_content)
      string(REPLACE "'${baseNAME}.out'" "'${baseNAME}_64.out'"
        test_input_content "${test_input_content}")
      set(test_input_64 "${CMAKE_CURRENT_BINARY_DIR}/${baseNAME}_64.in")
      file(WRITE "${test_input_64}" "${test_input_content}")
    endif()

    _add_blas_test(${name}_64 "${source_64}" "${test_input_64}")
    target_compile_options(${name}_64 PRIVATE ${FOPT_ILP64})
  endif()
endfunction()

if(BUILD_SINGLE)
  add_blas_test(xblat1s sblat1.f)
  add_blas_test(xblat2s sblat2.f)
  add_blas_test(xblat3s sblat3.f)
endif()

if(BUILD_DOUBLE)
  add_blas_test(xblat1d dblat1.f)
  add_blas_test(xblat2d dblat2.f)
  add_blas_test(xblat3d dblat3.f)
endif()

if(BUILD_COMPLEX)
  add_blas_test(xblat1c cblat1.f)
  add_blas_test(xblat2c cblat2.f)
  add_blas_test(xblat3c cblat3.f)
endif()

if(BUILD_COMPLEX16)
  add_blas_test(xblat1z zblat1.f)
  add_blas_test(xblat2z zblat2.f)
  add_blas_test(xblat3z zblat3.f)
endif()
