function(_add_blas_test name src test_input test_output)
  add_executable(${name} ${src})
  target_link_libraries(${name} PRIVATE ${BLASLIB})
  lapack_add_coverage(${name})

  # A driver that reads an input file takes the name of its summary file from
  # it and opens that itself on unit 6, so redirecting stdout there as well
  # would open the one file twice; xBLAT2 fails outright with "Cannot change
  # STATUS parameter in OPEN statement".  Only xBLAT1, which has no input and
  # writes to stdout, needs the redirect.
  if(EXISTS "${test_input}")
    lapack_runtest_command(test_command ${name} INPUT "${test_input}")
  else()
    lapack_runtest_command(test_command ${name} OUTPUT "${test_output}")
  endif()
  add_test(NAME BLAS-${name} COMMAND ${test_command})

  # The special values (Inf, NaN) that SXVALS and DXVALS return do not
  # survive nagfor folding the expressions that make them.
  lapack_nag_disable_constant_propagation_target(${name})
endfunction()

function(add_blas_test name source)
  get_filename_component(baseNAME ${source} NAME_WE)
  lapack_test_input(test_input
    "${CMAKE_CURRENT_SOURCE_DIR}/${baseNAME}.in")

  if(BUILD_DEFAULT_API)
    _add_blas_test(${name} ${source} "${test_input}"
      "${CMAKE_CURRENT_BINARY_DIR}/${baseNAME}.out")
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

    _add_blas_test(${name}_64 "${source_64}" "${test_input_64}"
      "${CMAKE_CURRENT_BINARY_DIR}/${baseNAME}_64.out")
    target_compile_options(${name}_64 PRIVATE ${FOPT_ILP64})
  endif()
endfunction()

# These programs check their own results and return a nonzero status.
# Keep their stdout in CTest; classical .out files belong to lapack_testing.py.
function(add_blas_self_test name source)
  if(BUILD_DEFAULT_API)
    add_executable(${name} ${source})
    target_link_libraries(${name} PRIVATE ${BLASLIB})
    lapack_add_coverage(${name})
    lapack_runtest_command(test_command ${name})
    add_test(NAME BLAS-${name} COMMAND ${test_command})
  endif()
  if(BUILD_INDEX64_EXT_API)
    include(ExtendedAPIHelpers)
    generate_64bit_suffixed_sources(${name} source source_64 NO_STRING_REPLACEMENTS)
    add_executable(${name}_64 ${source_64})
    target_compile_options(${name}_64 PRIVATE ${FOPT_ILP64})
    target_link_libraries(${name}_64 PRIVATE ${BLASLIB})
    lapack_add_coverage(${name}_64)
    lapack_runtest_command(test_command ${name}_64)
    add_test(NAME BLAS-${name}_64 COMMAND ${test_command})
  endif()
endfunction()

if(BUILD_SINGLE)
  add_blas_self_test(xnrm2overflow_s test_snrm2_overflow.f90)
  add_blas_test(xblat1s sblat1.f)
  add_blas_test(xblat2s sblat2.f)
  add_blas_test(xblat3s sblat3.f)
endif()

if(BUILD_DOUBLE)
  add_blas_self_test(xnrm2overflow_d test_dnrm2_overflow.f90)
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
