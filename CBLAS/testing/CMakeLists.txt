#######################################################################
#  This CMakeLists.txt creates the test programs for the CBLAS routines.
#
#######################################################################

function(add_cblas_ctest output input target)
  if(NOT "${input}" STREQUAL "")
    set(TEST_INPUT "${CMAKE_CURRENT_SOURCE_DIR}/${input}")
  endif()
  set(TEST_OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/${output}")
  set(testName "${target}")

  if(DEFINED TEST_INPUT AND EXISTS "${TEST_INPUT}")
    add_test(NAME CBLAS-${testName} COMMAND "${CMAKE_COMMAND}"
      -DTEST=$<TARGET_FILE:${target}>
      -DINPUT=${TEST_INPUT}
      -DOUTPUT=${TEST_OUTPUT}
      -DINTDIR=${CMAKE_CFG_INTDIR}
      -P "${LAPACK_SOURCE_DIR}/TESTING/runtest.cmake")
  else()
    add_test(NAME CBLAS-${testName} COMMAND "${CMAKE_COMMAND}"
      -DTEST=$<TARGET_FILE:${target}>
      -DOUTPUT=${TEST_OUTPUT}
      -DINTDIR=${CMAKE_CFG_INTDIR}
      -P "${LAPACK_SOURCE_DIR}/TESTING/runtest.cmake")
  endif()
endfunction()

function(add_cblas_test_target target)
  add_executable(${target} ${ARGN})

  if(HAS_ATTRIBUTE_WEAK_SUPPORT)
    target_compile_definitions(${target} PRIVATE HAS_ATTRIBUTE_WEAK_SUPPORT)
  endif()

  if(WIN32 AND BUILD_SHARED_LIBS)
    target_compile_definitions(${target} PRIVATE CBLAS_DLL_IMPORTS)
  endif()

  target_link_libraries(${target} PRIVATE ${CBLASLIB})
endfunction()

function(add_cblas_test target output input f_source)
  set(c_sources ${ARGN})

  if(BUILD_DEFAULT_API)
    add_cblas_test_target(${target} ${f_source} ${c_sources})
    add_cblas_ctest(${output} "${input}" ${target})
  endif()

  if(BUILD_INDEX64_EXT_API)
    include(ExtendedAPIHelpers)
    generate_64bit_suffixed_sources(${target} f_source f_source_64)

    add_cblas_test_target(${target}_64 ${f_source_64} ${c_sources})
    target_compile_definitions(${target}_64 PRIVATE
      "$<$<COMPILE_LANGUAGE:C>:WeirdNEC>" "$<$<COMPILE_LANGUAGE:C>:CBLAS_API64>")
    target_compile_options(${target}_64 PRIVATE
      "$<$<COMPILE_LANGUAGE:Fortran>:${FOPT_ILP64}>")

    get_filename_component(output_name "${output}" NAME_WLE)
    get_filename_component(output_ext "${output}" EXT)
    set(output_64 "${output_name}_64${output_ext}")
    add_cblas_ctest(${output_64} "${input}" ${target}_64)
  endif()
endfunction()

# Object files for single precision real
set(STESTL1O c_sblas1.c)
set(STESTL2O c_sblas2.c c_s2chke.c auxiliary.c c_xerbla.c)
set(STESTL3O c_sblas3.c c_s3chke.c auxiliary.c c_xerbla.c)

# Object files for double precision real
set(DTESTL1O c_dblas1.c)
set(DTESTL2O c_dblas2.c c_d2chke.c auxiliary.c c_xerbla.c)
set(DTESTL3O c_dblas3.c c_d3chke.c auxiliary.c c_xerbla.c)

# Object files for single precision complex
set(CTESTL1O c_cblas1.c)
set(CTESTL2O c_cblas2.c c_c2chke.c auxiliary.c c_xerbla.c)
set(CTESTL3O c_cblas3.c c_c3chke.c auxiliary.c c_xerbla.c)

# Object files for double precision complex
set(ZTESTL1O c_zblas1.c)
set(ZTESTL2O c_zblas2.c c_z2chke.c auxiliary.c c_xerbla.c)
set(ZTESTL3O c_zblas3.c c_z3chke.c auxiliary.c c_xerbla.c)

if(BUILD_SINGLE)
  add_cblas_test(xscblat1 stest1.out "" c_sblat1.f ${STESTL1O})
  add_cblas_test(xscblat2 stest2.out sin2 c_sblat2.f ${STESTL2O})
  add_cblas_test(xscblat3 stest3.out sin3 c_sblat3.f ${STESTL3O})
endif()

if(BUILD_DOUBLE)
  add_cblas_test(xdcblat1 dtest1.out "" c_dblat1.f ${DTESTL1O})
  add_cblas_test(xdcblat2 dtest2.out din2 c_dblat2.f ${DTESTL2O})
  add_cblas_test(xdcblat3 dtest3.out din3 c_dblat3.f ${DTESTL3O})
endif()

if(BUILD_COMPLEX)
  add_cblas_test(xccblat1 ctest1.out "" c_cblat1.f ${CTESTL1O})
  add_cblas_test(xccblat2 ctest2.out cin2 c_cblat2.f ${CTESTL2O})
  add_cblas_test(xccblat3 ctest3.out cin3 c_cblat3.f ${CTESTL3O})
endif()

if(BUILD_COMPLEX16)
  add_cblas_test(xzcblat1 ztest1.out "" c_zblat1.f ${ZTESTL1O})
  add_cblas_test(xzcblat2 ztest2.out zin2 c_zblat2.f ${ZTESTL2O})
  add_cblas_test(xzcblat3 ztest3.out zin3 c_zblat3.f ${ZTESTL3O})
endif()
