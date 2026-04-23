#######################################################################
#  This CMakeLists.txt creates the test programs for the CBLAS routines.
#
#######################################################################

macro(add_cblas_test output input target)
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
endmacro()

macro(configure_cblas_test_target target test_sources)
  add_executable(${target} ${test_sources})

  if(HAS_ATTRIBUTE_WEAK_SUPPORT)
    target_compile_definitions(${target} PRIVATE HAS_ATTRIBUTE_WEAK_SUPPORT)
  endif()

  if(WIN32 AND BUILD_SHARED_LIBS)
    target_compile_definitions(${target} PRIVATE CBLAS_DLL_IMPORTS)
  endif()

  target_link_libraries(${target} PRIVATE ${CBLASLIB})
endmacro()

# Object files for single precision real
set(STESTL1O c_sblat1.f c_sblas1.c)
set(STESTL2O c_sblat2.f c_sblas2.c c_s2chke.c auxiliary.c c_xerbla.c)
set(STESTL3O c_sblat3.f c_sblas3.c c_s3chke.c auxiliary.c c_xerbla.c)

# Object files for double precision real
set(DTESTL1O c_dblat1.f c_dblas1.c)
set(DTESTL2O c_dblat2.f c_dblas2.c c_d2chke.c auxiliary.c c_xerbla.c)
set(DTESTL3O c_dblat3.f c_dblas3.c c_d3chke.c auxiliary.c c_xerbla.c)

# Object files for single precision complex
set(CTESTL1O c_cblat1.f c_cblas1.c)
set(CTESTL2O c_cblat2.f c_cblas2.c c_c2chke.c auxiliary.c c_xerbla.c)
set(CTESTL3O c_cblat3.f c_cblas3.c c_c3chke.c auxiliary.c c_xerbla.c)

# Object files for double precision complex
set(ZTESTL1O c_zblat1.f c_zblas1.c)
set(ZTESTL2O c_zblat2.f c_zblas2.c c_z2chke.c auxiliary.c c_xerbla.c)
set(ZTESTL3O c_zblat3.f c_zblas3.c c_z3chke.c auxiliary.c c_xerbla.c)

if(BUILD_SINGLE)
  configure_cblas_test_target(xscblat1 "${STESTL1O}")
  configure_cblas_test_target(xscblat2 "${STESTL2O}")
  configure_cblas_test_target(xscblat3 "${STESTL3O}")

  add_cblas_test(stest1.out ""   xscblat1)
  add_cblas_test(stest2.out sin2 xscblat2)
  add_cblas_test(stest3.out sin3 xscblat3)
endif()

if(BUILD_DOUBLE)
  configure_cblas_test_target(xdcblat1 "${DTESTL1O}")
  configure_cblas_test_target(xdcblat2 "${DTESTL2O}")
  configure_cblas_test_target(xdcblat3 "${DTESTL3O}")

  add_cblas_test(dtest1.out ""   xdcblat1)
  add_cblas_test(dtest2.out din2 xdcblat2)
  add_cblas_test(dtest3.out din3 xdcblat3)
endif()

if(BUILD_COMPLEX)
  configure_cblas_test_target(xccblat1 "${CTESTL1O}")
  configure_cblas_test_target(xccblat2 "${CTESTL2O}")
  configure_cblas_test_target(xccblat3 "${CTESTL3O}")

  add_cblas_test(ctest1.out ""   xccblat1)
  add_cblas_test(ctest2.out cin2 xccblat2)
  add_cblas_test(ctest3.out cin3 xccblat3)
endif()

if(BUILD_COMPLEX16)
  configure_cblas_test_target(xzcblat1 "${ZTESTL1O}")
  configure_cblas_test_target(xzcblat2 "${ZTESTL2O}")
  configure_cblas_test_target(xzcblat3 "${ZTESTL3O}")

  add_cblas_test(ztest1.out ""   xzcblat1)
  add_cblas_test(ztest2.out zin2 xzcblat2)
  add_cblas_test(ztest3.out zin3 xzcblat3)
endif()
