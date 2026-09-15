# Dedicated tests for the high-level-only behavior the Fortran test suite
# cannot exercise: NaN checks actually detecting NaNs, and workspace /
# transposition allocation failures (via malloc failure injection).
#
# One test file per routine and precision in <x><name>.c, each defining one
# function via LAPACKE_TEST(<x><name>). There is one runner per precision,
# whose test table is generated from the three lists below: the routines
# LAPACKE provides in all four precisions, in the real ones only and in the
# complex ones only. To add a test, add the file for every precision the
# routine exists in and list its name in the matching list. The complex
# runners also link the real generator of their precision, for the real
# vectors the complex routines take.
set(LAPACKE_HIGH_LEVEL_ROUTINES
  bdsqr gbcon gbequ gbrfs gbsvx gbtrs gecon gecxx geequ gelq2 gelqf gels
  gelsd gelss gelsy geqlf geqp3 geqr geqr2 geqrf geqrfp geqrt geqrt2 geqrt3
  gerfs gerqf gesv gesvd gesvx getf2 getrf getri getrs getsls getsqrhrt
  gtcon gtrfs gtsv gtsvx gttrf gttrs lange lansy lantr lascl laset laswp
  matgen pbcon pbequ pbrfs pbsv pbsvx pbtrf pbtrs pocon poequ porfs posv
  posvx potrf potri potrs ppcon ppequ pprfs ppsv ppsvx pptrf pptri pptrs
  pstrf ptcon ptrfs ptsv ptsvx pttrf pttrs spcon sprfs spsv spsvx sptrf
  sptri sptrs sycon sycon_3 syrfs sysv sysv_rk sysv_rook sysvx sytrf
  sytrf_rk sytrf_rook sytri sytri2 sytri_3 sytrs sytrs2 sytrs_rook tbcon
  tbrfs tbtrs trcon trrfs trtri trtrs tzrzf)
set(LAPACKE_HIGH_LEVEL_ROUTINES_REAL
  orglq orgql orgqr orgrq orhr_col ormlq ormql ormqr sysv_aa sytrf_aa)
set(LAPACKE_HIGH_LEVEL_ROUTINES_COMPLEX
  hecon hecon_3 herfs hesv hesv_aa hesv_rk hesvx hetrf hetrf_aa hetrf_rk
  hetrf_rook hetri hetri2 hetri_3 hetrs hetrs2 hetrs_rook hpcon hprfs hpsv
  hpsvx hptrf hptri hptrs lacgv lanhe rot unglq ungql ungqr ungrq unhr_col
  unmlq unmql unmqr)

set(LAPACKE_HIGH_LEVEL_COMMON_SOURCES
  lapacke_test.c lapacke_test_alloc.c lapacke_test_nan_sweep.c)

foreach(_prec IN LISTS LAPACKE_TEST_PRECISIONS)
  set(_sources ${LAPACKE_HIGH_LEVEL_COMMON_SOURCES} lapacke_test_matgen_${_prec}.c)
  if(_prec STREQUAL "c")
    list(APPEND _sources lapacke_test_matgen_s.c)
  elseif(_prec STREQUAL "z")
    list(APPEND _sources lapacke_test_matgen_d.c)
  endif()

  # Test table for this precision, substituted into the runner main below.
  set(_routines ${LAPACKE_HIGH_LEVEL_ROUTINES})
  if(_prec STREQUAL "s" OR _prec STREQUAL "d")
    list(APPEND _routines ${LAPACKE_HIGH_LEVEL_ROUTINES_REAL})
  else()
    list(APPEND _routines ${LAPACKE_HIGH_LEVEL_ROUTINES_COMPLEX})
  endif()
  list(SORT _routines)

  set(LAPACKE_HIGH_LEVEL_EXTERNS "")
  set(LAPACKE_HIGH_LEVEL_ENTRIES "")
  foreach(_routine IN LISTS _routines)
    set(_test_name ${_prec}${_routine})
    list(APPEND _sources ${_test_name}.c)
    string(APPEND LAPACKE_HIGH_LEVEL_EXTERNS
      "extern void lapacke_test_fn_${_test_name}( void );\n")
    string(APPEND LAPACKE_HIGH_LEVEL_ENTRIES
      "    { \"${_test_name}\", lapacke_test_fn_${_test_name} },\n")
  endforeach()

  if(BUILD_DEFAULT_API)
    set(_runner ${_prec}lapacke_high_level)
    set(LAPACKE_HIGH_LEVEL_PROGRAM ${_runner})
    configure_file(lapacke_high_level_main.c.in
      ${CMAKE_CURRENT_BINARY_DIR}/${_runner}_main.c @ONLY)

    add_executable(${_runner} ${CMAKE_CURRENT_BINARY_DIR}/${_runner}_main.c ${_sources})
    target_include_directories(${_runner} PRIVATE
      ${LAPACK_BINARY_DIR}/include ${LAPACK_SOURCE_DIR}/LAPACKE/include)
    target_link_libraries(${_runner}
      ${LAPACKELIB} ${TMGLIB} ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})

    add_test(NAME LAPACKE-high-level-${_prec} COMMAND ${_runner})
  endif()

  if(BUILD_INDEX64_EXT_API)
    set(_runner ${_prec}lapacke_high_level_64)
    set(LAPACKE_HIGH_LEVEL_PROGRAM ${_runner})
    configure_file(lapacke_high_level_main.c.in
      ${CMAKE_CURRENT_BINARY_DIR}/${_runner}_main.c @ONLY)

    add_executable(${_runner} ${CMAKE_CURRENT_BINARY_DIR}/${_runner}_main.c ${_sources})
    target_compile_definitions(${_runner} PRIVATE
      LAPACK_ILP64 LAPACKE_API64 WeirdNEC CBLAS_API64)
    target_include_directories(${_runner} PRIVATE
      ${LAPACK_BINARY_DIR}/include ${LAPACK_SOURCE_DIR}/LAPACKE/include)
    target_link_libraries(${_runner}
      ${LAPACKELIB} ${TMGLIB} ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})

    add_test(NAME LAPACKE-high-level-${_prec}_64 COMMAND ${_runner})
  endif()
endforeach()
