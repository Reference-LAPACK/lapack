# The LAPACKE tests use the unsuffixed API, which only exists when the
# default API is built (mirrors LAPACKE/example).
if(NOT BUILD_DEFAULT_API)
  return()
endif()

# Precisions to test: those whose LAPACKE interface is built. A
# precision-pruned build tests only what it has.
set(LAPACKE_TEST_PRECISIONS)
foreach(_spec "s;SINGLE" "d;DOUBLE" "c;COMPLEX" "z;COMPLEX16")
  list(GET _spec 0 _prec)
  list(GET _spec 1 _build)
  if(BUILD_${_build} AND LAPACKE_BUILD_${_build})
    list(APPEND LAPACKE_TEST_PRECISIONS ${_prec})
  endif()
endforeach()
if(NOT LAPACKE_TEST_PRECISIONS)
  return()
endif()

add_subdirectory(wrappers)
add_subdirectory(high_level)

# The LIN tests, rerouted through LAPACKE.
#
# The test sources are rewritten so that calls to allowlisted routines
# resolve to the <NAME>_TEST wrappers above, which forward to LAPACKE; the
# unmodified numerical and error-exit checks then exercise LAPACKE. String
# replacements must stay off so that the SRNAMT names used by the error-exit
# tests keep matching the names the unrenamed LAPACK routines pass to the
# (interposed) testing XERBLA.
include(ExtendedAPIHelpers)

# An empty SYMBOL_ALLOWLIST would rename every symbol, so pass the option
# only when there is something to route.
set(_allowlist_args)
if(LAPACKE_TEST_ALLOWLISTED_ROUTINES)
  list(APPEND _allowlist_args
    SYMBOL_ALLOWLIST ${LAPACKE_TEST_ALLOWLISTED_ROUTINES})
endif()

# Outputs live in the TESTING/lapacke subdirectory, where lapack_testing.py
# analyzes them as their own LAPACKE section. Only the work-level
# column-major flavor runs the error-exit tests: the row-major interface has
# different leading-dimension semantics (LDA >= N instead of LDA >= M), and
# the high-level interface ignores the caller's workspace so LWORK error
# exits cannot fire. The other flavors read a generated input with the
# error-exit tests disabled.
set(_lapacke_test_out_dir ${LAPACK_BINARY_DIR}/TESTING/lapacke)
file(MAKE_DIRECTORY ${_lapacke_test_out_dir})

foreach(_prec IN LISTS LAPACKE_TEST_PRECISIONS)
  string(TOUPPER "${_prec}" _prec_upper)
  if(NOT LAPACKE_TEST_LIN_${_prec_upper}_SOURCES)
    continue()
  endif()

  generate_suffixed_sources(xlintst${_prec}_lapacke
    LAPACKE_TEST_LIN_${_prec_upper}_SOURCES _lin_sources_test
    SUFFIX "_TEST"
    ${_allowlist_args}
    NO_STRING_REPLACEMENTS)

  add_library(xlintst${_prec}_lapacke_core OBJECT ${_lin_sources_test})

  file(READ ${LAPACK_SOURCE_DIR}/TESTING/${_prec}test.in _test_in_content)
  string(REPLACE
    "T                      Put T to test the error exits"
    "F                      Put T to test the error exits"
    _test_in_noerr_content "${_test_in_content}")
  file(WRITE ${_lapacke_test_out_dir}/${_prec}test_noerr.in "${_test_in_noerr_content}")

  foreach(_layer work high)
    foreach(_layout cm rm)
      set(_driver xlintst${_prec}_lapacke_${_layer}_${_layout})

      add_executable(${_driver} $<TARGET_OBJECTS:xlintst${_prec}_lapacke_core>)
      set_target_properties(${_driver} PROPERTIES LINKER_LANGUAGE Fortran)
      target_link_libraries(${_driver}
        lapacke_test_wrappers_${_layer}_${_layout}
        ${TMGLIB} ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})

      if(_layer STREQUAL "work" AND _layout STREQUAL "cm")
        set(_input ${LAPACK_SOURCE_DIR}/TESTING/${_prec}test.in)
      else()
        set(_input ${_lapacke_test_out_dir}/${_prec}test_noerr.in)
      endif()
      add_test(NAME LAPACK-${_driver}_${_prec}test_in
        COMMAND "${CMAKE_COMMAND}"
        -DTEST=$<TARGET_FILE:${_driver}>
        -DINPUT=${_input}
        -DOUTPUT=${_lapacke_test_out_dir}/${_prec}test_${_layer}_${_layout}.out
        -DINTDIR=${CMAKE_CFG_INTDIR}
        -P "${LAPACK_SOURCE_DIR}/TESTING/runtest.cmake")
    endforeach()
  endforeach()
endforeach()
