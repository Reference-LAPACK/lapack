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

# LAPACKE test wrapper libraries.
#
# The wrappers map <routine>_TEST Fortran symbols onto LAPACKE calls. The
# LAPACK testing programs are rewritten by the extended-API source generator
# (see TESTING/LIN/CMakeLists.txt) so that calls to allowlisted routines
# resolve to these wrappers. The same wrapper sources are compiled once per
# (layer, layout) combination:
#   layer:  work -- middle-level (_work) interface, caller workspace passes
#                   through (faithful WORK/LWORK semantics, error exits)
#           high -- high-level interface (NaN checks, internal allocation)
#   layout: cm   -- column-major LAPACKE interface (faithful passthrough)
#           rm   -- row-major interface through row-major shadow copies of
#                   all matrix arguments

set(LAPACKE_TEST_WRAPPER_SOURCES wrappers/lapacke_test_wrappers.c)
foreach(_prec IN LISTS LAPACKE_TEST_PRECISIONS)
  list(APPEND LAPACKE_TEST_WRAPPER_SOURCES
    wrappers/lapacke_test_wrappers_${_prec}.c)
endforeach()

foreach(_layer work high)
  string(TOUPPER "${_layer}" _layer_upper)
  foreach(_layout cm rm)
    if(_layout STREQUAL "cm")
      set(_layout_define LAPACK_COL_MAJOR)
    else()
      set(_layout_define LAPACK_ROW_MAJOR)
    endif()

    set(_wrapper_lib lapacke_test_wrappers_${_layer}_${_layout})
    add_library(${_wrapper_lib} STATIC ${LAPACKE_TEST_WRAPPER_SOURCES})
    target_compile_definitions(${_wrapper_lib} PRIVATE
      LAPACKE_TEST_LAYOUT=${_layout_define}
      LAPACKE_TEST_LAYER=LAPACKE_TEST_LAYER_${_layer_upper})
    target_include_directories(${_wrapper_lib} PRIVATE
      ${LAPACK_SOURCE_DIR}/LAPACKE/include
      ${LAPACK_BINARY_DIR}/include)
    target_link_libraries(${_wrapper_lib} PUBLIC ${LAPACKELIB})
  endforeach()
endforeach()

# Dedicated tests for the high-level-only behavior the Fortran test suite
# cannot exercise: NaN checks actually detecting NaNs, and workspace /
# transposition allocation failures (via malloc failure injection).
#
# One test file per routine in high_level/<name>.c, each defining one
# function via LAPACKE_TEST(<name>). The runner main with the test table is
# generated from the list below; to add a test, add the file and list its
# name here.
set(LAPACKE_HIGH_LEVEL_ROUTINES
  matgen
  getrf
  getrs
  getri
  potrf
  potrs
  geqrf
  gels
  lange)

set(LAPACKE_HIGH_LEVEL_TESTS)
foreach(_prec IN LISTS LAPACKE_TEST_PRECISIONS)
  foreach(_routine IN LISTS LAPACKE_HIGH_LEVEL_ROUTINES)
    list(APPEND LAPACKE_HIGH_LEVEL_TESTS ${_prec}${_routine})
  endforeach()
endforeach()

set(LAPACKE_HIGH_LEVEL_SOURCES
  high_level/lapacke_test.c
  high_level/lapacke_test_alloc.c
  high_level/lapacke_test_nan_sweep.c)
foreach(_prec IN LISTS LAPACKE_TEST_PRECISIONS)
  list(APPEND LAPACKE_HIGH_LEVEL_SOURCES
    high_level/lapacke_test_matgen_${_prec}.c)
endforeach()
set(LAPACKE_HIGH_LEVEL_EXTERNS "")
set(LAPACKE_HIGH_LEVEL_ENTRIES "")
foreach(_test_name IN LISTS LAPACKE_HIGH_LEVEL_TESTS)
  list(APPEND LAPACKE_HIGH_LEVEL_SOURCES high_level/${_test_name}.c)
  string(APPEND LAPACKE_HIGH_LEVEL_EXTERNS
    "extern void lapacke_test_fn_${_test_name}( void );\n")
  string(APPEND LAPACKE_HIGH_LEVEL_ENTRIES
    "    { \"${_test_name}\", lapacke_test_fn_${_test_name} },\n")
endforeach()
configure_file(high_level/lapacke_high_level_main.c.in
  ${CMAKE_CURRENT_BINARY_DIR}/lapacke_high_level_main.c @ONLY)

add_executable(xlapacke_high_level
  ${CMAKE_CURRENT_BINARY_DIR}/lapacke_high_level_main.c
  ${LAPACKE_HIGH_LEVEL_SOURCES})
target_include_directories(xlapacke_high_level PRIVATE
  ${LAPACK_SOURCE_DIR}/LAPACKE/include
  ${LAPACK_BINARY_DIR}/include)
target_link_libraries(xlapacke_high_level
  ${LAPACKELIB} ${TMGLIB} ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
add_test(NAME LAPACKE-high-level COMMAND xlapacke_high_level)

# The LIN tests, rerouted through LAPACKE.
#
# The test sources are rewritten so that calls to allowlisted routines
# resolve to the <NAME>_TEST wrappers above, which forward to LAPACKE; the
# unmodified numerical and error-exit checks then exercise LAPACKE. String
# replacements must stay off so that the SRNAMT names used by the error-exit
# tests keep matching the names the unrenamed LAPACK routines pass to the
# (interposed) testing XERBLA.
include(ExtendedAPIHelpers)

# Generate the symbol allowlist: a routine is routed through LAPACKE if it
# has a LAPACKE implementation (a lapacke_<name>_work.c source) and a
# <NAME>_TEST wrapper in the LAPACKE test wrapper sources.
file(GLOB _lapacke_work_sources RELATIVE "${LAPACK_SOURCE_DIR}/LAPACKE/src"
  "${LAPACK_SOURCE_DIR}/LAPACKE/src/lapacke_*_work.c")
set(_lapacke_routines)
foreach(_work_source IN LISTS _lapacke_work_sources)
  string(REGEX REPLACE "^lapacke_(.+)_work\\.c$" "\\1" _routine "${_work_source}")
  list(APPEND _lapacke_routines "${_routine}")
endforeach()

set(_allowlisted_routines)
foreach(_prec IN LISTS LAPACKE_TEST_PRECISIONS)
  set(_wrapper_source
    "${CMAKE_CURRENT_SOURCE_DIR}/wrappers/lapacke_test_wrappers_${_prec}.c")
  # Rerun the configuration when wrappers are added or removed
  set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS "${_wrapper_source}")
  file(STRINGS "${_wrapper_source}" _wrapper_symbols
    REGEX "LAPACK_GLOBAL\\([a-z0-9_]+_test,")
  foreach(_wrapper_symbol IN LISTS _wrapper_symbols)
    string(REGEX MATCH "LAPACK_GLOBAL\\(([a-z0-9_]+)_test," _unused "${_wrapper_symbol}")
    set(_routine "${CMAKE_MATCH_1}")
    if(_routine IN_LIST _lapacke_routines)
      list(APPEND _allowlisted_routines "${_routine}")
    else()
      message(WARNING
        "LAPACKE test wrapper '${_routine}_test' has no matching LAPACKE"
        " routine and is not routed through LAPACKE")
    endif()
  endforeach()
endforeach()
list(REMOVE_DUPLICATES _allowlisted_routines)
list(SORT _allowlisted_routines)
list(LENGTH _allowlisted_routines _allowlisted_routine_count)
message(STATUS "LAPACKE testing: routing ${_allowlisted_routine_count} routines through LAPACKE")

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
    SYMBOL_ALLOWLIST ${_allowlisted_routines}
    NO_STRING_REPLACEMENTS)

  add_library(xlintst${_prec}_lapacke_core OBJECT ${_lin_sources_test})

  file(READ ${LAPACK_SOURCE_DIR}/TESTING/${_prec}test.in _test_in_content)
  string(REPLACE
    "T                      Put T to test the error exits"
    "F                      Put T to test the error exits"
    _test_in_noerr_content "${_test_in_content}")
  file(WRITE ${_lapacke_test_out_dir}/${_prec}test_noerr.in
    "${_test_in_noerr_content}")

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
