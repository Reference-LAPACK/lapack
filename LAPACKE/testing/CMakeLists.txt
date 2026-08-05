# The LAPACKE tests use the unsuffixed API, which only exists when the
# default API is built (mirrors LAPACKE/example), and cover the double
# precision routines, so precision-pruned LAPACKE builds must skip them
# (the LAPACKE_d* symbols and the DLATMS matrix generator would be missing).
if(NOT BUILD_DEFAULT_API OR NOT LAPACKE_BUILD_DOUBLE)
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

set(LAPACKE_TEST_WRAPPER_SOURCES
  wrappers/lapacke_test_wrappers.c
  wrappers/lapacke_test_wrappers_d.c)

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
set(LAPACKE_HIGH_LEVEL_TESTS
  matgen
  dgetrf
  dgetrs
  dgetri
  dpotrf
  dpotrs
  dgeqrf
  dgels
  dlange)

set(LAPACKE_HIGH_LEVEL_SOURCES
  high_level/lapacke_test.c
  high_level/lapacke_test_alloc.c
  high_level/lapacke_test_matgen.c
  high_level/lapacke_test_nan_sweep.c)
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
