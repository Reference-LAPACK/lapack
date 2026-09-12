# The test drivers ship their own XERBLA, which records the expected INFO and
# returns, and rely on it replacing the one in the library, which prints and
# stops.  A shared library on Windows always calls the XERBLA it was linked
# with, so the first deliberate illegal argument stops the driver before any
# test runs.  CBLAS skips its xerbla tests for the same reason, see
# CBLAS/testing.
if(WIN32 AND BUILD_SHARED_LIBS)
  set(LAPACK_SKIP_ERROR_EXIT_TESTS ON)
  message(STATUS
    "Disabling the error-exit tests: the test suite's XERBLA cannot replace "
    "the one in a shared library on Windows")
else()
  set(LAPACK_SKIP_ERROR_EXIT_TESTS OFF)
endif()

# nagfor folds the SQRT( -ONE ) with which the ?errcxx drivers manufacture
# a NaN and rejects it ("Invalid operand for intrinsic SQRT ... Errors found
# during constant propagation"), so those drivers are compiled with constant
# propagation disabled.  Call this with every source list such a driver is
# built from: the generated _64 and _TEST copies need it as much as the
# originals.
function(lapack_nag_disable_constant_propagation)
  if(NOT CMAKE_Fortran_COMPILER_ID STREQUAL "NAG")
    return()
  endif()
  foreach(source IN LISTS ARGN)
    get_filename_component(name "${source}" NAME)
    if(name MATCHES "^[scdz]errcxx(_[A-Za-z0-9]+)*\\.f$")
      set_source_files_properties("${source}"
        PROPERTIES COMPILE_OPTIONS "-Onopropagate")
    endif()
  endforeach()
endfunction()

# Set ${out_var} to the test input to feed to a test driver.  Every driver
# that tests error exits reads a TSTERR flag from its input, so where those
# tests cannot run, or where DISABLE_ERROR_EXIT_TESTS asks for it anyway,
# this is a copy with the flag turned off; everywhere else, and for the
# drivers that take no input, it is ${input} itself.
#
# The copy goes to OUTPUT, or to the name of ${input} in the current binary
# directory.  Write it where lapack_testing.py looks for that input, which is
# the directory holding the .out files of the tests that read it, so that
# running the drivers from the script uses the same input ctest does.
function(lapack_test_input out_var input)
  cmake_parse_arguments(ARG "DISABLE_ERROR_EXIT_TESTS" "OUTPUT" "" ${ARGN})

  if(NOT EXISTS "${input}" OR
      NOT (LAPACK_SKIP_ERROR_EXIT_TESTS OR ARG_DISABLE_ERROR_EXIT_TESTS))
    set(${out_var} "${input}" PARENT_SCOPE)
    return()
  endif()

  file(READ "${input}" content)
  string(REGEX REPLACE
    "\n[ \t]*(T|\\.TRUE\\.)([ \t]+[^\n]*(Put T to test the error exits|LOGICAL FLAG, T TO TEST ERROR EXITS))"
    "\nF\\2" content "${content}")

  set(rewritten "${ARG_OUTPUT}")
  if(NOT rewritten)
    get_filename_component(name "${input}" NAME)
    set(rewritten "${CMAKE_CURRENT_BINARY_DIR}/${name}")
  endif()
  file(WRITE "${rewritten}" "${content}")
  set(${out_var} "${rewritten}" PARENT_SCOPE)
endfunction()

# Build the ctest COMMAND that runs ${target}, optionally reading standard
# input from INPUT and writing standard output to OUTPUT.
#
# Under LAPACK_MEMORY_CHECK on Unix the command execs the executable, so the
# test runs as a single process.  The cmake wrapper used otherwise starts the
# executable with execute_process(), which leaves valgrind with two processes
# writing to the one --log-file it was given, corrupting their records; exec
# replaces the shell instead of forking, so only one process ever writes.
#
# Everything else keeps the cmake wrapper.  It needs no shell, it handles the
# per-configuration directory of multi-config generators, and it echoes the
# output of a failing test, which a bare exec cannot do.
function(lapack_runtest_command out_var target)
  cmake_parse_arguments(ARG "" "INPUT;OUTPUT" "" ${ARGN})

  if(UNIX AND LAPACK_MEMORY_CHECK)
    set(command "exec \"$<TARGET_FILE:${target}>\"")
    if(ARG_INPUT AND EXISTS "${ARG_INPUT}")
      string(APPEND command " < \"${ARG_INPUT}\"")
    endif()
    if(ARG_OUTPUT)
      string(APPEND command " > \"${ARG_OUTPUT}\" 2> \"${ARG_OUTPUT}.err\"")
    endif()
    set(${out_var} sh -c "${command}" PARENT_SCOPE)
  else()
    set(command "${CMAKE_COMMAND}" -DTEST=$<TARGET_FILE:${target}>)
    if(ARG_INPUT AND EXISTS "${ARG_INPUT}")
      list(APPEND command -DINPUT=${ARG_INPUT})
    endif()
    if(ARG_OUTPUT)
      list(APPEND command -DOUTPUT=${ARG_OUTPUT})
    endif()
    list(APPEND command -DINTDIR=${CMAKE_CFG_INTDIR}
      -P "${LAPACK_SOURCE_DIR}/TESTING/runtest.cmake")
    set(${out_var} ${command} PARENT_SCOPE)
  endif()
endfunction()
