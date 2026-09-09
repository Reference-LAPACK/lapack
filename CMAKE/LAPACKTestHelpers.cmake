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

# Set ${out_var} to the test input to feed to a test driver.  Almost every
# driver reads a TSTERR flag from its input, so where the error-exit tests
# cannot run this is a copy in the build tree with that flag turned off;
# everywhere else, and for the drivers that take no input, it is ${input}
# itself.  The drivers that do not read the flag are handled in the source,
# see LAPACK_SKIP_ERROR_EXIT_TESTS in TESTING/EIG/xchkee.F.
function(lapack_test_input out_var input)
  if(NOT LAPACK_SKIP_ERROR_EXIT_TESTS OR NOT EXISTS "${input}")
    set(${out_var} "${input}" PARENT_SCOPE)
    return()
  endif()

  file(READ "${input}" content)
  string(REGEX REPLACE
    "\n[ \t]*(T|\\.TRUE\\.)([ \t]+[^\n]*(Put T to test the error exits|LOGICAL FLAG, T TO TEST ERROR EXITS))"
    "\nF\\2" content "${content}")

  get_filename_component(name "${input}" NAME)
  set(rewritten "${CMAKE_CURRENT_BINARY_DIR}/${name}")
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
    if(ARG_INPUT)
      string(APPEND command " < \"${ARG_INPUT}\"")
    endif()
    if(ARG_OUTPUT)
      string(APPEND command " > \"${ARG_OUTPUT}\" 2> \"${ARG_OUTPUT}.err\"")
    endif()
    set(${out_var} sh -c "${command}" PARENT_SCOPE)
  else()
    set(command "${CMAKE_COMMAND}" -DTEST=$<TARGET_FILE:${target}>)
    if(ARG_INPUT)
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
