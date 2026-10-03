# The test drivers ship their own XERBLA, which records the expected INFO and
# returns, and rely on it replacing the one in the library, which prints and
# stops.  A shared library on Windows always calls the XERBLA it was linked
# with, so the first deliberate illegal argument stops the driver before any
# test runs.  CBLAS skips its xerbla tests for the same reason, see
# CBLAS/testing.  macOS binds a shared library with a two-level namespace by
# default, which keeps the override from taking effect in the same way;
# USE_FLAT_NAMESPACE restores it, see the top-level CMakeLists.txt.
#
# Recomputed on every configure, so that toggling BUILD_SHARED_LIBS in an
# existing build tree is picked up, but left alone when it was asked for on
# the command line: -D puts it in the cache, which defines it here.
if(NOT DEFINED LAPACK_SKIP_ERROR_EXIT_TESTS)
  if(WIN32 AND BUILD_SHARED_LIBS)
    set(LAPACK_SKIP_ERROR_EXIT_TESTS ON)
    message(STATUS
      "Disabling the error-exit tests: the test suite's XERBLA cannot replace "
      "the one in a shared library on Windows")
  elseif(APPLE AND BUILD_SHARED_LIBS AND NOT USE_FLAT_NAMESPACE)
    set(LAPACK_SKIP_ERROR_EXIT_TESTS ON)
    message(STATUS
      "Disabling the error-exit tests: the test suite's XERBLA cannot replace "
      "the one in a shared library bound with the default two-level namespace; "
      "configure with -D USE_FLAT_NAMESPACE=ON to run them")
  else()
    set(LAPACK_SKIP_ERROR_EXIT_TESTS OFF)
  endif()
endif()

# nagfor folds expressions that deliberately manufacture an Inf or a NaN
# and then rejects the result ("Invalid operand for intrinsic SQRT ... Errors
# found during constant propagation"), so the test drivers that do that are
# built with constant propagation disabled.  The compiler id and the flag it
# takes live here only.
#
# Read when this is called rather than when the file is included: the top
# level enables Fortran well after that, so CMAKE_Fortran_COMPILER_ID is
# still empty at include time.
function(lapack_no_constant_propagation_flag out_var)
  if(CMAKE_Fortran_COMPILER_ID STREQUAL "NAG")
    set(${out_var} "-Onopropagate" PARENT_SCOPE)
  else()
    set(${out_var} "" PARENT_SCOPE)
  endif()
endfunction()

# Disable it for the CXX error-exit and packed symmetric-indefinite tests
# among ARGN, which manufacture NaNs. Call this with every source list:
# the generated _64 and _TEST copies need it as much as the originals.
function(lapack_nag_disable_constant_propagation)
  lapack_no_constant_propagation_flag(flag)
  if(NOT flag)
    return()
  endif()
  foreach(source IN LISTS ARGN)
    get_filename_component(name "${source}" NAME)
    if(name MATCHES "^[scdz](errcxx|chksp|chkhp)(_[A-Za-z0-9]+)*\\.f$")
      set_source_files_properties("${source}"
        PROPERTIES COMPILE_OPTIONS "${flag}")
    endif()
  endforeach()
endfunction()

# The same for a whole target: the BLAS drivers take the special values they
# test from SXVALS and DXVALS, so there is no one source file to pick out by
# name here.
function(lapack_nag_disable_constant_propagation_target target)
  lapack_no_constant_propagation_flag(flag)
  if(NOT flag)
    return()
  endif()
  target_compile_options(${target} PRIVATE "${flag}")
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

  set(rewritten "${ARG_OUTPUT}")
  if(NOT rewritten)
    get_filename_component(name "${input}" NAME)
    set(rewritten "${CMAKE_CURRENT_BINARY_DIR}/${name}")
  endif()

  set(content "")
  set(disabled "")
  if(EXISTS "${input}" AND
      (LAPACK_SKIP_ERROR_EXIT_TESTS OR ARG_DISABLE_ERROR_EXIT_TESTS))
    # The rewrite happens at configure time, so ask CMake to re-run when
    # ${input} changes; otherwise the build tree goes on serving a copy
    # made from a version of it that no longer exists.
    set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS "${input}")
    file(READ "${input}" content)
    string(REGEX REPLACE
      "\n[ \t]*T([ \t]+[^\n]*(Put T to test the error exits|LOGICAL FLAG, T TO TEST ERROR EXITS))"
      "\nF\\1" disabled "${content}")
  endif()

  # Only a copy that differs is worth having.  Everywhere else -- a build
  # that runs the error exits, an input with no flag to turn off -- uses
  # ${input} itself, and any copy an earlier configure left behind goes:
  # lapack_testing.py --run prefers the copy in the build tree, so a stale
  # one would feed the drivers what ctest never sees.
  if(disabled STREQUAL content)
    file(REMOVE "${rewritten}")
    set(${out_var} "${input}" PARENT_SCOPE)
  else()
    file(WRITE "${rewritten}" "${disabled}")
    set(${out_var} "${rewritten}" PARENT_SCOPE)
  endif()
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
