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
