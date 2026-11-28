# Compile/link only, so timer selection also works while cross compiling.
function(lapack_check_time_function result)
  foreach(timer INT_ETIME EXT_ETIME_ EXT_ETIME INT_CPU_TIME NONE)
    try_compile(works
      SOURCES "${LAPACK_SOURCE_DIR}/INSTALL/second_${timer}.f"
              "${LAPACK_SOURCE_DIR}/INSTALL/secondtst.f"
      NO_CACHE
      LOG_DESCRIPTION "LAPACK timer ${timer}")
    if(works)
      message(STATUS "LAPACK timer: ${timer}")
      set(${result} "${timer}" PARENT_SCOPE)
      return()
    endif()
  endforeach()
  message(FATAL_ERROR "Cannot compile any LAPACK timing implementation.")
endfunction()
