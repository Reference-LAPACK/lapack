# This module perdorms several try-compiles to determine the default integer
# size being used by the fortran compiler
#
# After execution, the following variables are set.  If they are un set then
# size detection was not possible
# 
# SIZEOF_CHARACTER - Number of bytes used to store the default CHARACTER type
# SIZEOF_LOGICAL   - Number of bytes used to store the default LOGICAL type
# SIZEOF_INTEGER   - Number of bytes used to store the default INTEGER type
# SIZEOF_REAL      - Number of bytes used to store the default REAL type
# SIZEOF_COMPLEX   - Number of bytes used to store the default COMPLEX type
#  
#=============================================================================
# Author: Chuck Atkins
# Copyright 2010
#=============================================================================

macro( CHECK_FORTRAN_TYPE_SIZES )
  if( NOT CMAKE_Fortran_COMPILER_SUPPORTS_F90 )
    message( FATAL_ERROR "Type size tests require Fortran 90 support" )
  endif()
  foreach( _TEST_TYPE "CHARACTER" "LOGICAL" "INTEGER" "REAL" "COMPLEX" )
    string( REPLACE " " "_" _TEST_TYPE_VAR "${_TEST_TYPE}" )
    foreach( _TEST_SIZE 1 2 4 8 16 32 )
      set( _TEST_FILE ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/testFortran${_TEST_TYPE_VAR}Size${_TEST_SIZE}.f90 )
      file( WRITE ${_TEST_FILE} 
"
       PROGRAM check_size
         ${_TEST_TYPE}*${_TEST_SIZE}, TARGET :: a
         ${_TEST_TYPE}, POINTER :: pa
         pa => a
       END PROGRAM
")
      try_compile( SIZEOF_${_TEST_TYPE_VAR} ${CMAKE_BINARY_DIR} ${_TEST_FILE} )
      if( SIZEOF_${_TEST_TYPE_VAR} )
        message( STATUS "Testing default ${_TEST_TYPE}*${_TEST_SIZE} - found" )
        set( SIZEOF_${_TEST_TYPE_VAR} ${_TEST_SIZE} CACHE INTERNAL "Size of the default ${_TEST_TYPE} type" FORCE )
        break()
      endif()
    endforeach()
    if( NOT SIZEOF_${_TEST_TYPE_VAR} )
      message( WARNING "Unable to determine default size of type ${_TEST_TYPE}" )
      unset( SIZEOF_${_TEST_TYPE_VAR} CACHE )
    endif()
  endforeach()
endmacro()

