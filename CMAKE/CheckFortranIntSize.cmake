# This module perdorms several try-compiles to determine the default integer
# size being used by the fortran compiler
#
# After execution, the following variable is set:
# 
# SIZEOF_INTEGER - Number of bytes used to store the default integer type
#                  Will be set to 1, 2, 4, 8 or 16 if successful, otherwise it
#                  will be unset
#  
#=============================================================================
# Author: Chuck Atkins
# Copyright 2010
#=============================================================================

macro( CHECK_FORTRAN_INT_SIZE )
  if( NOT CMAKE_Fortran_COMPILER_SUPPORTS_F90 )
    message( FATAL_ERROR "Type size tests require Fortran 90 support" )
  endif()
  if( NOT SIZEOF_INTEGER )
    foreach( _TEST_SIZE 1 2 4 8 16 )
      message( STATUS "Testing default integer*${_TEST_SIZE} - " )
      set( _TEST_FILE ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/testFortranInteger${_TEST_SIZE}.f90 )
      file( WRITE ${_TEST_FILE} "
       module check_int_size
         interface int_size
           function int_size_${_TEST_SIZE}(a) result(int_size)
           integer*${_TEST_SIZE} a
             integer int_size
             end function int_size_${_TEST_SIZE}
         end interface
       end module

       program check_int_size_${_TEST_SIZE}
       use check_int_size
         integer a, b
         a = 7
         b = int_size(a)
       end program 
     
       function int_size_${_TEST_SIZE}(a) result(int_size)
         integer*4 a
         int_size = 7
       end function
      ")
      try_compile( SIZEOF_INTEGER ${CMAKE_BINARY_DIR} ${_TEST_FILE} )
      if( SIZEOF_INTEGER )
        message( STATUS "Testing default integer*${_TEST_SIZE} - found" )
        set( SIZEOF_INTEGER ${_TEST_SIZE} CACHE INTERNAL "Size of the default INTEGER type" FORCE )
        break()
      else()
        message( STATUS "Testing default integer*${_TEST_SIZE} - not found" )
      endif()
    endforeach()
  endif()
  
  if( NOT SIZEOF_INTEGER )
    unset( SIZEOF_INTEGER )
  endif()
endmacro()

