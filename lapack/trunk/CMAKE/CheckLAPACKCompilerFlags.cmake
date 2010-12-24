# This module checks against various known compilers and thier respective
# flags to determine any specific flags needing to be set.
# 
#  1.  If FPE traps are enabled either abort or disable them
#  2.  Specify fixed form if needed
#
#=============================================================================
# Author: Chuck Atkins
# Copyright 2010 Kitware, Inc.
#
# Distributed under the OSI-approved BSD License (the "License");
# see accompanying file Copyright.txt for details.
#
# This software is distributed WITHOUT ANY WARRANTY; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the License for more information.
#=============================================================================

macro( CheckLAPACKCompilerFlags )

set( FPE_EXIT FALSE )

# GNU Fortran
if( CMAKE_Fortran_COMPILER_ID STREQUAL "GNU" )
  if( "${CMAKE_Fortran_FLAGS}" MATCHES "-ffpe-trap=[izoupd]") 
    set( FPE_EXIT TRUE )
  endif()

# Intel Fortran
elseif( CMAKE_Fortran_COMPILER_ID STREQUAL "Intel" )
  if( ("${CMAKE_Fortran_FLAGS}" MATCHES "[-/]fpe0") OR
      ("${CMAKE_Fortran_FLAGS}" MATCHES "[-/]fpe-all=0") )
    set( FPE_EXIT TRUE )
  endif()

# SunPro F95
elseif( CMAKE_Fortran_COMPILER_ID STREQUAL "SunPro" )
  if( ("${CMAKE_Fortran_FLAGS}" MATCHES "-ftrap=") AND
      NOT ("${CMAKE_Fortran_FLAGS}" MATCHES "-ftrap=%none") )
    set( FPE_EXIT TRUE )
  elseif( NOT (CMAKE_Fortran_FLAGS MATCHES "-ftrap=") )
    message( STATUS "Disabling FPE trap handlers with -ftrap=%none" )
    set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ftrap=%none"
         CACHE STRING "Flags for Fortran compiler." FORCE )
  endif()

# IBM XL Fortran
elseif( (CMAKE_Fortran_COMPILER_ID STREQUAL "VisualAge" ) OR  # CMake 2.6
        (CMAKE_Fortran_COMPILER_ID STREQUAL "XL" ) )          # CMake 2.8
  if( "${CMAKE_Fortran_FLAGS}" MATCHES "-qflttrap=[a-zA-Z:]:enable" )
    set( FPE_EXIT TRUE )
  endif()

  if( NOT ("${CMAKE_Fortran_FLAGS}" MATCHES "-qfixed") )
    set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -qfixed"
         CACHE STRING "Flags for Fortran compiler." FORCE )
  endif()

else()
endif()

if( FPE_EXIT )
  message( FATAL_ERROR "Floating Point Exception (FPE) trap handlers are currently explicitly enabled in the compiler flags.  LAPACK is designed to check for and handle these cases internally and enabling these traps will likely cause LAPACK to crash.  Please re-configure with floating point exception trapping disabled." )
endif()

endmacro()
