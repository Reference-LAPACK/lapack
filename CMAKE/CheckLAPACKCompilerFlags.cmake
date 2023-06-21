# This module checks against various known compilers and their respective
# flags to determine any specific flags needing to be set.
#
#  1.  If FPE traps are enabled either abort or disable them
#  2.  Specify fixed form if needed
#  3.  Ensure that Release builds use O2 instead of O3
#
#=============================================================================
# Author: Chuck Atkins
# Copyright 2011
#=============================================================================

macro( CheckLAPACKCompilerFlags )

set( FPE_EXIT FALSE )

# FORTRAN ILP default
if ( FORTRAN_ILP )
    if( CMAKE_Fortran_COMPILER_ID STREQUAL "Intel" )
        if ( WIN32 )
            set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} /integer-size:64")
        else ()
            set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -integer-size 64")
        endif()
    elseif( (CMAKE_Fortran_COMPILER_ID STREQUAL "VisualAge" ) OR  # CMake 2.6
            (CMAKE_Fortran_COMPILER_ID STREQUAL "XL" ) )          # CMake 2.8
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -qintsize=8")
    elseif( CMAKE_Fortran_COMPILER_ID STREQUAL "NAG" )
        if ( WIN32 )
            set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} /i8")
        else ()
            set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -i8")
        endif()
    else()
        set(CPE_ENV $ENV{PE_ENV})
        if(CPE_ENV STREQUAL "CRAY")
          set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -sinteger64")
        elseif(CPE_ENV STREQUAL "NVIDIA")
          set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -i8")
        else()
          set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fdefault-integer-8")  
        endif()        
    endif()
endif()

# GNU Fortran
if( CMAKE_Fortran_COMPILER_ID STREQUAL "GNU" )
  if( "${CMAKE_Fortran_FLAGS}" MATCHES "-ffpe-trap=[izoupd]")
    set( FPE_EXIT TRUE )
  endif()

# Intel Fortran
elseif( CMAKE_Fortran_COMPILER_ID STREQUAL "Intel" )
  if( "${CMAKE_Fortran_FLAGS}" MATCHES "[-/]fpe(-all=|)0" )
    set( FPE_EXIT TRUE )
  endif()

# SunPro F95
elseif( CMAKE_Fortran_COMPILER_ID STREQUAL "SunPro" )
  if( ("${CMAKE_Fortran_FLAGS}" MATCHES "-ftrap=") AND
      NOT ("${CMAKE_Fortran_FLAGS}" MATCHES "-ftrap=(%|)none") )
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

# HP Fortran
elseif( CMAKE_Fortran_COMPILER_ID STREQUAL "HP" )
  if( "${CMAKE_Fortran_FLAGS}" MATCHES "\\+fp_exception" )
    set( FPE_EXIT TRUE )
  endif()

  if( NOT ("${CMAKE_Fortran_FLAGS}" MATCHES "\\+fltconst_strict") )
    message( STATUS "Enabling strict float conversion with +fltconst_strict" )
    set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} +fltconst_strict"
         CACHE STRING "Flags for Fortran compiler." FORCE )
  endif()

  # Most versions of cmake don't have good default options for the HP compiler
  set( CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -g"
       CACHE STRING "Flags used by the compiler during debug builds" FORCE )
  set( CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_MINSIZEREL} +Osize"
       CACHE STRING "Flags used by the compiler during release minsize builds" FORCE )
  set( CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_RELEASE} +O2"
       CACHE STRING "Flags used by the compiler during release builds" FORCE )
  set( CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_RELWITHDEBINFO} +O2 -g"
       CACHE STRING "Flags used by the compiler during release with debug info builds" FORCE )

# NAG Fortran
elseif( CMAKE_Fortran_COMPILER_ID STREQUAL "NAG" )
  if( "${CMAKE_Fortran_FLAGS}" MATCHES "[-/]ieee=(stop|nonstd)" )
    set( FPE_EXIT TRUE )
  endif()

  if( NOT ("${CMAKE_Fortran_FLAGS}" MATCHES "[-/]ieee=full") )
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ieee=full")
  endif()

  if( NOT ("${CMAKE_Fortran_FLAGS}" MATCHES "[-/]dcfuns") )
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -dcfuns")
  endif()

  if( NOT ("${CMAKE_Fortran_FLAGS}" MATCHES "[-/]thread_safe") )
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -thread_safe")
  endif()

  # Disable warnings
  if( NOT ("${CMAKE_Fortran_FLAGS}" MATCHES "[-/]w=obs") )
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -w=obs")
  endif()

  if( NOT ("${CMAKE_Fortran_FLAGS}" MATCHES "[-/]w=x77") )
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -w=x77")
  endif()

  if( NOT ("${CMAKE_Fortran_FLAGS}" MATCHES "[-/]w=ques") )
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -w=ques")
  endif()

  if( NOT ("${CMAKE_Fortran_FLAGS}" MATCHES "[-/]w=unused") )
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -w=unused")
  endif()

  # Suppress compiler banner and summary
  check_fortran_compiler_flag("-quiet" _quiet)
  if( _quiet AND NOT ("${CMAKE_Fortran_FLAGS}" MATCHES "[-/]quiet") )
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -quiet")
  endif()

else()
endif()

if( "${CMAKE_Fortran_FLAGS_RELEASE}" MATCHES "O[3-9]" )
  message( STATUS "Reducing RELEASE optimization level to O2" )
  string( REGEX REPLACE "O[3-9]" "O2" CMAKE_Fortran_FLAGS_RELEASE
          "${CMAKE_Fortran_FLAGS_RELEASE}" )
  set( CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}"
       CACHE STRING "Flags used by the compiler during release builds" FORCE )
endif()


if( FPE_EXIT )
  message( FATAL_ERROR "Floating Point Exception (FPE) trap handlers are currently explicitly enabled in the compiler flags.  LAPACK is designed to check for and handle these cases internally and enabling these traps will likely cause LAPACK to crash.  Please re-configure with floating point exception trapping disabled." )
endif()

endmacro()
