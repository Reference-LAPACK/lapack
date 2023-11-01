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
set(FOPT_ILP64)
if( CMAKE_Fortran_COMPILER_ID MATCHES "Intel" )
    if ( WIN32 )
        set(FOPT_ILP64 /integer-size:64)
    else ()
        set(FOPT_ILP64 "-integer-size 64")
    endif()
elseif( (CMAKE_Fortran_COMPILER_ID STREQUAL "VisualAge" ) OR  # CMake 2.6
        (CMAKE_Fortran_COMPILER_ID STREQUAL "XL" ) )          # CMake 2.8
    set(FOPT_ILP64 -qintsize=8)
elseif( CMAKE_Fortran_COMPILER_ID STREQUAL "NAG" )
    if ( WIN32 )
        set(FOPT_ILP64 /i8)
    else ()
        set(FOPT_ILP64 -i8)
    endif()
elseif( CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC" )
    if ( WIN32 )
        set(FOPT_ILP64 /i8)
    else ()
        set(FOPT_ILP64 -i8)
    endif()
else()
    set(CPE_ENV $ENV{PE_ENV})
    if(CPE_ENV STREQUAL "CRAY")
      set(FOPT_ILP64 -sinteger64)
    elseif(CPE_ENV STREQUAL "NVIDIA")
      set(FOPT_ILP64 -i8)
    else()
      set(FOPT_ILP64 -fdefault-integer-8)
    endif()
endif()
if ( FORTRAN_ILP )
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${FOPT_ILP64}")
endif()

# GNU Fortran
if( CMAKE_Fortran_COMPILER_ID STREQUAL "GNU" )
  if( "${CMAKE_Fortran_FLAGS}" MATCHES "-ffpe-trap=[izoupd]")
    set( FPE_EXIT TRUE )
  endif()
  if( NOT ("${CMAKE_Fortran_FLAGS}" MATCHES "-frecursive") )
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -frecursive"
      CACHE STRING "Recursive flag must be set" FORCE)
  endif()

# Intel Fortran
elseif( CMAKE_Fortran_COMPILER_ID MATCHES "Intel" )
  if( "${CMAKE_Fortran_FLAGS}" MATCHES "[-/]fpe(-all=|)0" )
    set( FPE_EXIT TRUE )
  endif()

  if( NOT ("${CMAKE_Fortran_FLAGS}" MATCHES "-recursive") )
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -recursive"
      CACHE STRING "Recursive flag must be set" FORCE)
  endif()

  if( UNIX AND NOT ("${CMAKE_Fortran_FLAGS}" MATCHES "-fp-model[ \t]strict") )
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fp-model strict")
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

  if(UNIX)
    # Delete libmtsk in linking sequence for Sun/Oracle Fortran Compiler.
    # This library is not present in the Sun package SolarisStudio12.3-linux-x86-bin
    string(REPLACE \;mtsk\; \; CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES "${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES}")
  endif()

# IBM XL Fortran
elseif( (CMAKE_Fortran_COMPILER_ID STREQUAL "VisualAge" ) OR  # CMake 2.6
        (CMAKE_Fortran_COMPILER_ID STREQUAL "XL" ) )          # CMake 2.8
  if( "${CMAKE_Fortran_FLAGS}" MATCHES "-qflttrap=[a-zA-Z:]:enable" )
    set( FPE_EXIT TRUE )
  endif()

  if( NOT ("${CMAKE_Fortran_FLAGS}" MATCHES "-qrecur") )
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -qrecur"
      CACHE STRING "Recursive flag must be set" FORCE)
  endif()

  if( UNIX AND NOT ("${CMAKE_Fortran_FLAGS}" MATCHES "-qnosave") )
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -qnosave")
  endif()


  if( UNIX AND NOT ("${CMAKE_Fortran_FLAGS}" MATCHES "-qstrict") )
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -qstrict")
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

  if( NOT ("${CMAKE_Fortran_FLAGS}" MATCHES "-recursive") )
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -recursive"
      CACHE STRING "Recursive flag must be set" FORCE)
  endif()

  # Suppress compiler banner and summary
  check_fortran_compiler_flag("-quiet" _quiet)
  if( _quiet AND NOT ("${CMAKE_Fortran_FLAGS}" MATCHES "[-/]quiet") )
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -quiet")
  endif()

# NVIDIA HPC SDK
elseif( CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC" )
  if( ("${CMAKE_Fortran_FLAGS}" MATCHES "-Ktrap=") AND
      NOT ("${CMAKE_Fortran_FLAGS}" MATCHES "-Ktrap=none") )
    set( FPE_EXIT TRUE )
  endif()

  if( NOT ("${CMAKE_Fortran_FLAGS}" MATCHES "[-/]Kieee") )
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Kieee")
  endif()

  if( NOT ("${CMAKE_Fortran_FLAGS}" MATCHES "-Mrecursive") )
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Mrecursive"
      CACHE STRING "Recursive flag must be set" FORCE)
  endif()

# Flang Fortran
elseif( CMAKE_Fortran_COMPILER_ID STREQUAL "Flang" )
  if( NOT ("${CMAKE_Fortran_FLAGS}" MATCHES "-Mrecursive") )
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Mrecursive"
      CACHE STRING "Recursive flag must be set" FORCE)
  endif()

# Compaq Fortran
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "Compaq")
  if(WIN32)
    if(CMAKE_GENERATOR STREQUAL "NMake Makefiles")
      get_filename_component(CMAKE_Fortran_COMPILER_CMDNAM ${CMAKE_Fortran_COMPILER} NAME_WE)
      message(STATUS "Using Compaq Fortran compiler with command name ${CMAKE_Fortran_COMPILER_CMDNAM}")
      set(cmd ${CMAKE_Fortran_COMPILER_CMDNAM})
      string(TOLOWER "${cmd}" cmdlc)
      if(cmdlc STREQUAL "df")
        message(STATUS "Assume the Compaq Visual Fortran Compiler is being used")
        set(CMAKE_Fortran_USE_RESPONSE_FILE_FOR_OBJECTS 1)
        set(CMAKE_Fortran_USE_RESPONSE_FILE_FOR_INCLUDES 1)
        #This is a workaround that is needed to avoid forward-slashes in the
        #filenames listed in response files from incorrectly being interpreted as
        #introducing compiler command options
        if(${BUILD_SHARED_LIBS})
          message(FATAL_ERROR "Making of shared libraries with CVF has not been tested.")
        endif()
        set(str "NMake version 9 or later should be used. NMake version 6.0 which is\n")
        set(str "${str}   included with the CVF distribution fails to build Lapack because\n")
        set(str "${str}   the number of source files exceeds the limit for NMake v6.0\n")
        message(STATUS ${str})
        set(CMAKE_Fortran_LINK_EXECUTABLE "LINK /out:<TARGET> <LINK_FLAGS> <LINK_LIBRARIES> <OBJECTS>")
      endif()
    endif()
  endif()

else()
  message(WARNING "Fortran local arrays should be allocated on the stack."
    " Please use a compiler which guarantees that feature."
    " See https://github.com/Reference-LAPACK/lapack/pull/188 and references therein.")
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
