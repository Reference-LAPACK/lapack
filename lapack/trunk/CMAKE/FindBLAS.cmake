# - Find BLAS library
# This module finds an installed fortran library that implements the BLAS
# linear-algebra interface (see http://www.netlib.org/blas/).
#
# Note: This is loosely based off of but a complete re-write of cmake's
# included FindBLAS.cmake module
#
# This module sets the following variables:
#  BLAS_FOUND             - set to true if an appropriate BLAS library is found.
#  BLAS_LINKER_FLAGS      - list of required linker flags
#  BLAS_LIBRARIES         - list of libraries (with paths) to link against
#  BLAS_STATIC            - if set then static libraries will be searched for.
#  BLAS_VENDORS           - A list of specific vendors implemented BLAS to check
#                           for, if not set, checks all known vendors.
#  BLAS_VENDORS_FOUND     - A list of located BLAS vendors
#  BLAS_${VENDOR}_LIB_DIR - An additional library dir to search for:
#                           Ex: BLAS_ACML_LIB_DIR=/opt/acml4.4.0/gfortran64/lib
#                               BLAS_MKL_LP64_LIB_DIR=/opt/intel/mkl/lib/intel64
##########
#
# Valid values for the BLAS_VENDOR setting are:
#   ACML      - Single threaded version of the AMD Core Math Library
#   ACML_MP   - Multithreaded version of the AMD Core Math Library using OpenMP
#             See http://developer.amd.com/cpu/Libraries/acml
#   ACCELERATE - Apple's Accelerate library
#              See http://developer.apple.com/performance/accelerateframework
#   ATLAS      - Automatically Tuned Linear Algebra Software
#             See http://math-atlas.sourceforge.net/
#   GOTO      - Goto BLAS v2
#             See http://www.tacc.utexas.edu/tacc-projects/gotoblas2
#   GENERIC   - Search for a generic libblas
#   VECLIB    - HP's Math Library: VECLIB
#   VECLIB8   - HP's Math Library: VECLIB8 (64 bit integers)
#             See http://www.hp.com/go/mlib
#   ESSL      - IBM's Engineering and Scientific Subroutine Library
#   ESSL_6464 - IBM's Engineering and Scientific Subroutine Library (int64)
#   ESSL_SMP  - IBM's Engineering and Scientific Subroutine Library (smp)
#   ESSL_SMP_6464 - (smp + int64)
#             See http://www-03.ibm.com/systems/software/essl/
#   MKL       - Intel Math Kernel Library (dynamic interface)
#             Use MKL_THREADING_LAYER and MKL_INTERFACE_LAYER environment vars
#
#   MKL_IA32  - Intel Math Kernel Library x86
#   MKL_LP64  - Intel Math Kernel Library emt64/ia64
#   MKL_ILP64 - Intel Math Kernel Library (emt64/ia64 + int64)
#             See http://software.intel.com/en-us/intel-mkl
#   PERFLIB   - Oracle Performance Library (formerly Sun Performance Library)
#   SUNPERF   - Oracle Performance Library (formerly Sun Performance Library)
#             See http://www.oracle.com/technetwork/server-storage/solarisstudio
#   SCCL      - SGI's Scientific Computing Software Library
#             See http://www.sgi.com/products/software/irix/scsl.html
#
#   OTHER     - Any other unsupported BLAS library
#             In order to specify other BLAS libraries, set the following:
#             BLAS_VENDORS         = OTHER
#             BLAS_OTHER_LIB_NAMES = A list of libraries to link to
#             BLAS_OTHER_LIB_DIR   = Search path for the libraries
#             Ex:
#             set(BLAS_VENDORS OTHER)
#             set(BLAS_OTHER_LIB_NAMES "my_blas;my_blas_support")
#             set(BLAS_OTHER_LIB_DIR /home/chuck/lib)
#             find_package(BLAS)
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
#

# Early exit if already found
if( BLAS_FOUND )
  return()
endif()

include( CheckFortranFunctionExists )
include( CheckLibraryExists )

# Check the language being used
get_property( _LANGUAGES_ GLOBAL PROPERTY ENABLED_LANGUAGES )
if( _LANGUAGES_ MATCHES Fortran )
  set( _CHECK_FORTRAN TRUE )
elseif( (_LANGUAGES_ MATCHES C) OR (_LANGUAGES_ MATCHES CXX) )
  set( _CHECK_FORTRAN FALSE )
else()
  if(BLAS_FIND_REQUIRED)
    message(FATAL_ERROR "FindBLAS requires Fortran, C, or C++ to be enabled.")
  else(BLAS_FIND_REQUIRED)
    message(STATUS "Looking for BLAS... - NOT found (Unsupported languages)")
    return()
  endif(BLAS_FIND_REQUIRED)
endif( )


# Set the library suffix to look for
if( BLAS_STATIC )
  if( WIN32 )
    set( CMAKE_FIND_LIBRARY_SUFFIXES ".lib" )
  else()
    set( CMAKE_FIND_LIBRARY_SUFFIXES ".a" )
  endif()
else()
  if( WIN32 )
    set( CMAKE_FIND_LIBRARY_SUFFIXES ".dll;.lib" )
  elseif(APPLE)
    set( CMAKE_FIND_LIBRARY_SUFFIXES ".dylib" )
  else()
    set( CMAKE_FIND_LIBRARY_SUFFIXES ".so" )
  endif()
endif()

# Set extra library dirs
if( WIN32 )
  set( _BLAS_EXTRA_LIB_DIRS $ENV{LIB} )
elseif( APPLE )
  string( REPLACE ":" ";" _BLAS_EXTRA_LIB_DIRS "$ENV{DYLD_LIBRARY_PATH}" )
else()
  string( REPLACE ":" ";" _BLAS_EXTRA_LIB_DIRS "$ENV{LD_LIBRARY_PATH}" )
endif()

# Macro to locate a library and check for a specified symbol
macro( _BLAS_LOCATE_AND_TEST __BLAS_VENDOR __BLAS_LIBNAMES __BLAS_FLAGS )
  set( BLAS_${__BLAS_VENDOR}_LIBRARIES )
  foreach( __BLAS_LIBNAME ${__BLAS_LIBNAMES} )
    message( STATUS "FindBLAS: Searching for ${__BLAS_VENDOR} ${__BLAS_LIBNAME} - " )
    find_library( BLAS_${__BLAS_VENDOR}_${__BLAS_LIBNAME}_LIBRARY
      NAMES ${__BLAS_LIBNAME}
      PATHS ${BLAS_${__BLAS_VENDOR}_LIB_DIR} ${_BLAS_EXTRA_LIB_DIRS}
    )
    message( STATUS "FindBLAS: Searching for ${__BLAS_VENDOR} ${__BLAS_LIBNAME} - ${BLAS_${__BLAS_VENDOR}_${__BLAS_LIBNAME}_LIBRARY}" )
    if( NOT BLAS_${__BLAS_VENDOR}_${__BLAS_LIBNAME}_LIBRARY )
      unset( BLAS_${__BLAS_VENDOR}_LIBRARIES )
      break()
    endif()
    set( BLAS_${__BLAS_VENDOR}_LIBRARIES
      ${BLAS_${__BLAS_VENDOR}_LIBRARIES}
      ${BLAS_${__BLAS_VENDOR}_${__BLAS_LIBNAME}_LIBRARY}
    )
  endforeach()

  if( BLAS_${__BLAS_VENDOR}_LIBRARIES OR ("${__BLAS_LIBNAMES}" STREQUAL "") )

    # Check the library as Fortran
    set( BLAS_${__BLAS_VENDOR}_LINKER_FLAGS "${__BLAS_FLAGS}" )
    set( CMAKE_REQUIRED_LIBRARIES "${BLAS_${__BLAS_VENDOR}_LIBRARIES}")
    if( BLAS_${__BLAS_VENDOR}_LINKER_FLAGS )
      set( CMAKE_REQUIRED_LIBRARIES 
        "${BLAS_${__BLAS_VENDOR}_LINKER_FLAGS} ${CMAKE_REQUIRED_LIBRARIES}")
    endif()
    CHECK_FORTRAN_FUNCTION_EXISTS( "dgemm" BLAS_${__BLAS_VENDOR}_DGEMM )
    unset( CMAKE_REQUIRED_LIBRARIES )

    if( BLAS_${__BLAS_VENDOR}_DGEMM )
      set( BLAS_${__BLAS_VENDOR}_FOUND TRUE )
    endif()
  endif()
endmacro()

# Loop through the BLAS vendors looking for specific libraries
if( NOT BLAS_VENDORS )
  # If not specified, we will search through the list of known suppliers
  # Note that for libs that contains both 4 and 8 byte int versions, only the
  # 4 byte versions are searched for.
  set( BLAS_VENDORS ACML ACCELERATE ATLAS GOTO VECLIB ESSL MKL PERFLIB SCCL GENERIC)
endif()
set( BLAS_VENDORS_FOUND )
foreach( _BLAS_VENDOR ${BLAS_VENDORS} )
  
  # Other BLAS Library
  if( _BLAS_VENDOR STREQUAL "OTHER" )
    if( NOT BLAS_${_BLAS_VENDOR}_FLAGS )
        set( BLAS_${_BLAS_VENDOR}_FLAGS )
    endif()
    message( STATUS "FindBLAS: Searching for user specified BLAS" )
    _BLAS_LOCATE_AND_TEST( 
      ${_BLAS_VENDOR}
      "${BLAS_${_BLAS_VENDOR}_LIB_NAMES}"
      "${BLAS_${_BLAS_VENDOR}_FLAGS}"
    )
  
  # Single threaded ACML
  elseif( _BLAS_VENDOR STREQUAL "ACML" )
    message( STATUS "FindBLAS: Searching for AMD ACML" )
    _BLAS_LOCATE_AND_TEST( ${_BLAS_VENDOR} "acml" "" )
  
  # Multithreaded threaded ACML
  elseif( _BLAS_VENDOR STREQUAL "ACML_MP" )
    message( STATUS "FindBLAS: Searching for AMD ACML MP" )
    _BLAS_LOCATE_AND_TEST( ${_BLAS_VENDOR} "acml_mp" "" )
  
  # Apple Accelerate 
  elseif( _BLAS_VENDOR STREQUAL "ACCELERATE" )
    message( STATUS "FindBLAS: Searching for Apple Accelerate" )
    _BLAS_LOCATE_AND_TEST( 
      ${_BLAS_VENDOR} "Accelerate" "-framework Accelerate" 
    )
  
  # ATLAS
  elseif( _BLAS_VENDOR STREQUAL "ATLAS" )
    message( STATUS "FindBLAS: Searching for ATLAS BLAS" )
    _BLAS_LOCATE_AND_TEST( ${_BLAS_VENDOR} "f77blas;atlas" "" )
  
  # GotoBLAS2
  elseif( _BLAS_VENDOR STREQUAL "GOTO" )
    message( STATUS "FindBLAS: Searching for GotoBLAS2" )
    _BLAS_LOCATE_AND_TEST( ${_BLAS_VENDOR} "goto2" "" )
  
  # VECLIB
  elseif( _BLAS_VENDOR STREQUAL "VECLIB" )
    message( STATUS "FindBLAS: Searching for VECLIB" )
    _BLAS_LOCATE_AND_TEST( ${_BLAS_VENDOR} "veclib" "" )
  
  # VECLIB8
  elseif( _BLAS_VENDOR STREQUAL "VECLIB8" )
    message( STATUS "FindBLAS: Searching for VECLIB8" )
    _BLAS_LOCATE_AND_TEST( ${_BLAS_VENDOR} "veclib8" "" )
  
  # IBM ESSL
  elseif( _BLAS_VENDOR STREQUAL "ESSL" )
    message( STATUS "FindBLAS: Searching for IBM ESSL" )
    _BLAS_LOCATE_AND_TEST( ${_BLAS_VENDOR} "essl" "" )
  
  # IBM ESSL (SMP Version)
  elseif( _BLAS_VENDOR STREQUAL "ESSL_SMP" )
    message( STATUS "FindBLAS: Searching for IBM ESSL (SMP)" )
    _BLAS_LOCATE_AND_TEST( ${_BLAS_VENDOR} "esslsmp" "" )
  
  # IBM ESSL int64
  elseif( _BLAS_VENDOR STREQUAL "ESSL_6464" )
    message( STATUS "FindBLAS: Searching for IBM ESSL (int64)" )
    _BLAS_LOCATE_AND_TEST( ${_BLAS_VENDOR} "essl6464" "" )
  
  # IBM ESSL (SMP + int64)
  elseif( _BLAS_VENDOR STREQUAL "ESSL_SMP_6464" )
    message( STATUS "FindBLAS: Searching for IBM ESSL (SMP + int64)" )
    _BLAS_LOCATE_AND_TEST( ${_BLAS_VENDOR} "esslsmp6464" "" )
  
  # Intel MKL (dynamic)
  elseif( _BLAS_VENDOR STREQUAL "MKL" )
    message( STATUS "FindBLAS: Searching for Intel MKL (dynamic)" )
    _BLAS_LOCATE_AND_TEST( ${_BLAS_VENDOR} "mkl_rt" "" )
  
  # Intel MKL (x86)
  elseif( _BLAS_VENDOR STREQUAL "MKL_IA32" )
    message( STATUS "FindBLAS: Searching for Intel MKL (x86)" )
    if( CMAKE_Fortran_COMPILER_ID STREQUAL "Intel" )
      _BLAS_LOCATE_AND_TEST( 
        ${_BLAS_VENDOR} 
        "mkl_intel;mkl_intel_thread;mkl_core;iomp5;pthread" 
        "" 
      )
    elseif( CMAKE_Fortran_COMPILER_ID STREQUAL "PGI" )
      _BLAS_LOCATE_AND_TEST( 
        ${_BLAS_VENDOR} "mkl_intel_ia32;mkl_pgi_thread;mkl_core;pgmp" "" 
      )
    else() #if( CMAKE_Fortran_COMPILER_ID STREQUAL "GNU" )
      _BLAS_LOCATE_AND_TEST( 
        ${_BLAS_VENDOR} "mkl_intel_ia32;mkl_gnu_thread;mkl_core" "-fopenmp" 
      )
    endif()
  
  # Intel MKL (emt64 / ia64)
  elseif( _BLAS_VENDOR STREQUAL "MKL_LP64" )
    message( STATUS "FindBLAS: Searching for Intel MKL (emt64/ia64)" )
    if( CMAKE_Fortran_COMPILER_ID STREQUAL "Intel" )
      _BLAS_LOCATE_AND_TEST( 
        ${_BLAS_VENDOR} 
        "mkl_intel_lp64;mkl_intel_thread;mkl_core;iomp5;pthread" 
        "" 
      )
    elseif( CMAKE_Fortran_COMPILER_ID STREQUAL "PGI" )
      _BLAS_LOCATE_AND_TEST( 
        ${_BLAS_VENDOR} "mkl_intel_lp64;mkl_pgi_thread;mkl_core;pgmp" "" 
       )
    else() #if( CMAKE_Fortran_COMPILER_ID STREQUAL "GNU" )
      _BLAS_LOCATE_AND_TEST( 
        ${_BLAS_VENDOR} "mkl_intel_lp64;mkl_gnu_thread;mkl_core" "-fopenmp" 
       )
    endif()
  
  # Intel MKL (emt64/ia64 + int64)
  elseif( _BLAS_VENDOR STREQUAL "MKL_ILP64" )
    message( STATUS "FindBLAS: Searching for Intel MKL (emt64/ia64 + int64)" )
    if( CMAKE_Fortran_COMPILER_ID STREQUAL "Intel" )
      _BLAS_LOCATE_AND_TEST( 
        ${_BLAS_VENDOR} 
        "mkl_intel_ilp64;mkl_intel_thread;mkl_core;iomp5;pthread" 
        "" 
      )
    elseif( CMAKE_Fortran_COMPILER_ID STREQUAL "PGI" )
      _BLAS_LOCATE_AND_TEST( 
        ${_BLAS_VENDOR} "mkl_intel_ilp64;mkl_pgi_thread;mkl_core;pgmp" "" 
       )
    else() #if( CMAKE_Fortran_COMPILER_ID STREQUAL "GNU" )
      _BLAS_LOCATE_AND_TEST( 
        ${_BLAS_VENDOR} "mkl_intel_ilp64;mkl_gnu_thread;mkl_core" "-fopenmp" 
      )
    endif()
  
  # Generic BLAS
  elseif( _BLAS_VENDOR STREQUAL "GENERIC" )
    message( STATUS "FindBLAS: Searching for generic BLAS" )
    _BLAS_LOCATE_AND_TEST( ${_BLAS_VENDOR} "blas" "" )
  
  # SGI
  elseif( _BLAS_VENDOR STREQUAL "SCCL" )
    message( STATUS "FindBLAS: Searching for SGI SCCL" )
    _BLAS_LOCATE_AND_TEST( ${_BLAS_VENDOR} "sccl" "" )
  
  # Sun / Oracle PerfLib
  elseif( (_BLAS_VENDOR STREQUAL "PERFLIB") OR 
          (_BLAS_VENDOR STREQUAL "SUNPERF") )
    message( STATUS "FindBLAS: Searching for Sun PerfLib" )
    _BLAS_LOCATE_AND_TEST( ${_BLAS_VENDOR} "" "-xlic_lib=sunperf" )
  
  else()
  endif()

  if( BLAS_${_BLAS_VENDOR}_FOUND )
    set( BLAS_VENDORS_FOUND ${BLAS_VENDORS_FOUND} ${_BLAS_VENDOR} )
  endif()
endforeach()

# Parse the search results
message( STATUS "FindBLAS: BLAS vendors found: ${BLAS_VENDORS_FOUND}" )
list( LENGTH BLAS_VENDORS_FOUND _BLAS_VENDORS_FOUND_LENGTH )
if( _BLAS_VENDORS_FOUND_LENGTH EQUAL 0 )
  message( STATUS "FindBLAS: BLAS library not found" )
  return()
endif()
list( GET BLAS_VENDORS_FOUND 0 BLAS_VENDOR_FOUND )
message( STATUS "FindBLAS: BLAS vendor selected: ${BLAS_VENDOR_FOUND}" )
set( BLAS_LIBRARIES ${BLAS_${BLAS_VENDOR_FOUND}_LIBRARIES} CACHE PATH "")
set( BLAS_LINKER_FLAGS ${BLAS_${BLAS_VENDOR_FOUND}_LINKER_FLAGS} CACHE PATH "" )
set( BLAS_FOUND TRUE CACHE OPTION "")

