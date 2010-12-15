# - Find BLAS library
# This module finds an installed fortran library that implements the BLAS
# linear-algebra interface (see http://www.netlib.org/blas/).
#
# Note: This is loosely based off of but a complete re-write of cmake's
# included FindBLAS.cmake module
#
# This module sets the following variables:
#  BLAS_FOUND        - set to true if an appropriate BLAS library is found.
#  BLAS_LINKER_FLAGS - list of required linker flags for a given BLAS library
#  BLAS_LIBRARIES    - list of libraries (using full path name) to link against
#  BLAS_STATIC       - if set then static libraries will be searched for.
#  BLAS_VENDORS      - A list of specific vendors to check for.  if not set,
#                      checks all known vendors.
#  BLAS_VENDORS_FOUND - A list of located BLAS vendors
#  BLAS_${VENDOR}_LIB_DIR - An additional library dir to search for:
#                           Ex: BLAS_AMD_LIB_DIR=/opt/acml4.4.0/gfortran64/lib
#                               BLAS_INTEL_LIB_DIR=/opt/intel/mkl/lib/intel64
##########
#
# Valid values for the BLAS_VENDOR setting are:
#   AMD      - Single threaded version of the AMD Core Math Library
#   AMD_MP   - Multithreaded version of the AMD Core Math Library using OpenMP
#            See http://developer.amd.com/cpu/Libraries/acml
#   APPLE    - Apple's Accelerate library
#            See http://developer.apple.com/performance/accelerateframework.html
#   ATLAS    - Automatically Tuned Linear Algebra Software
#            See http://math-atlas.sourceforge.net/
#   GOTO     - Goto BLAS v2
#            See http://www.tacc.utexas.edu/tacc-projects/gotoblas2
#   HP       - HP's Math Library: VECLIB
#   HP_INT64 - HP's Math Library: VECLIB8 (64 bit integers)
#            See http://www.hp.com/go/mlib
#   IBM      - IBM's Engineering and Scientific Subroutine Library
#            See http://www-03.ibm.com/systems/software/essl/
#   INTEL32       - Intel Math Kernel Library x86
#   INTEL64       - Intel Math Kernel Library x86_64
#   INTEL64_INT64 - Intel Math Kernel Library x86_64 (64 bit integers)
#               See http://software.intel.com/en-us/intel-mkl
#   NETLIB   - Reference BLAS implementation
#            See http://www.netlib.org/blas
#   ORACLE   - Oracle Performance Library (formerly Sun Performance Library)
#            See http://www.oracle.com/technetwork/server-storage/solarisstudio
#   SGI      - SGI's Scientific Computing Software Library
#            See http://www.sgi.com/products/software/irix/scsl.html
#   SUN      - Sun Performance Library (now Oracle)
#            See http://www.oracle.com/technetwork/server-storage/solarisstudio
#
#   OTHER    - Any other unsupported BLAS library
#            In order to specify other BLAS libraries, set the following:
#            BLAS_VENDORS         = OTHER
#            BLAS_OTHER_LIB_NAMES = A list of libraries to link to
#            BLAS_OTHER_LIB_DIR   = Search path for the libraries
#            Ex:
#            set(BLAS_VENDORS OTHER)
#            set(BLAS_OTHER_LIB_NAMES "my_blas;my_blas_support")
#            set(BLAS_OTHER_LIB_DIR /home/chuck/lib)
#            find_package(BLAS)
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

  if( BLAS_${__BLAS_VENDOR}_LIBRARIES )

    # Check the library as Fortran
    set( BLAS_${__BLAS_VENDOR}_LINKER_FLAGS "${__BLAS_FLAGS}" )
    if( _CHECK_FORTRAN )
      set( CMAKE_REQUIRED_LIBRARIES ${BLAS_${__BLAS_VENDOR}_LIBRARIES} )
      set( CMAKE_REQUIRED_FLAGS "${BLAS_${__BLAS_VENDOR}_LINKER_FLAGS}" )
      CHECK_FORTRAN_FUNCTION_EXISTS( "dgemm" BLAS_${__BLAS_VENDOR}_DGEMM )
      unset( CMAKE_REQUIRED_LIBRARIES )
      unset( CMAKE_REQUIRED_FLAGS )

    # Check the library as C
    #else()
    #  message( STATUS "Checking ${__BLAS_LIBNAME} for dgemm_" )
    #  get_filename_component(
    #    __BLAS_${VENDOR}_PATH
    #    ${BLAS_${__BLAS_VENDOR}_LIBRARY}
    #    PATH
    #  )
    #  CHECK_LIBRARY_EXISTS( 
    #    ${__BLAS_LIBNAME} 
    #    "dgemm_"
    #    "${__BLAS_${VENDOR}_PATH}"
    #    BLAS_${_BLAS_VENDOR}_DGEMM
    #  )
    endif()
    if( BLAS_${__BLAS_VENDOR}_DGEMM )
      set( BLAS_${__BLAS_VENDOR}_FOUND TRUE )
    endif()
  endif()
endmacro()

# Loop through the BLAS vendors looking for specific libraries
if( NOT BLAS_VENDORS )
  set( BLAS_VENDORS AMD AMD_MP ACCELERATE ATLAS GOTO VECLIB VECLIB8 IBM INTEL32 INTEL64 INTEL64_INT64 PRE-INSTALLED SGI SUN)
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
  elseif( _BLAS_VENDOR STREQUAL "AMD" )
    message( STATUS "FindBLAS: Searching for AMD ACML" )
    _BLAS_LOCATE_AND_TEST( ${_BLAS_VENDOR} "acml" "" )
  
  # Multithreaded threaded ACML
  elseif( _BLAS_VENDOR STREQUAL "AMD_MP" )
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
  elseif( _BLAS_VENDOR STREQUAL "IBM" )
    message( STATUS "FindBLAS: Searching for IBM ESSL" )
    _BLAS_LOCATE_AND_TEST( ${_BLAS_VENDOR} "essl" "" )
  
  # IBM ESSL
  elseif( _BLAS_VENDOR STREQUAL "IBM_INT64" )
    message( STATUS "FindBLAS: Searching for IBM ESSL int64" )
    _BLAS_LOCATE_AND_TEST( ${_BLAS_VENDOR} "essl8" "" )
  
  # Intel MKL
  elseif( _BLAS_VENDOR STREQUAL "INTEL32" )
    message( STATUS "FindBLAS: Searching for Intel MKL" )
    if( CMAKE_Fortran_COMPILER_ID STREQUAL "Intel" )
      _BLAS_LOCATE_AND_TEST( 
        ${_BLAS_VENDOR} "iomp5;mkl_core;mkl_intel_thread;mkl_intel_ia32" "" 
      )
    elseif( CMAKE_Fortran_COMPILER_ID STREQUAL "PGI" )
      _BLAS_LOCATE_AND_TEST( 
        ${_BLAS_VENDOR} "pgmp;mkl_core;mkl_pgi_thread;mkl_intel_ia32" "" 
      )
    else() #if( CMAKE_Fortran_COMPILER_ID STREQUAL "GNU" )
      _BLAS_LOCATE_AND_TEST( 
        ${_BLAS_VENDOR} "gomp;mkl_core;mkl_gnu_thread;mkl_intel_ia32" "" 
      )
    endif()
  
  # Intel MKL
  elseif( _BLAS_VENDOR STREQUAL "INTEL64" )
    message( STATUS "FindBLAS: Searching for Intel MKL" )
    if( CMAKE_Fortran_COMPILER_ID STREQUAL "Intel" )
      _BLAS_LOCATE_AND_TEST( 
        ${_BLAS_VENDOR} "iomp5;mkl_core;mkl_intel_thread;mkl_intel_lp64" "" 
      )
    elseif( CMAKE_Fortran_COMPILER_ID STREQUAL "PGI" )
      _BLAS_LOCATE_AND_TEST( 
        ${_BLAS_VENDOR} "pgmp;mkl_core;mkl_pgi_thread;mkl_intel_lp64" "" 
      )
    else() #if( CMAKE_Fortran_COMPILER_ID STREQUAL "GNU" )
      _BLAS_LOCATE_AND_TEST( 
        ${_BLAS_VENDOR} "gomp;mkl_core;mkl_gnu_thread;mkl_intel_lp64" "" 
      )
    endif()
  
  # Intel MKL
  elseif( _BLAS_VENDOR STREQUAL "INTEL64_INT64" )
    message( STATUS "FindBLAS: Searching for Intel MKL int64" )
    if( CMAKE_Fortran_COMPILER_ID STREQUAL "Intel" )
      _BLAS_LOCATE_AND_TEST( 
        ${_BLAS_VENDOR} "iomp5;mkl_core;mkl_intel_thread;mkl_intel_ilp64" "" 
      )
    elseif( CMAKE_Fortran_COMPILER_ID STREQUAL "PGI" )
      _BLAS_LOCATE_AND_TEST( 
        ${_BLAS_VENDOR} "pgmp;mkl_core;mkl_pgi_thread;mkl_intel_ilp64" "" 
      )
    else() #if( CMAKE_Fortran_COMPILER_ID STREQUAL "GNU" )
      _BLAS_LOCATE_AND_TEST( 
        ${_BLAS_VENDOR} "gomp;mkl_core;mkl_gnu_thread;mkl_intel_ilp64" "" 
      )
    endif()
  
  # Pre-installed BLAS
  elseif( _BLAS_VENDOR STREQUAL "PRE-INSTALLED" )
    message( STATUS "FindBLAS: Searching for pre-installed BLAS" )
    _BLAS_LOCATE_AND_TEST( ${_BLAS_VENDOR} "blas" "" )
  
  # SGI
  elseif( _BLAS_VENDOR STREQUAL "SGI" )
    message( STATUS "FindBLAS: Searching for SGI SCCL" )
    _BLAS_LOCATE_AND_TEST( ${_BLAS_VENDOR} "sccl" "" )
  
  # Sun / Oracle PerfLib
  elseif( (_BLAS_VENDOR STREQUAL "SUN") OR (_BLAS_VENDOR MATCHES "ORACLE") )
    message( STATUS "FindBLAS: Searching for Sun PerfLib" )
    _BLAS_LOCATE_AND_TEST( ${_BLAS_VENDOR} "perflib" "" )
  
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
message( STATUS "FindBLAS: BLAS Vendor selected - ${BLAS_VENDOR_FOUND}" )
set( BLAS_LIBRARIES ${BLAS_${BLAS_VENDOR_FOUND}_LIBRARIES} CACHE PATH "")
set( BLAS_LINKER_FLAGS ${BLAS_${BLAS_VENDOR_FOUND}_LINKER_FLAGS} CACHE PATH "" )
set( BLAS_FOUND TRUE CACHE OPTION "")

