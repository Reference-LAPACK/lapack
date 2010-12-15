cmake_minimum_required(VERSION 2.6)
project(LAPACK Fortran)

# Add the CMake directory for custon CMake modules
set(CMAKE_MODULE_PATH "${LAPACK_SOURCE_DIR}/CMAKE" ${CMAKE_MODULE_PATH})

if (UNIX)
   if ( "${CMAKE_Fortran_COMPILER}" MATCHES "ifort" )
  set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fltconsistency -fp_port" )
   endif ()
endif ()

macro(lapack_install_library lib)
  install(TARGETS ${lib} EXPORT lapack-targets
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
    RUNTIME DESTINATION bin
    )
endmacro()

# --------------------------------------------------
# Testing

enable_testing()
include(CTest)
enable_testing()
# --------------------------------------------------

set(CMAKE_MODULE_PATH ${LAPACK_SOURCE_DIR} ${CMAKE_MODULE_PATH})

# --------------------------------------------------
# Check second function

include(CheckTimeFunction)
set(TIME_FUNC NONE ${TIME_FUNC})
CHECK_TIME_FUNCTION(NONE TIME_FUNC)
CHECK_TIME_FUNCTION(INT_CPU_TIME TIME_FUNC)
CHECK_TIME_FUNCTION(EXT_ETIME TIME_FUNC)
CHECK_TIME_FUNCTION(EXT_ETIME_ TIME_FUNC)
CHECK_TIME_FUNCTION(INT_ETIME TIME_FUNC)
message(STATUS "--> Will use second_${TIME_FUNC}.f and dsecnd_${TIME_FUNC}.f as timing function.")

set(SECOND_SRC  ${LAPACK_SOURCE_DIR}/INSTALL/second_${TIME_FUNC}.f)
set(DSECOND_SRC  ${LAPACK_SOURCE_DIR}/INSTALL/dsecnd_${TIME_FUNC}.f)

# --------------------------------------------------
# Subdirectories that need to be processed

option(LOOK_FOR_OPTIMZED_BLAS "Whether or not to use an optimized BLAS library or the included netlib BLAS" OFF)
if( LOOK_FOR_OPTIMZED_BLAS )
  find_package( BLAS )
endif(LOOK_FOR_OPTIMZED_BLAS)
if(NOT BLAS_FOUND)
  add_subdirectory(BLAS)
  set( BLAS_LIBRARIES blas )
else()
  set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${BLAS_LINKER_FLAGS}" 
    CACHE STRING "Flags for Fortran compiler" FORCE
  )
endif( NOT BLAS_FOUND )


add_subdirectory(SRC)
add_subdirectory(TESTING)

# --------------------------------------------------
# CPACK Packaging 

SET(CPACK_PACKAGE_NAME "LAPACK")
SET(CPACK_PACKAGE_VENDOR "University of Tennessee, Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd")
SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "LAPACK- Linear Algebra Package")
set(LAPACK_VERSION 3.2.2)
set(CPACK_PACKAGE_VERSION_MAJOR 3)
set(CPACK_PACKAGE_VERSION_MINOR 2)
set(CPACK_PACKAGE_VERSION_PATCH 2)
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE")
SET(CPACK_PACKAGE_INSTALL_DIRECTORY "LAPACK")
IF(WIN32 AND NOT UNIX)
  # There is a bug in NSI that does not handle full unix paths properly. Make
  # sure there is at least one set of four (4) backlasshes.
  SET(CPACK_NSIS_HELP_LINK "http:\\\\\\\\http://icl.cs.utk.edu/lapack-forum")
  SET(CPACK_NSIS_URL_INFO_ABOUT "http:\\\\\\\\www.netlib.org/lapack")
  SET(CPACK_NSIS_CONTACT "lapack@eecs.utk.edu")
  SET(CPACK_NSIS_MODIFY_PATH ON)
  SET(CPACK_NSIS_DISPLAY_NAME "LAPACK-${LAPACK_VERSION}")
  set(CPACK_PACKAGE_RELOCATABLE "true")
ELSE(WIN32 AND NOT UNIX)
  SET(CPACK_GENERATOR "TGZ")
  SET(CPACK_SOURCE_GENERATOR TGZ)
  SET(CPACK_SOURCE_PACKAGE_FILE_NAME "lapack-${LAPACK_VERSION}" )
  SET(CPACK_SOURCE_IGNORE_FILES ~$ .svn ${CPACK_SOURCE_IGNORE_FILES} )
ENDIF(WIN32 AND NOT UNIX)
INCLUDE(CPack)


# --------------------------------------------------
# By default static library
OPTION(BUILD_SHARED_LIBS "Build shared libraries" OFF )
#OPTION(BUILD_SHARED_LIBS "Build shared libraries" ON )

if( NOT BLAS_FOUND )
  export(TARGETS blas lapack FILE lapack-targets.cmake)
else( NOT BLAS_FOUND )
  export(TARGETS lapack FILE lapack-targets.cmake)
endif( NOT BLAS_FOUND )
configure_file(${LAPACK_SOURCE_DIR}/lapack-config-version.cmake.in
  ${LAPACK_BINARY_DIR}/lapack-config-version.cmake @ONLY)
configure_file(${LAPACK_SOURCE_DIR}/lapack-config.cmake.in
  ${LAPACK_BINARY_DIR}/lapack-config.cmake @ONLY)

configure_file(${LAPACK_SOURCE_DIR}/lapack-config-install.cmake.in
  ${LAPACK_BINARY_DIR}/CMakeFiles/lapack-config.cmake @ONLY)
install(FILES
  ${LAPACK_BINARY_DIR}/CMakeFiles/lapack-config.cmake
  ${LAPACK_BINARY_DIR}/lapack-config-version.cmake
  DESTINATION lib/cmake/lapack-${LAPACK_VERSION}
  )
install(EXPORT lapack-targets
  DESTINATION lib/cmake/lapack-${LAPACK_VERSION})
