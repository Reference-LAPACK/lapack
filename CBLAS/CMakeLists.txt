message(STATUS "CBLAS enabled")
enable_language(C)

include(CheckLanguage)
check_language(Fortran)
if(CMAKE_Fortran_COMPILER)
  enable_language(Fortran)

  # Check for any necessary platform specific compiler flags
  include(CheckLAPACKCompilerFlags)
  CheckLAPACKCompilerFlags()
endif()

set(LAPACK_INSTALL_EXPORT_NAME ${CBLASLIB}-targets)

include(CheckCSourceCompiles)
check_c_source_compiles("int __attribute__((weak)) main() {};"
                        HAS_ATTRIBUTE_WEAK_SUPPORT)

include_directories(include ${LAPACK_BINARY_DIR}/include)
add_subdirectory(src)

macro(append_subdir_files variable dirname)
get_directory_property(holder DIRECTORY ${dirname} DEFINITION ${variable})
foreach(depfile ${holder})
  list(APPEND ${variable} "${dirname}/${depfile}")
endforeach()
endmacro()

if(BUILD_TESTING)
  add_subdirectory(testing)
  add_subdirectory(examples)
endif()

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cblas.pc.in
  ${CMAKE_CURRENT_BINARY_DIR}/${CBLASLIB}.pc @ONLY)
install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/${CBLASLIB}.pc
  DESTINATION ${PKG_CONFIG_DIR}
  COMPONENT Development
  )

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cmake/cblas-config-version.cmake.in
  ${LAPACK_BINARY_DIR}/${CBLASLIB}-config-version.cmake @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cmake/cblas-config-build.cmake.in
  ${LAPACK_BINARY_DIR}/${CBLASLIB}-config.cmake @ONLY)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cmake/cblas-config-install.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/${CBLASLIB}-config.cmake @ONLY)
install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/${CBLASLIB}-config.cmake
  ${LAPACK_BINARY_DIR}/${CBLASLIB}-config-version.cmake
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${CBLASLIB}-${LAPACK_VERSION}
  COMPONENT Development
  )

install(EXPORT ${CBLASLIB}-targets
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${CBLASLIB}-${LAPACK_VERSION}
  COMPONENT Development
  )
