message(STATUS "CBLAS enable")
enable_language(C)

set(LAPACK_INSTALL_EXPORT_NAME ${CBLASLIB}-targets)

# Create a header file cblas.h for the routines called in my C programs
include(FortranCInterface)
## Ensure that the fortran compiler and c compiler specified are compatible
FortranCInterface_VERIFY()
FortranCInterface_HEADER(${LAPACK_BINARY_DIR}/include/cblas_mangling.h
                         MACRO_NAMESPACE "F77_"
                         SYMBOL_NAMESPACE "F77_")
if(NOT FortranCInterface_GLOBAL_FOUND OR NOT FortranCInterface_MODULE_FOUND)
  message(WARNING "Reverting to pre-defined include/lapacke_mangling.h")
    configure_file(include/lapacke_mangling_with_flags.h.in
                  ${LAPACK_BINARY_DIR}/include/lapacke_mangling.h)
    configure_file(include/cblas_mangling_with_flags.h.in
                 ${LAPACK_BINARY_DIR}/include/cblas_mangling.h)
endif()

include(CheckCSourceCompiles)
check_c_source_compiles("void __attribute__((weak)) main() {};"
                        HAS_ATTRIBUTE_WEAK_SUPPORT)

include_directories(include ${LAPACK_BINARY_DIR}/include)
add_subdirectory(include)
add_subdirectory(src)

macro(append_subdir_files variable dirname)
get_directory_property(holder DIRECTORY ${dirname} DEFINITION ${variable})
foreach(depfile ${holder})
  list(APPEND ${variable} "${dirname}/${depfile}")
endforeach()
endmacro()

append_subdir_files(CBLAS_INCLUDE "include")
install(FILES ${CBLAS_INCLUDE} ${LAPACK_BINARY_DIR}/include/cblas_mangling.h
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  COMPONENT Development
  )

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
