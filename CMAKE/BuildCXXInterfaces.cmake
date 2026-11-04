# Keep the released opt-in C++ integrations separate from the reference build.
if(NOT BUILD_DEFAULT_API)
  message(FATAL_ERROR "BLAS++ and LAPACK++ require BUILD_DEFAULT_API.")
endif()
if((NOT lapack_external_blas OR (LAPACK++ AND NOT lapack_external_lapack))
    AND NOT (BUILD_SINGLE AND BUILD_DOUBLE AND
             BUILD_COMPLEX AND BUILD_COMPLEX16))
  message(FATAL_ERROR "BLAS++ and LAPACK++ require all four provider precisions.")
endif()
include(ExternalProject)
enable_language(CXX)
set(cxx_stage "${LAPACK_BINARY_DIR}/cxx-install")
set(cxx_modules "${LAPACK_BINARY_DIR}/cxx-modules")
file(MAKE_DIRECTORY "${cxx_modules}")
foreach(package BLAS LAPACK)
  string(TOLOWER "${package}" name)
  if(lapack_external_${name})
    set(provider_code "include(\"\${CMAKE_ROOT}/Modules/Find${package}.cmake\")")
  else()
    set(library "${${package}LIB}")
    set(provider_code "find_package(${library} CONFIG REQUIRED)\nset(${package}_FOUND TRUE)\nset(${package}_LIBRARIES ${package}::${library})")
  endif()
  file(WRITE "${cxx_modules}/Find${package}.cmake" "${provider_code}\n")
endforeach()
set(cxx_integer int)
if(BUILD_INDEX64)
  set(cxx_integer int64_t)
endif()
set(cxx_dependencies)
foreach(target ${lapack_blas_target} ${lapack_lapack_target})
  get_target_property(imported ${target} IMPORTED)
  if(NOT imported)
    list(APPEND cxx_dependencies ${target})
  endif()
endforeach()
set(cxx_packages blaspp)
if(LAPACK++)
  list(APPEND cxx_packages lapackpp)
endif()
foreach(package IN LISTS cxx_packages)
  # The dependency on reference targets precedes *configuration*: upstream
  # performs link probes while configuring its own project.
  ExternalProject_Add(${package}
    URL https://bitbucket.org/icl/${package}/downloads/${package}-2020.10.02.tar.gz
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
    LIST_SEPARATOR "|"
    CMAKE_ARGS
      "-DCMAKE_INSTALL_PREFIX=${cxx_stage}"
      "-DCMAKE_INSTALL_LIBDIR=${CMAKE_INSTALL_LIBDIR}"
      "-DCMAKE_INSTALL_INCLUDEDIR=${CMAKE_INSTALL_INCLUDEDIR}"
      "-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}"
      "-DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}"
      "-DCMAKE_PREFIX_PATH=${LAPACK_BINARY_DIR}|${cxx_stage}"
      "-DCMAKE_MODULE_PATH=${cxx_modules}"
      "-DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS}"
      "-DBLA_SIZEOF_INTEGER=${lapack_integer_size}"
      "-DBLA_VENDOR=${BLA_VENDOR}"
      "-Dblas_int=${cxx_integer}"
      "-Dbuild_tests=OFF"
      "-Duse_cmake_find_blas=ON"
      "-Duse_cmake_find_lapack=ON"
    DEPENDS ${cxx_dependencies})
  set(cxx_dependency BLAS)
  if(package STREQUAL "lapackpp")
    set(cxx_dependency LAPACK)
  endif()
  string(TOLOWER "${cxx_dependency}" dependency_name)
  if(lapack_external_${dependency_name})
    set(cxx_find_dependency "find_dependency(${cxx_dependency} MODULE)")
  else()
    set(library "${${cxx_dependency}LIB}")
    set(cxx_find_dependency "find_dependency(${library} CONFIG HINTS \"\${CMAKE_CURRENT_LIST_DIR}/../cmake/${library}\")")
  endif()
  configure_file("${LAPACK_SOURCE_DIR}/CMAKE/LAPACKCXXInstall.cmake.in"
    "${LAPACK_BINARY_DIR}/${package}-install-config.cmake" @ONLY)
  ExternalProject_Add_Step(${package} dependencies
    COMMAND ${CMAKE_COMMAND} -P "${LAPACK_BINARY_DIR}/${package}-install-config.cmake"
    DEPENDEES install)
  list(APPEND cxx_dependencies ${package})
endforeach()
# Let upstream own its headers, platform-specific library names and exports.
install(DIRECTORY "${cxx_stage}/" DESTINATION . COMPONENT Development)
