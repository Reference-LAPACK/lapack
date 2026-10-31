cmake_minimum_required(VERSION 3.26...4.0)
project(LAPACK VERSION 3.12.1 LANGUAGES C)

list(PREPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/CMAKE")
include(GNUInstallDirs)
include(CMakeDependentOption)
include(CMakePackageConfigHelpers)

# Build options. Keep the established cache names together here.
string(TOLOWER "${CMAKE_BUILD_TYPE}" build_type)
set(testing_default OFF)
if(build_type STREQUAL "coverage")
  set(testing_default ON)
endif()
option(BUILD_SHARED_LIBS "Build shared libraries" OFF)
option(BUILD_TESTING "Build tests" ${testing_default})

option(BUILD_SINGLE "Build single precision real" ON)
option(BUILD_DOUBLE "Build double precision real" ON)
option(BUILD_COMPLEX "Build single precision complex" ON)
option(BUILD_COMPLEX16 "Build double precision complex" ON)
option(BUILD_DEPRECATED "Build deprecated routines" OFF)

option(USE_OPTIMIZED_BLAS "Use an external BLAS implementation" OFF)
option(USE_OPTIMIZED_LAPACK "Use an external LAPACK implementation" OFF)
option(USE_XBLAS "Build extra precision routines using XBLAS" OFF)
# External providers do not normally implement Netlib's suffixed ILP64 API.
set(lapack_external_blas OFF)
set(lapack_external_lapack OFF)
if(USE_OPTIMIZED_BLAS OR BLAS_LIBRARIES)
  set(lapack_external_blas ON)
endif()
if(USE_OPTIMIZED_LAPACK OR LAPACK_LIBRARIES)
  set(lapack_external_lapack ON)
endif()
set(extended_default ON)
if(lapack_external_blas OR lapack_external_lapack)
  set(extended_default OFF)
endif()
option(BUILD_DEFAULT_API "Build the unsuffixed API" ON)
option(BUILD_INDEX64 "Use 64-bit integers in the unsuffixed API" OFF)
option(BUILD_INDEX64_EXT_API "Build the suffixed _64 API" ${extended_default})

option(CBLAS "Build the C BLAS interface" OFF)
option(LAPACKE "Build the C LAPACK interface" OFF)
option(LAPACKE_WITH_TMG "Include matrix generators in LAPACKE" OFF)
cmake_dependent_option(LAPACKE_BUILD_SINGLE "Build LAPACKE single precision real" ON "BUILD_SINGLE" OFF)
cmake_dependent_option(LAPACKE_BUILD_DOUBLE "Build LAPACKE double precision real" ON "BUILD_DOUBLE" OFF)
cmake_dependent_option(LAPACKE_BUILD_COMPLEX "Build LAPACKE single precision complex" ON "BUILD_COMPLEX" OFF)
cmake_dependent_option(LAPACKE_BUILD_COMPLEX16 "Build LAPACKE double precision complex" ON "BUILD_COMPLEX16" OFF)

option(TEST_FORTRAN_COMPILER "Test Fortran complex arithmetic" OFF)
option(USE_FLAT_NAMESPACE "Use a flat namespace on macOS" OFF)
option(LAPACK_TESTING_USE_PYTHON "Use Python for the numerical test summary" ON)
option(LAPACK_MEMORY_CHECK "Configure the tests for valgrind" OFF)
set(LAPACK_TESTING_JUNIT_JOB_ID "" CACHE STRING "Job identifier in JUnit reports")
mark_as_advanced(LAPACK_TESTING_JUNIT_JOB_ID)
option(BUILD_HTML_DOCUMENTATION "Build the HTML API documentation" OFF)
option(BUILD_MAN_DOCUMENTATION "Build the manual pages" OFF)
option(BLAS++ "Build the optional BLAS++ package" OFF)
option(LAPACK++ "Build the optional LAPACK++ package" OFF)
set(LIBRARY_PREFIX "" CACHE STRING "Prefix added to library base names")
set(LAPACK_BINARY_PATH_SUFFIX "" CACHE STRING "Subdirectory for installed libraries and binaries")

if(NOT (BUILD_SINGLE OR BUILD_DOUBLE OR BUILD_COMPLEX OR BUILD_COMPLEX16))
  message(FATAL_ERROR "Enable at least one of BUILD_SINGLE, BUILD_DOUBLE, BUILD_COMPLEX or BUILD_COMPLEX16.")
endif()
if(NOT BUILD_DEFAULT_API AND NOT BUILD_INDEX64_EXT_API)
  message(FATAL_ERROR "Enable BUILD_DEFAULT_API or BUILD_INDEX64_EXT_API.")
endif()
if(LAPACKE_WITH_TMG)
  set(LAPACKE ON)
endif()
# Recompute the platform default unless the caller explicitly supplied it.
if(NOT DEFINED LAPACK_SKIP_ERROR_EXIT_TESTS)
  set(LAPACK_SKIP_ERROR_EXIT_TESTS OFF)
  if(BUILD_SHARED_LIBS AND (WIN32 OR (APPLE AND NOT USE_FLAT_NAMESPACE)))
    set(LAPACK_SKIP_ERROR_EXIT_TESTS ON)
    message(STATUS "Disabling error-exit tests: this shared-library platform cannot override XERBLA")
  endif()
endif()

# Build layout and integer ABI.
if(CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_CURRENT_BINARY_DIR)
  message(FATAL_ERROR "Use an out-of-source build: cmake -S . -B build")
endif()
if(CMAKE_INSTALL_PREFIX STREQUAL CMAKE_CURRENT_BINARY_DIR)
  message(FATAL_ERROR "The installation prefix must differ from the build directory.")
endif()
if(PROJECT_IS_TOP_LEVEL AND NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build configuration" FORCE)
endif()
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${LAPACK_BINARY_DIR}/bin")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${LAPACK_BINARY_DIR}/lib")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${LAPACK_BINARY_DIR}/lib")
set(lapack_integer_size 4)
set(lapack_library_suffix "")
if(BUILD_INDEX64)
  set(lapack_integer_size 8)
  set(lapack_library_suffix 64)
endif()
set(BLASLIB blas${lapack_library_suffix})
set(CBLASLIB cblas${lapack_library_suffix})
set(LAPACKLIB lapack${lapack_library_suffix})
set(LAPACKELIB lapacke${lapack_library_suffix})
set(TMGLIB tmglib${lapack_library_suffix})

# Only LAPACKE with both external providers can use a C-only toolchain.
if(NOT lapack_external_blas OR NOT lapack_external_lapack OR
    BUILD_TESTING OR LAPACKE_WITH_TMG OR
    TEST_FORTRAN_COMPILER)
  enable_language(Fortran)
else()
  include(CheckLanguage)
  check_language(Fortran)
  if(CMAKE_Fortran_COMPILER)
    enable_language(Fortran)
  endif()
endif()

# Private requirements shared by libraries, object libraries and test drivers.
add_library(lapack-build-options INTERFACE)
include(CheckLAPACKCompilerFlags)
lapack_compiler_flags(lapack-build-options)
target_include_directories(lapack-build-options INTERFACE "${LAPACK_BINARY_DIR}/include")
if(BUILD_INDEX64)
  target_compile_definitions(lapack-build-options INTERFACE
    "$<$<COMPILE_LANGUAGE:C>:WeirdNEC;LAPACK_ILP64;HAVE_LAPACK_CONFIG_H>")
endif()
if(APPLE AND USE_FLAT_NAMESPACE)
  target_link_options(lapack-build-options INTERFACE "LINKER:-flat_namespace")
endif()
if(build_type STREQUAL "coverage")
  find_package(codecov REQUIRED)
endif()
function(lapack_configure_target target)
  target_link_libraries(${target} PRIVATE $<BUILD_LOCAL_INTERFACE:lapack-build-options>)
  set_target_properties(${target} PROPERTIES POSITION_INDEPENDENT_CODE ON
    C_STANDARD 99 C_STANDARD_REQUIRED ON C_EXTENSIONS OFF)
  if(build_type STREQUAL "coverage")
    add_coverage(${target})
    if(CMAKE_C_COMPILER_ID STREQUAL "GNU" OR CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
      target_link_libraries(${target} PRIVATE gcov)
    endif()
  endif()
endfunction()
include(ExtendedAPIHelpers)

# External providers. Raw library lists are confined to this compatibility boundary.
set(lapack_blas_target BLAS::${BLASLIB})
set(lapack_lapack_target LAPACK::${LAPACKLIB})
set(lapack_blas_package ${BLASLIB})
set(lapack_lapack_package ${LAPACKLIB})
set(probe_precision z)
foreach(spec "COMPLEX;c" "SINGLE;s" "DOUBLE;d")
  list(GET spec 0 precision)
  list(GET spec 1 letter)
  if(BUILD_${precision})
    set(probe_precision ${letter})
  endif()
endforeach()
foreach(provider BLAS LAPACK)
  string(TOLOWER "${provider}" name)
  if(NOT lapack_external_${name})
    continue()
  endif()
  if(NOT TARGET ${provider}::${provider})
    if(${provider}_LIBRARIES)
      add_library(${provider}::${provider} INTERFACE IMPORTED)
      set_target_properties(${provider}::${provider} PROPERTIES
        INTERFACE_LINK_LIBRARIES "${${provider}_LIBRARIES}")
      separate_arguments(link_flags NATIVE_COMMAND "${${provider}_LINKER_FLAGS}")
      set_property(TARGET ${provider}::${provider} PROPERTY INTERFACE_LINK_OPTIONS "${link_flags}")
    else()
      block(PROPAGATE ${provider}_LIBRARIES ${provider}_LINKER_FLAGS)
        set(BLA_SIZEOF_INTEGER ${lapack_integer_size})
        find_package(${provider} REQUIRED MODULE)
      endblock()
    endif()
  endif()
  set(symbol ${probe_precision}gemm)
  if(provider STREQUAL "LAPACK")
    set(symbol ${probe_precision}geqrt)
  endif()
  set(symbols)
  if(BUILD_DEFAULT_API)
    list(APPEND symbols ${symbol})
  endif()
  if(BUILD_INDEX64_EXT_API)
    list(APPEND symbols ${symbol}_64)
  endif()
  foreach(symbol IN LISTS symbols)
    if(CMAKE_Fortran_COMPILER_LOADED)
      set(probe_file check.f90)
      set(probe "program check\nexternal ${symbol}\ncall ${symbol}()\nend program")
    else()
      set(probe_file check.c)
      set(probe "extern void ${symbol}_(void);\nint main(void) { ${symbol}_(); return 0; }")
    endif()
    try_compile(works SOURCE_FROM_CONTENT ${probe_file} "${probe}"
      LINK_LIBRARIES ${provider}::${provider} NO_CACHE)
    if(NOT works)
      message(FATAL_ERROR "The requested ${provider} cannot link ${symbol}; check the provider, integer ABI and BUILD_INDEX64_EXT_API.")
    endif()
  endforeach()
  set(lapack_${name}_target ${provider}::${provider})
  set(lapack_${name}_package ${provider})
endforeach()
if(USE_XBLAS)
  find_library(XBLAS_LIBRARY NAMES xblas REQUIRED)
  add_library(LAPACK::XBLAS UNKNOWN IMPORTED)
  set_target_properties(LAPACK::XBLAS PROPERTIES IMPORTED_LOCATION "${XBLAS_LIBRARY}")
endif()

# Canonical targets. Each source directory contributes to these libraries.
# IntelLLVM otherwise expands a large LAPACKE response file past link.exe's limit.
if(WIN32 AND BUILD_SHARED_LIBS AND CMAKE_C_COMPILER_ID STREQUAL "IntelLLVM"
    AND CMAKE_GENERATOR MATCHES "^Ninja")
  set(CMAKE_C_RESPONSE_FILE_LINK_FLAG "/Qoption,link,@")
  set(CMAKE_C_USE_RESPONSE_FILE_FOR_LIBRARIES OFF)
endif()
set(lapack_libraries)
if(NOT lapack_external_blas)
  add_library(${BLASLIB})
  add_library(BLAS::${BLASLIB} ALIAS ${BLASLIB})
  if(NOT TARGET BLAS::BLAS)
    add_library(BLAS::BLAS ALIAS ${BLASLIB})
  endif()
  list(APPEND lapack_libraries ${BLASLIB})
endif()
if(NOT lapack_external_lapack)
  add_library(${LAPACKLIB})
  add_library(LAPACK::${LAPACKLIB} ALIAS ${LAPACKLIB})
  if(NOT TARGET LAPACK::LAPACK)
    add_library(LAPACK::LAPACK ALIAS ${LAPACKLIB})
  endif()
  target_link_libraries(${LAPACKLIB} PRIVATE ${lapack_blas_target})
  if(USE_XBLAS)
    target_link_libraries(${LAPACKLIB} PRIVATE LAPACK::XBLAS)
  endif()
  list(APPEND lapack_libraries ${LAPACKLIB})
endif()
if(CBLAS)
  add_library(${CBLASLIB})
  add_library(CBLAS::${CBLASLIB} ALIAS ${CBLASLIB})
  target_link_libraries(${CBLASLIB} PRIVATE ${lapack_blas_target})
  list(APPEND lapack_libraries ${CBLASLIB})
endif()
if(BUILD_TESTING OR LAPACKE_WITH_TMG)
  add_library(${TMGLIB})
  add_library(LAPACK::${TMGLIB} ALIAS ${TMGLIB})
  target_link_libraries(${TMGLIB} PRIVATE ${lapack_lapack_target} ${lapack_blas_target})
  list(APPEND lapack_libraries ${TMGLIB})
endif()
if(LAPACKE)
  add_library(${LAPACKELIB})
  add_library(LAPACKE::${LAPACKELIB} ALIAS ${LAPACKELIB})
  target_link_libraries(${LAPACKELIB} PRIVATE ${lapack_lapack_target})
  if(LAPACKE_WITH_TMG)
    target_link_libraries(${LAPACKELIB} PRIVATE LAPACK::${TMGLIB})
  endif()
  list(APPEND lapack_libraries ${LAPACKELIB})
endif()
foreach(target IN LISTS lapack_libraries)
  set_target_properties(${target} PROPERTIES VERSION ${PROJECT_VERSION}
    SOVERSION ${PROJECT_VERSION_MAJOR} OUTPUT_NAME "${LIBRARY_PREFIX}${target}"
    WINDOWS_EXPORT_ALL_SYMBOLS ON)
  target_include_directories(${target} PUBLIC
    $<BUILD_INTERFACE:${LAPACK_BINARY_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
  if(BUILD_INDEX64)
    target_compile_definitions(${target} INTERFACE
      "$<$<COMPILE_LANGUAGE:C,CXX>:WeirdNEC;LAPACK_ILP64;HAVE_LAPACK_CONFIG_H>")
  endif()
  if(NOT BUILD_SHARED_LIBS)
    target_link_libraries(${target} PRIVATE ${lapack_fortran_runtime})
  endif()
  if(TEST_FORTRAN_COMPILER)
    add_dependencies(${target} run_test_zcomplexabs run_test_zcomplexdiv
      run_test_zcomplexmult run_test_zminMax)
  endif()
endforeach()

# Sources and tests.
if(BUILD_TESTING)
  enable_testing()
  if(LAPACK_MEMORY_CHECK)
    find_program(MEMORYCHECK_COMMAND valgrind REQUIRED)
    set(MEMORYCHECK_COMMAND_OPTIONS "--leak-check=full --show-leak-kinds=definite,indirect --track-origins=yes --trace-children=yes")
    set(LAPACK_TESTING_USE_PYTHON OFF)
    include(CTest) # Generates the configuration needed by ctest -T memcheck.
  endif()
  include(LAPACKTestHelpers)
  if(LAPACK_TESTING_USE_PYTHON)
    find_package(Python3 COMPONENTS Interpreter)
  endif()
  configure_file(CTestCustom.cmake.in CTestCustom.cmake @ONLY)
endif()
add_subdirectory(CBLAS/include)
add_subdirectory(LAPACKE/include)
if(NOT lapack_external_blas)
  add_subdirectory(BLAS)
endif()
if(NOT lapack_external_lapack)
  include(CheckTimeFunction)
  lapack_check_time_function(TIME_FUNC)
  set(SECOND_SRC "${LAPACK_SOURCE_DIR}/INSTALL/second_${TIME_FUNC}.f")
  set(DSECOND_SRC "${LAPACK_SOURCE_DIR}/INSTALL/dsecnd_${TIME_FUNC}.f")
  add_subdirectory(SRC)
endif()
if(TARGET ${TMGLIB})
  add_subdirectory(TESTING/MATGEN)
endif()
if(CBLAS)
  add_subdirectory(CBLAS)
endif()
if(BUILD_TESTING OR TEST_FORTRAN_COMPILER)
  add_subdirectory(TESTING)
endif()
if(LAPACKE)
  add_subdirectory(LAPACKE)
endif()
foreach(target IN LISTS lapack_libraries)
  lapack_configure_target(${target})
endforeach()
if(BUILD_TESTING AND Python3_EXECUTABLE AND LAPACK_TESTING_USE_PYTHON)
  message(STATUS "Enabling LAPACK test summary (see TESTING/testing_results.txt and lapack_testing_junit.xml)")

  # Every test registered anywhere below the source root, so that adding a
  # test directory never needs a matching edit here. Collected before the
  # summary itself is added, so that it does not depend on itself.
  function(_collect_tests dir out)
    get_property(_tests DIRECTORY "${dir}" PROPERTY TESTS)
    get_property(_subdirs DIRECTORY "${dir}" PROPERTY SUBDIRECTORIES)
    foreach(_subdir IN LISTS _subdirs)
      _collect_tests("${_subdir}" _subdir_tests)
      list(APPEND _tests ${_subdir_tests})
    endforeach()
    set(${out} "${_tests}" PARENT_SCOPE)
  endfunction()

  _collect_tests("${LAPACK_SOURCE_DIR}" _all_tests)

  # Without a job identifier the report keeps the plain suite and class names.
  set(_summary_job_option "")
  if(LAPACK_TESTING_JUNIT_JOB_ID)
    set(_summary_job_option "--junit-job" "${LAPACK_TESTING_JUNIT_JOB_ID}")
  endif()

  file(COPY ${LAPACK_SOURCE_DIR}/lapack_testing.py DESTINATION ${LAPACK_BINARY_DIR})
  # The test drivers report a numerical failure in their output and still exit
  # with a status of 0, so this is the test that fails the run when one of them
  # found something: --fail-on-error turns any failure or error they reported
  # into a nonzero exit status of the summary.
  add_test(
    NAME LAPACK_Test_Summary
    WORKING_DIRECTORY ${LAPACK_BINARY_DIR}
    COMMAND ${Python3_EXECUTABLE} "lapack_testing.py"
      "--fail-on-error" "--fail-if-empty" "--fail-on-empty-output"
      "--fail-on-unrecognized" "--merge-apis"
      "--junit-xml" "lapack_testing_junit.xml" ${_summary_job_option}
    )

  if(_all_tests)
    set_property(TEST LAPACK_Test_Summary APPEND PROPERTY DEPENDS ${_all_tests})
  endif()
endif()

# Each package keeps its own config and export; both trees use the same template.
# The only nonstandard binary layout is the existing optional path suffix.
string(REGEX REPLACE "^/+" "" binary_suffix "${LAPACK_BINARY_PATH_SUFFIX}")
function(lapack_install_library target namespace config_template)
  set(package_dir "${CMAKE_INSTALL_LIBDIR}/cmake/${target}")
  if(binary_suffix)
    string(APPEND CMAKE_INSTALL_LIBDIR "/${binary_suffix}")
    string(APPEND CMAKE_INSTALL_BINDIR "/${binary_suffix}")
  endif()
  install(TARGETS ${target} EXPORT ${target}-targets
    ARCHIVE COMPONENT Development LIBRARY COMPONENT RuntimeLibraries
    RUNTIME COMPONENT RuntimeLibraries)
  configure_package_config_file(${config_template} "${target}-config.cmake"
    INSTALL_DESTINATION "${package_dir}")
  write_basic_package_version_file("${target}-config-version.cmake"
    VERSION ${LAPACK_VERSION} COMPATIBILITY SameMajorVersion)
  export(EXPORT ${target}-targets NAMESPACE ${namespace}:: FILE "${target}-targets.cmake")
  install(EXPORT ${target}-targets NAMESPACE ${namespace}::
    DESTINATION "${package_dir}" COMPONENT Development)
  install(FILES "${LAPACK_BINARY_DIR}/${target}-config.cmake"
    "${LAPACK_BINARY_DIR}/${target}-config-version.cmake"
    DESTINATION "${package_dir}" COMPONENT Development)
endfunction()
if(NOT lapack_external_blas)
  lapack_install_library(${BLASLIB} BLAS CMAKE/blas-config.cmake.in)
endif()
if(NOT lapack_external_lapack)
  lapack_install_library(${LAPACKLIB} LAPACK CMAKE/lapack-config.cmake.in)
endif()
if(CBLAS)
  lapack_install_library(${CBLASLIB} CBLAS CBLAS/cmake/cblas-config.cmake.in)
endif()
if(TARGET ${TMGLIB})
  lapack_install_library(${TMGLIB} LAPACK CMAKE/tmglib-config.cmake.in)
endif()
if(LAPACKE)
  lapack_install_library(${LAPACKELIB} LAPACKE LAPACKE/cmake/lapacke-config.cmake.in)
endif()

# pkg-config metadata, independent of the CMake target exports above.
set(${LAPACKLIB}_dependencies ${lapack_blas_package})
set(${CBLASLIB}_dependencies ${lapack_blas_package})
set(${LAPACKELIB}_dependencies ${lapack_lapack_package})
set(${TMGLIB}_dependencies ${lapack_lapack_package} ${lapack_blas_package})
if(USE_XBLAS)
  list(APPEND ${LAPACKLIB}_dependencies XBLAS)
endif()
if(LAPACKE_WITH_TMG)
  list(APPEND ${LAPACKELIB}_dependencies ${TMGLIB})
endif()
set(install_libdir "${CMAKE_INSTALL_LIBDIR}")
if(binary_suffix)
  string(APPEND install_libdir "/${binary_suffix}")
endif()
set(pc_libdir "${install_libdir}")
set(pc_includedir "${CMAKE_INSTALL_INCLUDEDIR}")
foreach(dir pc_libdir pc_includedir)
  if(NOT IS_ABSOLUTE "${${dir}}")
    set(${dir} "\${prefix}/${${dir}}")
  endif()
endforeach()
file(RELATIVE_PATH pc_prefix "${CMAKE_INSTALL_FULL_LIBDIR}/pkgconfig" "${CMAKE_INSTALL_PREFIX}")
set(pc_definitions "")
if(BUILD_INDEX64)
  set(pc_definitions " -DLAPACK_ILP64 -DWeirdNEC -DHAVE_LAPACK_CONFIG_H")
endif()
foreach(target IN LISTS lapack_libraries)
  # pkg-config consumes the raw lists provided by FindBLAS/FindLAPACK; CMake
  # targets above never depend on this conversion.
  set(pc_requires "")
  set(pc_private "")
  set(pc_libraries)
  foreach(dependency IN LISTS ${target}_dependencies)
    if(dependency MATCHES "^(BLAS|LAPACK)$" AND ${dependency}_LIBRARIES)
      list(APPEND pc_libraries ${${dependency}_LIBRARIES})
      string(APPEND pc_private " ${${dependency}_LINKER_FLAGS}")
    elseif(dependency STREQUAL "XBLAS")
      list(APPEND pc_libraries "${XBLAS_LIBRARY}")
    else()
      string(TOLOWER "${dependency}" pc_dependency)
      string(APPEND pc_requires " ${pc_dependency}")
    endif()
  endforeach()
  list(APPEND pc_libraries ${lapack_fortran_runtime})
  foreach(library IN LISTS pc_libraries)
    if(IS_ABSOLUTE "${library}")
      string(APPEND pc_private " \"${library}\"")
    elseif(library MATCHES "^-")
      string(APPEND pc_private " ${library}")
    else()
      string(APPEND pc_private " -l${library}")
    endif()
  endforeach()
  configure_file(CMAKE/lapack.pc.in "${target}.pc" @ONLY)
  install(FILES "${LAPACK_BINARY_DIR}/${target}.pc"
    DESTINATION "${CMAKE_INSTALL_LIBDIR}/pkgconfig" COMPONENT Development)
endforeach()

if(BUILD_HTML_DOCUMENTATION OR BUILD_MAN_DOCUMENTATION)
  add_subdirectory(DOCS)
endif()
# Keep the existing C++ integration until its dependency handling is migrated.
block()
  if(NOT lapack_external_blas)
    set(BLAS_LIBRARIES ${BLASLIB})
  endif()
  if(NOT lapack_external_lapack)
    set(LAPACK_LIBRARIES ${LAPACKLIB})
  endif()
function(_display_cpp_implementation_msg name)
  string(TOLOWER ${name} name_lc)
  message(STATUS "${name}++ enable")
  message(STATUS "----------------")
  message(STATUS "Thank you for your interest in ${name}++, a newly developed C++ API for ${name} library")
  message(STATUS "The objective of ${name}++ is to provide a convenient, performance oriented API for development in the C++ language, that, for the most part, preserves established conventions, while, at the same time, takes advantages of modern C++ features, such as: namespaces, templates, exceptions, etc.")
  message(STATUS "For support ${name}++ related question, please email: slate-user@icl.utk.edu")
  message(STATUS "----------------")
endfunction()
if (BLAS++)
  _display_cpp_implementation_msg("BLAS")
  include(ExternalProject)
  ExternalProject_Add(blaspp
	  URL https://bitbucket.org/icl/blaspp/downloads/blaspp-2020.10.02.tar.gz
	  CONFIGURE_COMMAND ${CMAKE_COMMAND} -E env LIBRARY_PATH=$ENV{LIBRARY_PATH}:${CMAKE_BINARY_DIR}/lib LD_LIBRARY_PATH=$ENV{LD_LIBRARY_PATH}:${PROJECT_BINARY_DIR}/lib ${CMAKE_COMMAND} -DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR} -DCMAKE_INSTALL_LIBDIR=lib -DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS} ${PROJECT_BINARY_DIR}/blaspp-prefix/src/blaspp
	  BUILD_COMMAND ${CMAKE_COMMAND} -E env LIBRARY_PATH=$ENV{LIBRARY_PATH}:${PROJECT_BINARY_DIR}/lib LIB_SUFFIX="" ${CMAKE_COMMAND} --build .
	  INSTALL_COMMAND ${CMAKE_COMMAND} -E env PREFIX=${PROJECT_BINARY_DIR} LIB_SUFFIX="" ${CMAKE_COMMAND} --install .
	  )
  ExternalProject_Add_StepDependencies(blaspp build ${BLAS_LIBRARIES})
endif()
if (LAPACK++)
	message (STATUS "linking lapack++ against ${LAPACK_LIBRARIES}")
  _display_cpp_implementation_msg("LAPACK")
  include(ExternalProject)
  if (BUILD_SHARED_LIBS)
  ExternalProject_Add(lapackpp
	  URL https://bitbucket.org/icl/lapackpp/downloads/lapackpp-2020.10.02.tar.gz
	  CONFIGURE_COMMAND ${CMAKE_COMMAND} -E env LIBRARY_PATH=$ENV{LIBRARY_PATH}:${CMAKE_BINARY_DIR}/lib LD_LIBRARY_PATH=$ENV{LD_LIBRARY_PATH}:${PROJECT_BINARY_DIR}/lib ${CMAKE_COMMAND} -DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR} -DCMAKE_INSTALL_LIBDIR=lib -DLAPACK_LIBRARIES=${LAPACK_LIBRARIES} -DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS}  ${PROJECT_BINARY_DIR}/lapackpp-prefix/src/lapackpp
	  BUILD_COMMAND ${CMAKE_COMMAND} -E env LIBRARY_PATH=$ENV{LIBRARY_PATH}:${PROJECT_BINARY_DIR}/lib LIB_SUFFIX="" ${CMAKE_COMMAND} --build .
	  INSTALL_COMMAND ${CMAKE_COMMAND} -E env PREFIX=${PROJECT_BINARY_DIR} LIB_SUFFIX="" ${CMAKE_COMMAND} --install .
	  )
  else ()
# FIXME this does not really work as the libraries list gets converted to a semicolon-separated list somewhere in the lapack++ build files
  ExternalProject_Add(lapackpp
	  URL https://bitbucket.org/icl/lapackpp/downloads/lapackpp-2020.10.02.tar.gz
	  CONFIGURE_COMMAND env LIBRARY_PATH=$ENV{LIBRARY_PATH}:${CMAKE_BINARY_DIR}/lib LD_LIBRARY_PATH=$ENV{LD_LIBRARY_PATH}:${PROJECT_BINARY_DIR}/lib ${CMAKE_COMMAND} -DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR} -DCMAKE_INSTALL_LIBDIR=lib -DLAPACK_LIBRARIES="${PROJECT_BINARY_DIR}/lib/liblapack.a -lgfortran"  -DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS}  ${PROJECT_BINARY_DIR}/lapackpp-prefix/src/lapackpp
	  BUILD_COMMAND env LIBRARY_PATH=$ENV{LIBRARY_PATH}:${PROJECT_BINARY_DIR}/lib LIB_SUFFIX="" ${CMAKE_COMMAND} --build .
	  INSTALL_COMMAND ${CMAKE_COMMAND} -E env PREFIX=${PROJECT_BINARY_DIR} LIB_SUFFIX="" ${CMAKE_COMMAND} --install .
  )
  endif()
  ExternalProject_Add_StepDependencies(lapackpp build blaspp ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
endif()

# --------------------------------------------------
# CPACK Packaging
if (LAPACK++)
  install(
  DIRECTORY "${LAPACK_BINARY_DIR}/lib/"
  DESTINATION "${CMAKE_INSTALL_LIBDIR}${LAPACK_BINARY_PATH_SUFFIX}"
  FILES_MATCHING REGEX "liblapackpp.(a|so)$"
  )
  install(
    DIRECTORY "${PROJECT_BINARY_DIR}/lapackpp-prefix/src/lapackpp/include/"
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
    FILES_MATCHING REGEX "\\.(h|hh)$"
  )
  write_basic_package_version_file(
    "lapackppConfigVersion.cmake"
    VERSION 2020.10.02
    COMPATIBILITY AnyNewerVersion
  )
  install(
    FILES "${CMAKE_CURRENT_BINARY_DIR}/lib/lapackpp/lapackppConfig.cmake"
          "${CMAKE_CURRENT_BINARY_DIR}/lib/lapackpp/lapackppConfigVersion.cmake"
    DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/"
  )

endif()
if (BLAS++)
  write_basic_package_version_file(
    "blasppConfigVersion.cmake"
    VERSION 2020.10.02
    COMPATIBILITY AnyNewerVersion
  )
  install(
    FILES "${CMAKE_CURRENT_BINARY_DIR}/lib/blaspp/blasppConfig.cmake"
          "${CMAKE_CURRENT_BINARY_DIR}/lib/blaspp/blasppConfigVersion.cmake"
    DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/"
  )
  install(
  DIRECTORY "${LAPACK_BINARY_DIR}/lib/"
  DESTINATION "${CMAKE_INSTALL_LIBDIR}${LAPACK_BINARY_PATH_SUFFIX}"
  FILES_MATCHING REGEX "libblaspp.(a|so)$"
  )
  install(
    DIRECTORY "${PROJECT_BINARY_DIR}/blaspp-prefix/src/blaspp/include/"
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
    FILES_MATCHING REGEX "\\.(h|hh)$"
  )
endif()
endblock()

# Standard CPack packaging, matching the ScaLAPACK build.
set(CPACK_PACKAGE_NAME "LAPACK")
set(CPACK_PACKAGE_VENDOR "University of Tennessee, University of California Berkeley, University of Colorado Denver and NAG Ltd")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "LAPACK - Linear Algebra Package")
set(CPACK_PACKAGE_VERSION ${LAPACK_VERSION})
set(CPACK_RESOURCE_FILE_LICENSE "${LAPACK_SOURCE_DIR}/LICENSE")
set(CPACK_PACKAGE_INSTALL_DIRECTORY "LAPACK")
set(CPACK_SOURCE_PACKAGE_FILE_NAME "lapack-${LAPACK_VERSION}")
if(NOT WIN32)
  set(CPACK_GENERATOR TGZ)
  set(CPACK_SOURCE_GENERATOR TGZ)
endif()
include(CPack)
