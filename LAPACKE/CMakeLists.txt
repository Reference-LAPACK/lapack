message(STATUS "LAPACKE enabled")
enable_language(C)

set(LAPACK_INSTALL_EXPORT_NAME ${LAPACKELIB}-targets)

include_directories(${LAPACK_BINARY_DIR}/include include)
add_subdirectory(src)
add_subdirectory(utils)

option(LAPACKE_BUILD_SINGLE "Build LAPACKE single precision real" ON)
option(LAPACKE_BUILD_DOUBLE "Build LAPACKE double precision real" ON)
option(LAPACKE_BUILD_COMPLEX "Build LAPACKE single precision complex" ON)
option(LAPACKE_BUILD_COMPLEX16 "Build LAPACKE double precision complex" ON)

foreach(precision SINGLE DOUBLE COMPLEX COMPLEX16)
  if(LAPACKE_BUILD_${precision} AND NOT BUILD_${precision})
    message(
      SEND_ERROR
      "LAPACKE_BUILD_${precision} requires BUILD_${precision} to be ON"
    )
  endif()
endforeach()

macro(append_subdir_files variable dirname)
  get_directory_property(holder DIRECTORY ${dirname} DEFINITION ${variable})
  foreach(depfile ${holder})
    list(APPEND ${variable} "${dirname}/${depfile}")
  endforeach()
endmacro()

message(STATUS "Build LAPACKE single precision real: ${LAPACKE_BUILD_SINGLE}")
message(STATUS "Build LAPACKE double precision real: ${LAPACKE_BUILD_DOUBLE}")
message(STATUS "Build LAPACKE single precision complex: ${LAPACKE_BUILD_COMPLEX}")
message(STATUS "Build LAPACKE double precision complex: ${LAPACKE_BUILD_COMPLEX16}")

append_subdir_files(COMMON_SOURCES "src")
append_subdir_files(SOURCES "src")
if (LAPACKE_BUILD_SINGLE)
  append_subdir_files(SOURCES_SINGLE "src")
  list(APPEND SOURCES ${SOURCES_SINGLE})
endif()
if (LAPACKE_BUILD_DOUBLE)
  append_subdir_files(SOURCES_DOUBLE "src")
  list(APPEND SOURCES ${SOURCES_DOUBLE})
endif()
if (LAPACKE_BUILD_COMPLEX)
  append_subdir_files(SOURCES_COMPLEX "src")
  list(APPEND SOURCES ${SOURCES_COMPLEX})
endif()
if (LAPACKE_BUILD_COMPLEX16)
  append_subdir_files(SOURCES_COMPLEX16 "src")
  list(APPEND SOURCES ${SOURCES_COMPLEX16})
endif()
append_subdir_files(DEPRECATED "src")
append_subdir_files(EXTENDED "src")
append_subdir_files(MATGEN "src")
append_subdir_files(UTILS "utils")
append_subdir_files(COMMON_UTILS "utils")

if(BUILD_DEPRECATED)
  list(APPEND SOURCES ${DEPRECATED})
endif()

set(XSOURCES)
if(USE_XBLAS)
  list(APPEND XSOURCES ${EXTENDED})
endif()
if(LAPACKE_WITH_TMG)
  list(APPEND SOURCES ${MATGEN})
endif()
list(APPEND SOURCES ${UTILS})

if(BUILD_DEFAULT_API)
  add_library(${LAPACKELIB}_obj OBJECT
    ${COMMON_SOURCES} ${COMMON_UTILS} ${SOURCES} ${XSOURCES})
else()
  add_library(${LAPACKELIB}_obj OBJECT ${COMMON_SOURCES} ${COMMON_UTILS})
endif()

if(BUILD_INDEX64_EXT_API)
  add_library(${LAPACKELIB}_64_obj OBJECT ${SOURCES})
  target_compile_definitions(${LAPACKELIB}_64_obj PRIVATE
    LAPACK_ILP64 LAPACKE_API64 WeirdNEC CBLAS_API64)
endif()

add_library(${LAPACKELIB}
  $<TARGET_OBJECTS:${LAPACKELIB}_obj>
  $<$<BOOL:${BUILD_INDEX64_EXT_API}>: $<TARGET_OBJECTS:${LAPACKELIB}_64_obj>>)

set_target_properties(
  ${LAPACKELIB} PROPERTIES
  LINKER_LANGUAGE C
  VERSION ${LAPACK_VERSION}
  SOVERSION ${LAPACK_MAJOR_VERSION}
  )
target_include_directories(${LAPACKELIB} PUBLIC
    $<BUILD_INTERFACE:${LAPACK_BINARY_DIR}/include>
  $<INSTALL_INTERFACE:include>
)

if(LAPACKE_WITH_TMG)
  target_link_libraries(${LAPACKELIB} PRIVATE ${TMGLIB})
endif()
target_link_libraries(${LAPACKELIB} PRIVATE ${LAPACK_LIBRARIES})

lapack_install_library(${LAPACKELIB})

if(BUILD_TESTING)
  add_subdirectory(example)
endif()


configure_file(${CMAKE_CURRENT_SOURCE_DIR}/lapacke.pc.in ${CMAKE_CURRENT_BINARY_DIR}/${LAPACKELIB}.pc @ONLY)
install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/${LAPACKELIB}.pc
  DESTINATION ${PKG_CONFIG_DIR}
  COMPONENT Development
  )

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cmake/lapacke-config-version.cmake.in
  ${LAPACK_BINARY_DIR}/${LAPACKELIB}-config-version.cmake @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cmake/lapacke-config-build.cmake.in
  ${LAPACK_BINARY_DIR}/${LAPACKELIB}-config.cmake @ONLY)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cmake/lapacke-config-install.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/${LAPACKELIB}-config.cmake @ONLY)
install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/${LAPACKELIB}-config.cmake
  ${LAPACK_BINARY_DIR}/${LAPACKELIB}-config-version.cmake
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${LAPACKELIB}-${LAPACK_VERSION}
  COMPONENT Development
  )

install(EXPORT ${LAPACKELIB}-targets
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${LAPACKELIB}-${LAPACK_VERSION}
  COMPONENT Development
  )
