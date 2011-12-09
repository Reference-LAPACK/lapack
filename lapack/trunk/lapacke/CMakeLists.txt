message(STATUS "LAPACKE enable")
enable_language(C)
if (WIN32 AND NOT UNIX)
 ADD_DEFINITIONS(-DHAVE_LAPACK_CONFIG_H -DLAPACK_COMPLEX_STRUCTURE -DLAPACK_NAME_PATTERN_UC)
 MESSAGE (STATUS "Windows BUILD")
endif (WIN32 AND NOT UNIX)

include_directories( include )
#add_subdirectory(include)
add_subdirectory(src)
add_subdirectory(utils)

macro(append_subdir_files variable dirname)
get_directory_property(holder DIRECTORY ${dirname} DEFINITION ${variable})
foreach(depfile ${holder})
  list(APPEND ${variable} "${depfile}")
endforeach()
endmacro()


append_subdir_files(SRC_OBJ "src")
append_subdir_files(UTILS_OBJ "utils")

add_library(lapacke ${SRC_OBJ} ${UTILS_OBJ})
target_link_libraries(lapacke ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${XBLAS_LIBRARY})
lapack_install_library(lapacke)

if(BUILD_TESTING)
   add_subdirectory(example)
endif(BUILD_TESTING)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/lapacke.pc.in ${CMAKE_CURRENT_BINARY_DIR}/lapacke.pc)
 install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/lapacke.pc
  DESTINATION ${PKG_CONFIG_DIR}
   )
