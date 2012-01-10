message(STATUS "LAPACKE enable")
enable_language(C)

include(FortranMangling)
FORTRAN_MANGLING(CDEFS)
set(CDEFS ${CDEFS} CACHE STRING "Fortran Mangling" FORCE)
MESSAGE(STATUS "=========")

# --------------------------------------------------
# Compiler Flags
ADD_DEFINITIONS( "-D${CDEFS}")

if (WIN32 AND NOT UNIX)
 ADD_DEFINITIONS(-DHAVE_LAPACK_CONFIG_H -DLAPACK_COMPLEX_STRUCTURE)
 MESSAGE (STATUS "Windows BUILD")
endif (WIN32 AND NOT UNIX)

get_directory_property( DirDefs COMPILE_DEFINITIONS )
message( "COMPILE_DEFINITIONS = ${DirDefs}" )

include_directories( include )
add_subdirectory(include)
add_subdirectory(src)
add_subdirectory(utils)

macro(append_subdir_files variable dirname)
get_directory_property(holder DIRECTORY ${dirname} DEFINITION ${variable})
foreach(depfile ${holder})
  list(APPEND ${variable} "${dirname}/${depfile}")
endforeach()
endmacro()

append_subdir_files(LAPACKE_INCLUDE "include")
append_subdir_files(SRC_OBJ "src")
append_subdir_files(SRCX_OBJ "src")
append_subdir_files(MATGEN_OBJ "src")
append_subdir_files(UTILS_OBJ "utils")

if (USE_XBLAS)
   add_library(lapacke ${SRC_OBJ} ${MATGEN_OBJ} ${SRCX_OBJ} ${UTILS_OBJ})
   target_link_libraries(lapacke tmglib lapack ${BLAS_LIBRARIES}  ${XBLAS_LIBRARY})
else (USE_XBLAS)
   add_library(lapacke ${SRC_OBJ} ${MATGEN_OBJ} ${UTILS_OBJ})
   target_link_libraries(lapacke tmglib lapack ${BLAS_LIBRARIES} )
endif(USE_XBLAS)

lapack_install_library(lapacke)
INSTALL( FILES ${LAPACKE_INCLUDE} DESTINATION include )

if(BUILD_TESTING)
   add_subdirectory(example)
endif(BUILD_TESTING)


configure_file(${CMAKE_CURRENT_SOURCE_DIR}/lapacke.pc.in ${CMAKE_CURRENT_BINARY_DIR}/lapacke.pc)
 install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/lapacke.pc
  DESTINATION ${PKG_CONFIG_DIR}
   )
