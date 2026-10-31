# Replace FORTRAN_STRLEN definition in lapack.h if needed
if(NOT FORTRAN_STRLEN_TYPE)
  set(FORTRAN_STRLEN_TYPE "size_t")
endif()
configure_file(lapack.h.in ${LAPACK_BINARY_DIR}/include/lapack.h)
install(
  FILES ${LAPACK_BINARY_DIR}/include/lapack.h
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  COMPONENT Development
)

configure_file(lapacke_config.h ${LAPACK_BINARY_DIR}/include/lapacke_config.h COPYONLY)
install(FILES lapacke_config.h DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} COMPONENT Development)

if(LAPACKE)
  list(APPEND LAPACKE_INCLUDE lapacke.h lapacke_64.h lapacke_utils.h)
  foreach(FILE IN LISTS LAPACKE_INCLUDE)
    configure_file(${FILE} ${LAPACK_BINARY_DIR}/include/${FILE} COPYONLY)
  endforeach()
  install(
    FILES ${LAPACKE_INCLUDE}
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    COMPONENT Development
  )
endif()

# Create a header file for LAPACK routine mangling (lapacke_mangling.h)
if(CMAKE_Fortran_COMPILER_LOADED)
  include(FortranCInterface)
  ## Ensure that the fortran compiler and c compiler specified are compatible
  FortranCInterface_VERIFY()
  FortranCInterface_HEADER(${LAPACK_BINARY_DIR}/include/lapacke_mangling.h
    MACRO_NAMESPACE "LAPACK_"
    SYMBOL_NAMESPACE "LAPACK_")
endif()
if(NOT FortranCInterface_GLOBAL_FOUND)
  message(WARNING "Reverting to pre-defined include/lapacke_mangling.h")
  configure_file(lapacke_mangling_with_flags.h.in
                 ${LAPACK_BINARY_DIR}/include/lapacke_mangling.h)
endif()

install(
  FILES ${LAPACK_BINARY_DIR}/include/lapacke_mangling.h
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  COMPONENT Development
)
