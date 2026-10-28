# Replace FORTRAN_STRLEN definition in cblas_f77.h if needed
if(NOT FORTRAN_STRLEN_TYPE)
  set(FORTRAN_STRLEN_TYPE "size_t")
endif()
configure_file(cblas_f77.h.in ${LAPACK_BINARY_DIR}/include/cblas_f77.h)
install(
  FILES ${LAPACK_BINARY_DIR}/include/cblas_f77.h
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  COMPONENT Development)

if(CBLAS)
  set(CBLAS_INCLUDE cblas.h cblas_64.h)
  configure_file(cblas.h ${LAPACK_BINARY_DIR}/include/cblas.h COPYONLY)
  configure_file(cblas_64.h ${LAPACK_BINARY_DIR}/include/cblas_64.h COPYONLY)
  install(FILES ${CBLAS_INCLUDE}
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    COMPONENT Development)
endif()

# Create a header file for CBLAS routine mangling (cblas_mangling.h)
if(CMAKE_Fortran_COMPILER_LOADED)
  include(FortranCInterface)
  ## Ensure that the fortran compiler and c compiler specified are compatible
  FortranCInterface_VERIFY()
  FortranCInterface_HEADER(${LAPACK_BINARY_DIR}/include/cblas_mangling.h
                           MACRO_NAMESPACE "F77_"
                           SYMBOL_NAMESPACE "F77_")

endif()
if(NOT FortranCInterface_GLOBAL_FOUND)
  message(WARNING "Reverting to pre-defined include/cblas_mangling.h")
    configure_file(
      cblas_mangling_with_flags.h.in
      ${LAPACK_BINARY_DIR}/include/cblas_mangling.h)
endif()

install(FILES ${LAPACK_BINARY_DIR}/include/cblas_mangling.h
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  COMPONENT Development)
