# Replace FORTRAN_STRLEN definition in cblas_f77.h if needed
file(READ cblas_f77.h cblas_f77_h_content)
if(FORTRAN_STRLEN_TYPE)
  string(REPLACE
    "#define FORTRAN_STRLEN size_t" "#define FORTRAN_STRLEN ${FORTRAN_STRLEN_TYPE}"
    cblas_f77_h_content "${cblas_f77_h_content}")
endif()
file(WRITE ${LAPACK_BINARY_DIR}/include/cblas_f77.h "${cblas_f77_h_content}")
install(
  FILES ${LAPACK_BINARY_DIR}/include/cblas_f77.h
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  COMPONENT Development)

if(CBLAS)
  set(CBLAS_INCLUDE cblas.h cblas_64.h)
  file(COPY ${CBLAS_INCLUDE} DESTINATION ${LAPACK_BINARY_DIR}/include)
  install(FILES ${CBLAS_INCLUDE}
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    COMPONENT Development)

  file(COPY cblas_test.h DESTINATION ${LAPACK_BINARY_DIR}/test_include)
endif()

# Create a header file for CBLAS routine mangling (cblas_mangling.h)
include(CheckLanguage)
check_language(Fortran)
check_language(C)
if(CMAKE_Fortran_COMPILER AND CMAKE_C_COMPILER)
  enable_language(Fortran)
  enable_language(C)
  include(FortranCInterface)
  ## Ensure that the fortran compiler and c compiler specified are compatible
  FortranCInterface_VERIFY()
  FortranCInterface_HEADER(${LAPACK_BINARY_DIR}/include/cblas_mangling.h
                           MACRO_NAMESPACE "F77_"
                           SYMBOL_NAMESPACE "F77_")

  # Check for any necessary platform specific compiler flags
  include(CheckLAPACKCompilerFlags)
  CheckLAPACKCompilerFlags()
endif()
if(NOT FortranCInterface_GLOBAL_FOUND OR NOT FortranCInterface_MODULE_FOUND)
  message(WARNING "Reverting to pre-defined include/cblas_mangling.h")
    configure_file(
      include/cblas_mangling_with_flags.h.in
      ${LAPACK_BINARY_DIR}/include/cblas_mangling.h)
endif()

install(FILES ${LAPACK_BINARY_DIR}/include/cblas_mangling.h
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  COMPONENT Development)
