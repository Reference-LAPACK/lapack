# Replace FORTRAN_STRLEN definition in lapack.h if needed
file(READ lapack.h lapack_h_content)
if(FORTRAN_STRLEN_TYPE)
  string(REPLACE
    "#define FORTRAN_STRLEN size_t" "#define FORTRAN_STRLEN ${FORTRAN_STRLEN_TYPE}"
    lapack_h_content "${lapack_h_content}")
endif()
file(WRITE ${LAPACK_BINARY_DIR}/include/lapack.h "${lapack_h_content}")
install(
  FILES ${LAPACK_BINARY_DIR}/include/lapack.h
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  COMPONENT Development
)

if(LAPACKE)
  list(APPEND LAPACKE_INCLUDE lapacke.h lapacke_64.h lapacke_config.h lapacke_utils.h)
  file(COPY ${LAPACKE_INCLUDE} DESTINATION ${LAPACK_BINARY_DIR}/include)
  install(
    FILES ${LAPACKE_INCLUDE}
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    COMPONENT Development
  )
endif()

# Create a header file for LAPACK routine mangling (lapacke_mangling.h)
include(CheckLanguage)
check_language(Fortran)
check_language(C)
if(CMAKE_Fortran_COMPILER AND CMAKE_C_COMPILER)
  enable_language(Fortran)
  enable_language(C)
  include(FortranCInterface)
  ## Ensure that the fortran compiler and c compiler specified are compatible
  FortranCInterface_VERIFY()
  FortranCInterface_HEADER(${LAPACK_BINARY_DIR}/include/lapacke_mangling.h
    MACRO_NAMESPACE "LAPACK_"
    SYMBOL_NAMESPACE "LAPACK_")
endif()
if(NOT FortranCInterface_GLOBAL_FOUND OR NOT FortranCInterface_MODULE_FOUND)
  message(WARNING "Reverting to pre-defined include/lapacke_mangling.h")
  configure_file(lapacke_mangling_with_flags.h.in
                 ${LAPACK_BINARY_DIR}/include/lapacke_mangling.h)
endif()

install(
  FILES ${LAPACK_BINARY_DIR}/include/lapacke_mangling.h
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  COMPONENT Development
)
