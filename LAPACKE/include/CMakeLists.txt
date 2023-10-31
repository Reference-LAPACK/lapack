set(LAPACKE_INCLUDE lapack.h)
IF(LAPACKE)
    list(APPEND LAPACKE_INCLUDE lapacke.h lapacke_config.h lapacke_utils.h)
endif()

file(COPY ${LAPACKE_INCLUDE} DESTINATION ${LAPACK_BINARY_DIR}/include)

# Create a header file lapacke_mangling.h for the routines called in my C programs
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