# This module checks against various known compilers and their respective
# flags to determine any specific flags needing to be set.
#
#  1.  If FPE traps are enabled either abort or disable them
#  2.  Specify fixed form if needed
#  3.  Ensure that Release builds use O2 instead of O3
#
#=============================================================================
# Author: Chuck Atkins
# Copyright 2011
#=============================================================================

function(lapack_compiler_flags target)
  set(FORTRAN_STRLEN_TYPE size_t)

  # FORTRAN ILP default
  set(FOPT_ILP64)
  if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    if(WIN32)
      set(FOPT_ILP64 /integer-size:64)
    else()
      set(FOPT_ILP64 "SHELL:-integer-size 64")
    endif()
  elseif((CMAKE_Fortran_COMPILER_ID STREQUAL "VisualAge") OR  # CMake 2.6
         (CMAKE_Fortran_COMPILER_ID STREQUAL "XL"))           # CMake 2.8
    set(FOPT_ILP64 -qintsize=8)
  elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "NAG")
    set(FOPT_ILP64 -i8)
  elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC")
    if(WIN32)
      set(FOPT_ILP64 /i8)
    else()
      set(FOPT_ILP64 -i8)
    endif()
  else()
    set(CPE_ENV $ENV{PE_ENV})
    if(CPE_ENV STREQUAL "CRAY" AND NOT CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
      set(FOPT_ILP64 -sinteger64)
    elseif(CPE_ENV STREQUAL "NVIDIA" AND NOT CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
      set(FOPT_ILP64 -i8)
    else()
      set(FOPT_ILP64 -fdefault-integer-8)
    endif()
  endif()
  if(BUILD_INDEX64)
    target_compile_options(${target} INTERFACE "$<$<COMPILE_LANGUAGE:Fortran>:${FOPT_ILP64}>")
  endif()

  # GNU Fortran
  if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    set(FPE_EXIT_FLAG "-ffpe-trap=[izoupd]")

    target_compile_options(${target} INTERFACE "$<$<COMPILE_LANGUAGE:Fortran>:-frecursive>")

    if(CMAKE_Fortran_COMPILER_VERSION VERSION_LESS "8")
      target_compile_definitions(${target} INTERFACE "$<$<COMPILE_LANGUAGE:C>:FORTRAN_STRLEN=int>")
      set(FORTRAN_STRLEN_TYPE "int")
    endif()

    # Disabling loop vectorization for GNU Fortran versions affected by
    # https://gcc.gnu.org/bugzilla/show_bug.cgi?id=122408. See issue
    # https://github.com/Reference-LAPACK/lapack/issues/1160 as well.
    if(CMAKE_SYSTEM_PROCESSOR MATCHES "arm|arm64|aarch64" AND
       ((CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL "14.0" AND
         CMAKE_Fortran_COMPILER_VERSION VERSION_LESS "14.5") OR
        (CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL "15.0" AND
         CMAKE_Fortran_COMPILER_VERSION VERSION_LESS "15.3")))
      message(WARNING
        "Disabling loop vectorization for GNU Fortran (14.0-14.4, 15.0-15.2) on ARM "
        "due to a compiler bug (https://gcc.gnu.org/bugzilla/show_bug.cgi?id=122408). "
        "For full performance, consider changing to a different compiler or compiler version.")
      target_compile_options(${target} INTERFACE "$<$<COMPILE_LANGUAGE:Fortran>:-fno-tree-loop-vectorize>")
    endif()

    # Disabling SLP vectorization for GNU Fortran versions affected by
    # https://gcc.gnu.org/bugzilla/show_bug.cgi?id=124802.
    if(CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL "15.0" AND
       CMAKE_Fortran_COMPILER_VERSION VERSION_LESS "15.3")
      message(WARNING
        "Disabling SLP vectorization for GNU Fortran (15.0-15.2) "
        "due to a compiler bug (https://gcc.gnu.org/bugzilla/show_bug.cgi?id=124802). "
        "For full performance, consider changing to a different compiler or compiler version.")
      target_compile_options(${target} INTERFACE "$<$<COMPILE_LANGUAGE:Fortran>:-fno-tree-slp-vectorize>")
    endif()

  # Intel Fortran
  elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    set(FPE_EXIT_FLAG "[-/]fpe(-all=|)0")

    target_compile_options(${target} INTERFACE "$<$<COMPILE_LANGUAGE:Fortran>:-recursive>")
    if(WIN32)
      target_compile_options(${target} INTERFACE "$<$<COMPILE_LANGUAGE:Fortran>:/fp:strict>")
    else()
      target_compile_options(${target} INTERFACE "$<$<COMPILE_LANGUAGE:Fortran>:-fp-model=strict>")
    endif()

    # disable: The Intel(R) Fortran Compiler Classic (ifort) is deprecated
    if(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
      if(WIN32)
        target_compile_options(${target} INTERFACE "$<$<COMPILE_LANGUAGE:Fortran>:/Qdiag-disable:10448>")
        target_link_options(${target} INTERFACE "$<$<LINK_LANGUAGE:Fortran>:/Qdiag-disable:10448>")

        # Bad preprocessor line bogus warning
        target_compile_options(${target} INTERFACE "$<$<COMPILE_LANGUAGE:Fortran>:/Qdiag-disable:5117>")
      else()
        target_compile_options(${target} INTERFACE "$<$<COMPILE_LANGUAGE:Fortran>:-diag-disable:10448>")
        target_link_options(${target} INTERFACE "$<$<LINK_LANGUAGE:Fortran>:-diag-disable:10448>")
      endif()
    endif()

  # SunPro F95
  elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "SunPro")
    set(FPE_EXIT_FLAG "-ftrap=")
    set(FPE_DISABLE_FLAG "-ftrap=(%|)none")

    message(STATUS "Disabling FPE trap handlers with -ftrap=%none")
    target_compile_options(${target} INTERFACE "$<$<COMPILE_LANGUAGE:Fortran>:-ftrap=%none>")

    if(UNIX)
      # Delete libmtsk in linking sequence for Sun/Oracle Fortran Compiler.
      # This library is not present in the Sun package SolarisStudio12.3-linux-x86-bin
      list(REMOVE_ITEM CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES mtsk)
      # Compiler metadata is needed by CMake's mixed-language linker as well
      # as the explicit runtime requirements below. Keep this in LAPACK's
      # directory scope when embedded in another project.
      set(CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES
        "${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES}" PARENT_SCOPE)
    endif()

  # IBM XL Fortran
  elseif((CMAKE_Fortran_COMPILER_ID STREQUAL "VisualAge") OR  # CMake 2.6
         (CMAKE_Fortran_COMPILER_ID STREQUAL "XL"))           # CMake 2.8
    set(FPE_EXIT_FLAG "-qflttrap=[a-zA-Z:]:enable")

    target_compile_options(${target} INTERFACE "$<$<COMPILE_LANGUAGE:Fortran>:-qrecur>")
    if(UNIX)
      target_compile_options(${target} INTERFACE "$<$<COMPILE_LANGUAGE:Fortran>:-qnosave>")
      target_compile_options(${target} INTERFACE "$<$<COMPILE_LANGUAGE:Fortran>:-qstrict>")
    endif()

  # HP Fortran
  elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "HP")
    set(FPE_EXIT_FLAG "\\+fp_exception")

    message(STATUS "Enabling strict float conversion with +fltconst_strict")
    target_compile_options(${target} INTERFACE "$<$<COMPILE_LANGUAGE:Fortran>:+fltconst_strict>")

    # Most versions of cmake don't have good default options for the HP compiler
    target_compile_options(${target} INTERFACE "$<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:DEBUG>>:-g>")
    target_compile_options(${target} INTERFACE "$<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:MINSIZEREL>>:+Osize>")
    target_compile_options(${target} INTERFACE "$<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:RELEASE>>:+O2>")
    target_compile_options(${target} INTERFACE "$<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:RELWITHDEBINFO>>:+O2 -g>")

  # NAG Fortran
  elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "NAG")
    set(FPE_EXIT_FLAG "[-/]ieee=(stop|nonstd)")

    target_compile_options(${target} INTERFACE "$<$<COMPILE_LANGUAGE:Fortran>:-ieee=full>")
    target_compile_options(${target} INTERFACE "$<$<COMPILE_LANGUAGE:Fortran>:-dcfuns>")
    target_compile_options(${target} INTERFACE "$<$<COMPILE_LANGUAGE:Fortran>:-thread_safe>")
    target_link_options(${target} INTERFACE "$<$<LINK_LANGUAGE:Fortran>:-thread_safe>")
    target_compile_options(${target} INTERFACE "$<$<COMPILE_LANGUAGE:Fortran>:-recursive>")

    # By default NAG Fortran uses 32bit integers as hidden STRLEN arguments
    if(UNIX)
      if(APPLE)
        target_compile_definitions(${target} INTERFACE "$<$<COMPILE_LANGUAGE:C>:FORTRAN_STRLEN=int>")
        set(FORTRAN_STRLEN_TYPE "int")
      else()
        # Get all flags added via `target_compile_options(${target} INTERFACE ...)`
        get_target_property(COMP_OPTIONS ${target} INTERFACE_COMPILE_OPTIONS)

        if(NOT("${CMAKE_Fortran_FLAGS};${COMP_OPTIONS}" MATCHES "-abi=64c"))
          target_compile_definitions(${target} INTERFACE "$<$<COMPILE_LANGUAGE:C>:FORTRAN_STRLEN=int>")
          set(FORTRAN_STRLEN_TYPE "int")
        endif()
      endif()
    endif()

    # Disable warnings
    target_compile_options(${target} INTERFACE "$<$<COMPILE_LANGUAGE:Fortran>:-w=obs>")
    target_compile_options(${target} INTERFACE "$<$<COMPILE_LANGUAGE:Fortran>:-w=x77>")
    target_compile_options(${target} INTERFACE "$<$<COMPILE_LANGUAGE:Fortran>:-w=ques>")
    target_compile_options(${target} INTERFACE "$<$<COMPILE_LANGUAGE:Fortran>:-w=unused>")

    # Suppress compiler banner and summary
    include(CheckFortranCompilerFlag)
    check_fortran_compiler_flag("-quiet" _quiet)
    target_compile_options(${target} INTERFACE "$<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<BOOL:${_quiet}>>:-quiet>")
    target_link_options(${target} INTERFACE "$<$<AND:$<LINK_LANGUAGE:Fortran>,$<BOOL:${_quiet}>>:-quiet>")

  # NVIDIA HPC SDK
  elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC")
    set(FPE_EXIT_FLAG "-Ktrap=")
    set(FPE_DISABLE_FLAG "-Ktrap=none")

    target_compile_options(${target} INTERFACE "$<$<COMPILE_LANGUAGE:Fortran>:-Kieee>")
    target_compile_options(${target} INTERFACE "$<$<COMPILE_LANGUAGE:Fortran>:-Mrecursive>")

  # Flang Fortran
  elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "Flang")
    target_compile_options(${target} INTERFACE "$<$<COMPILE_LANGUAGE:Fortran>:-Mrecursive>")

  # LLVM Flang
  elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "LLVMFlang")
    # Nothing to do here for now, but this is a placeholder for future
    # LLVM Flang specific flags


  elseif(CMAKE_Fortran_COMPILER_LOADED)
    message(WARNING "Fortran local arrays should be allocated on the stack."
      " Please use a compiler which guarantees that feature."
      " See https://github.com/Reference-LAPACK/lapack/pull/188 and references therein.")
  endif()

  if(CMAKE_C_COMPILER_ID MATCHES "Intel")
    if(WIN32)
      target_compile_options(${target} INTERFACE "$<$<COMPILE_LANGUAGE:C>:/fp:strict>")
    else()
      target_compile_options(${target} INTERFACE "$<$<COMPILE_LANGUAGE:C>:-fp-model=strict>")
    endif()

    # disable: The Intel(R) C++ Compiler Classic (ICC) is deprecated
    if(CMAKE_C_COMPILER_ID STREQUAL "Intel")
      if(WIN32)
        target_compile_options(${target} INTERFACE "$<$<COMPILE_LANGUAGE:C>:/Qdiag-disable:10441>")
        target_link_options(${target} INTERFACE "$<$<LINK_LANGUAGE:C>:/Qdiag-disable:10441>")
      else()
        target_compile_options(${target} INTERFACE "$<$<COMPILE_LANGUAGE:C>:-diag-disable:10441>")
        target_link_options(${target} INTERFACE "$<$<LINK_LANGUAGE:C>:-diag-disable:10441>")
      endif()
    endif()
  endif()

  if(CMAKE_Fortran_FLAGS_RELEASE MATCHES "-O[3-9]|-Ofast")
    target_compile_options(${target} INTERFACE
      "$<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:Release>>:-O2>")
  endif()

  if(FPE_EXIT_FLAG AND "${CMAKE_Fortran_FLAGS};${CMAKE_Fortran_FLAGS_DEBUG};${CMAKE_Fortran_FLAGS_RELEASE};${CMAKE_Fortran_FLAGS_RELWITHDEBINFO};${CMAKE_Fortran_FLAGS_MINSIZEREL}" MATCHES "${FPE_EXIT_FLAG}"
      AND (NOT FPE_DISABLE_FLAG OR NOT CMAKE_Fortran_FLAGS MATCHES "${FPE_DISABLE_FLAG}"))
    message( FATAL_ERROR "Floating Point Exception (FPE) trap handlers are"
      " currently explicitly enabled in the compiler flags.  LAPACK is designed"
      " to check for and handle these cases internally and enabling these traps"
      " will likely cause LAPACK to crash.  Please re-configure with floating"
      " point exception trapping disabled.")
  endif()

  # A consumer need not enable Fortran or use the same compiler installation.
  # Resolve only the runtime libraries missing from the C toolchain's link line.
  set(runtime_libraries ${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES})
  list(REMOVE_ITEM runtime_libraries ${CMAKE_C_IMPLICIT_LINK_LIBRARIES})
  list(REMOVE_DUPLICATES runtime_libraries)
  set(lapack_fortran_runtime)
  foreach(runtime IN LISTS runtime_libraries)
    unset(runtime_path)
    find_library(runtime_path NAMES ${runtime}
      HINTS ${CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES} NO_CACHE)
    if(runtime_path)
      list(APPEND lapack_fortran_runtime "${runtime_path}")
    else()
      list(APPEND lapack_fortran_runtime "${runtime}")
    endif()
  endforeach()
  set(lapack_fortran_runtime "${lapack_fortran_runtime}" PARENT_SCOPE)
  set(FOPT_ILP64 "${FOPT_ILP64}" PARENT_SCOPE)
  set(FORTRAN_STRLEN_TYPE "${FORTRAN_STRLEN_TYPE}" PARENT_SCOPE)
endfunction()
