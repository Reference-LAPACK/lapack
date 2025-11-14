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

macro(CheckLAPACKCompilerFlags)

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
    if(WIN32)
      set(FOPT_ILP64 /i8)
    else()
      set(FOPT_ILP64 -i8)
    endif()
  elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC")
    if(WIN32)
      set(FOPT_ILP64 /i8)
    else()
      set(FOPT_ILP64 -i8)
    endif()
  else()
    set(CPE_ENV $ENV{PE_ENV})
    if(CPE_ENV STREQUAL "CRAY")
      set(FOPT_ILP64 -sinteger64)
    elseif(CPE_ENV STREQUAL "NVIDIA")
      set(FOPT_ILP64 -i8)
    else()
      set(FOPT_ILP64 -fdefault-integer-8)
    endif()
  endif()
  if(FORTRAN_ILP)
    add_compile_options("$<$<COMPILE_LANGUAGE:Fortran>:${FOPT_ILP64}>")
  endif()

  # GNU Fortran
  if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    set(FPE_EXIT_FLAG "-ffpe-trap=[izoupd]")

    add_compile_options("$<$<COMPILE_LANGUAGE:Fortran>:-frecursive>")

    if(CMAKE_Fortran_COMPILER_VERSION VERSION_LESS "8")
      add_compile_definitions("$<$<COMPILE_LANGUAGE:C>:FORTRAN_STRLEN=int>")
    endif()

    # Disabling loop vectorization for GNU Fortran versions affected by
    # https://gcc.gnu.org/bugzilla/show_bug.cgi?id=122408. See issue
    # https://github.com/Reference-LAPACK/lapack/issues/1160 as well.
    if(CMAKE_HOST_SYSTEM_PROCESSOR MATCHES "arm|arm64|aarch64")
      if((CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL "14.0" AND
          CMAKE_Fortran_COMPILER_VERSION VERSION_LESS_EQUAL "14.4") OR
         (CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL "15.0" AND
          CMAKE_Fortran_COMPILER_VERSION VERSION_LESS_EQUAL "15.2"))
        message(WARNING
          "Disabling loop vectorization for GNU Fortran (14.0-14.4, 15.0-15.2) on ARM "
          "due to a compiler bug (https://gcc.gnu.org/bugzilla/show_bug.cgi?id=122408). "
          "For full performance, consider changing to a different compiler or compiler version.")
        add_compile_options("$<$<COMPILE_LANGUAGE:Fortran>:-fno-tree-loop-vectorize>")
      endif()
    endif()

  # Intel Fortran
  elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    set(FPE_EXIT_FLAG "[-/]fpe(-all=|)0")

    add_compile_options("$<$<COMPILE_LANGUAGE:Fortran>:-recursive>")
    if(UNIX)
      add_compile_options("$<$<COMPILE_LANGUAGE:Fortran>:SHELL:-fp-model strict>")
    endif()

  # SunPro F95
  elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "SunPro")
    set(FPE_EXIT_FLAG "-ftrap=")
    set(FPE_DISABLE_FLAG "-ftrap=(%|)none")

    message(STATUS "Disabling FPE trap handlers with -ftrap=%none")
    add_compile_options("$<$<COMPILE_LANGUAGE:Fortran>:-ftrap=%none>")

    if(UNIX)
      # Delete libmtsk in linking sequence for Sun/Oracle Fortran Compiler.
      # This library is not present in the Sun package SolarisStudio12.3-linux-x86-bin
      string(REPLACE \;mtsk\; \; CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES
             "${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES}")
    endif()

  # IBM XL Fortran
  elseif((CMAKE_Fortran_COMPILER_ID STREQUAL "VisualAge") OR  # CMake 2.6
         (CMAKE_Fortran_COMPILER_ID STREQUAL "XL"))           # CMake 2.8
    set(FPE_EXIT_FLAG "-qflttrap=[a-zA-Z:]:enable")

    add_compile_options("$<$<COMPILE_LANGUAGE:Fortran>:-qrecur>")
    if(UNIX)
      add_compile_options("$<$<COMPILE_LANGUAGE:Fortran>:-qnosave>")
      add_compile_options("$<$<COMPILE_LANGUAGE:Fortran>:-qstrict>")
    endif()

  # HP Fortran
  elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "HP")
    set(FPE_EXIT_FLAG "\\+fp_exception")

    message(STATUS "Enabling strict float conversion with +fltconst_strict")
    add_compile_options("$<$<COMPILE_LANGUAGE:Fortran>:+fltconst_strict>")

    # Most versions of cmake don't have good default options for the HP compiler
    add_compile_options("$<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:DEBUG>>:-g>")
    add_compile_options("$<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:MINSIZEREL>>:+Osize>")
    add_compile_options("$<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:RELEASE>>:+O2>")
    add_compile_options("$<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:RELWITHDEBINFO>>:+O2 -g>")

  # NAG Fortran
  elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "NAG")
    set(FPE_EXIT_FLAG "[-/]ieee=(stop|nonstd)")

    add_compile_options("$<$<COMPILE_LANGUAGE:Fortran>:-ieee=full>")
    add_compile_options("$<$<COMPILE_LANGUAGE:Fortran>:-dcfuns>")
    add_compile_options("$<$<COMPILE_LANGUAGE:Fortran>:-thread_safe>")
    add_link_options("$<$<COMPILE_LANGUAGE:Fortran>:-thread_safe>")
    add_compile_options("$<$<COMPILE_LANGUAGE:Fortran>:-recursive>")

    # By default NAG Fortran uses 32bit integers as hidden STRLEN arguments
    if(UNIX)
      if(APPLE)
        add_compile_definitions("$<$<COMPILE_LANGUAGE:C>:FORTRAN_STRLEN=int>")
      else()
        # Get all flags added via `add_compile_options(...)`
        get_directory_property(COMP_OPTIONS COMPILE_OPTIONS)

        if(NOT("${CMAKE_Fortran_FLAGS};${COMP_OPTIONS}" MATCHES "-abi=64c"))
          add_compile_definitions("$<$<COMPILE_LANGUAGE:C>:FORTRAN_STRLEN=int>")
        endif()
      endif()
    endif()

    # Disable warnings
    add_compile_options("$<$<COMPILE_LANGUAGE:Fortran>:-w=obs>")
    add_compile_options("$<$<COMPILE_LANGUAGE:Fortran>:-w=x77>")
    add_compile_options("$<$<COMPILE_LANGUAGE:Fortran>:-w=ques>")
    add_compile_options("$<$<COMPILE_LANGUAGE:Fortran>:-w=unused>")

    # Suppress compiler banner and summary
    include(CheckFortranCompilerFlag)
    check_fortran_compiler_flag("-quiet" _quiet)
    add_compile_options("$<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<BOOL:${_quiet}>>:-quiet>")
    add_link_options("$<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<BOOL:${_quiet}>>:-quiet>")

  # NVIDIA HPC SDK
  elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC")
    set(FPE_EXIT_FLAG "-Ktrap=")
    set(FPE_DISABLE_FLAG "-Ktrap=none")

    add_compile_options("$<$<COMPILE_LANGUAGE:Fortran>:-Kieee>")
    add_compile_options("$<$<COMPILE_LANGUAGE:Fortran>:-Mrecursive>")

  # Flang Fortran
  elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "Flang")
    add_compile_options("$<$<COMPILE_LANGUAGE:Fortran>:-Mrecursive>")

  # Compaq Fortran
  elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "Compaq")
    if(WIN32)
      if(CMAKE_GENERATOR STREQUAL "NMake Makefiles")
        get_filename_component(CMAKE_Fortran_COMPILER_CMDNAM ${CMAKE_Fortran_COMPILER} NAME_WE)
        message(STATUS "Using Compaq Fortran compiler with command name ${CMAKE_Fortran_COMPILER_CMDNAM}")
        set(cmd ${CMAKE_Fortran_COMPILER_CMDNAM})
        string(TOLOWER "${cmd}" cmdlc)
        if(cmdlc STREQUAL "df")
          message(STATUS "Assume the Compaq Visual Fortran Compiler is being used")
          set(CMAKE_Fortran_USE_RESPONSE_FILE_FOR_OBJECTS 1)
          set(CMAKE_Fortran_USE_RESPONSE_FILE_FOR_INCLUDES 1)
          #This is a workaround that is needed to avoid forward-slashes in the
          #filenames listed in response files from incorrectly being interpreted as
          #introducing compiler command options
          if(${BUILD_SHARED_LIBS})
            message(FATAL_ERROR "Making of shared libraries with CVF has not been tested.")
          endif()
          set(str "NMake version 9 or later should be used. NMake version 6.0 which is\n")
          set(str "${str}   included with the CVF distribution fails to build Lapack because\n")
          set(str "${str}   the number of source files exceeds the limit for NMake v6.0\n")
          message(STATUS ${str})
          set(CMAKE_Fortran_LINK_EXECUTABLE "LINK /out:<TARGET> <LINK_FLAGS> <LINK_LIBRARIES> <OBJECTS>")
        endif()
      endif()
    endif()

  else()
    message(WARNING "Fortran local arrays should be allocated on the stack."
      " Please use a compiler which guarantees that feature."
      " See https://github.com/Reference-LAPACK/lapack/pull/188 and references therein.")
  endif()

  if("${CMAKE_Fortran_FLAGS_RELEASE}" MATCHES "O[3-9]")
    message(STATUS "Reducing RELEASE optimization level to O2")
    string(REGEX REPLACE "O[3-9]" "O2" CMAKE_Fortran_FLAGS_RELEASE
           "${CMAKE_Fortran_FLAGS_RELEASE}")
  endif()

  # Get all flags added via `add_compile_options(...)`
  get_directory_property(COMP_OPTIONS COMPILE_OPTIONS)

  if(("${CMAKE_Fortran_FLAGS};${COMP_OPTIONS}" MATCHES "${FPE_EXIT_FLAG}") AND NOT
     ("${CMAKE_Fortran_FLAGS};${COMP_OPTIONS}" MATCHES "${FPE_DISABLE_FLAG}"))
    message( FATAL_ERROR "Floating Point Exception (FPE) trap handlers are"
      " currently explicitly enabled in the compiler flags.  LAPACK is designed"
      " to check for and handle these cases internally and enabling these traps"
      " will likely cause LAPACK to crash.  Please re-configure with floating"
      " point exception trapping disabled.")
  endif()

endmacro()
