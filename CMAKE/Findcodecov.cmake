# This file is part of CMake-codecov.
#
# https://github.com/RWTH-ELP/CMake-codecov
#
# Copyright (c)
#   2015-2016 RWTH Aachen University, Federal Republic of Germany
#
# LICENSE : BSD 3-Clause License
#
# Written by Alexander Haase, alexander.haase@rwth-aachen.de
# Updated by Guillaume Jacquenot, guillaume.jacquenot@gmail.com

include(CheckCCompilerFlag)
include(CheckCXXCompilerFlag)
include(CheckFortranCompilerFlag)

set(COVERAGE_FLAG_CANDIDATES
  # gcc and clang
  "-O0 -g -fprofile-arcs -ftest-coverage"

  # gcc and clang fallback
  "-O0 -g --coverage"
)

# The languages this module knows how to probe for coverage flags.
set(COVERAGE_LANGUAGES C CXX Fortran)


# Find the coverage flags for language ${LANG} and cache them in
# COVERAGE_<compiler id>_FLAGS. Coverage flags are not dependent on the
# language, but on the compiler, so languages sharing a compiler are probed
# only once. A language may be enabled after this module has been found, so
# this is called again from add_coverage_target() rather than only here.
function(codecov_find_flags LANG)
  set(COMPILER ${CMAKE_${LANG}_COMPILER_ID})
  if(NOT COMPILER OR DEFINED CACHE{COVERAGE_${COMPILER}_FLAGS})
    return()
  endif()

  set(CMAKE_REQUIRED_QUIET ${codecov_FIND_QUIETLY})

  foreach(FLAG IN LISTS COVERAGE_FLAG_CANDIDATES)
    if(NOT CMAKE_REQUIRED_QUIET)
      message(STATUS "Try ${COMPILER} code coverage flag = [${FLAG}]")
    endif()

    set(CMAKE_REQUIRED_FLAGS "${FLAG}")
    unset(COVERAGE_FLAG_DETECTED CACHE)

    if(LANG STREQUAL "C")
      check_c_compiler_flag("${FLAG}" COVERAGE_FLAG_DETECTED)
    elseif(LANG STREQUAL "CXX")
      check_cxx_compiler_flag("${FLAG}" COVERAGE_FLAG_DETECTED)
    elseif(LANG STREQUAL "Fortran")
      check_fortran_compiler_flag("${FLAG}" COVERAGE_FLAG_DETECTED)
    endif()

    if(COVERAGE_FLAG_DETECTED)
      # Cache the flags as a list, so they can be handed to
      # target_compile_options() and target_link_options() unmodified.
      string(REPLACE " " ";" FLAG_LIST "${FLAG}")
      set(COVERAGE_${COMPILER}_FLAGS "${FLAG_LIST}"
        CACHE STRING "${COMPILER} flags for code coverage.")
      mark_as_advanced(COVERAGE_${COMPILER}_FLAGS)
      return()
    endif()
  endforeach()

  # Remember that this compiler has no usable coverage flags, so that it is not
  # probed again for every target.
  set(COVERAGE_${COMPILER}_FLAGS "COVERAGE_${COMPILER}_FLAGS-NOTFOUND"
    CACHE STRING "${COMPILER} flags for code coverage.")
  mark_as_advanced(COVERAGE_${COMPILER}_FLAGS)
endfunction()


# Probe the languages that are already enabled, so that the check appears in
# the configure output next to the other compiler checks.
get_property(ENABLED_LANGUAGES GLOBAL PROPERTY ENABLED_LANGUAGES)
foreach(LANG IN LISTS ENABLED_LANGUAGES)
  if(LANG IN_LIST COVERAGE_LANGUAGES)
    codecov_find_flags(${LANG})
  endif()
endforeach()


# Helper function to get the language of a source file.
function(codecov_lang_of_source FILE RETURN_VAR)
  get_filename_component(FILE_EXT "${FILE}" EXT)
  string(TOLOWER "${FILE_EXT}" FILE_EXT)
  if(FILE_EXT STREQUAL "")
    set(${RETURN_VAR} "" PARENT_SCOPE)
    return()
  endif()
  string(SUBSTRING "${FILE_EXT}" 1 -1 FILE_EXT)

  get_property(ENABLED_LANGUAGES GLOBAL PROPERTY ENABLED_LANGUAGES)
  foreach(LANG IN LISTS ENABLED_LANGUAGES)
    if(FILE_EXT IN_LIST CMAKE_${LANG}_SOURCE_FILE_EXTENSIONS)
      set(${RETURN_VAR} "${LANG}" PARENT_SCOPE)
      return()
    endif()
  endforeach()

  set(${RETURN_VAR} "" PARENT_SCOPE)
endfunction()


# Helper function to get the relative path of the source file destination path.
# This path is needed by FindGcov and FindLcov cmake files to locate the
# captured data.
function(codecov_path_of_source FILE RETURN_VAR)
  # Generator expressions cannot be resolved to a path here, so they have no
  # coverage data of their own. In particular a $<TARGET_OBJECTS:...> source
  # belongs to an object library and is evaluated together with that library.
  # Note that a conditional such as
  #   $<$<BOOL:${COND}>: $<TARGET_OBJECTS:foo>>
  # reaches us split at the whitespace, so only one of its fragments mentions
  # TARGET_OBJECTS; skip anything that looks like a generator expression.
  if(FILE MATCHES "\\$<")
    set(${RETURN_VAR} "" PARENT_SCOPE)
    return()
  endif()

  string(REPLACE "${CMAKE_CURRENT_BINARY_DIR}/" "" FILE "${FILE}")
  if(IS_ABSOLUTE ${FILE})
    file(RELATIVE_PATH FILE ${CMAKE_CURRENT_SOURCE_DIR} ${FILE})
  endif()

  # get the right path for file
  string(REPLACE ".." "__" PATH "${FILE}")

  set(${RETURN_VAR} "${PATH}" PARENT_SCOPE)
endfunction()

# Add coverage support for target ${TNAME} and register target for coverage
# evaluation.
function(add_coverage_target TNAME)
  # A target may mix languages that are built by different compilers, which may
  # use incompatible coverage implementations. Guard each compiler's flags with
  # the language they were detected for, so that every source is compiled with
  # the flags of the compiler that actually builds it.
  set(COMPILE_OPTIONS "")
  set(LINK_OPTIONS "")
  foreach(LANG IN LISTS COVERAGE_LANGUAGES)
    codecov_find_flags(${LANG})

    set(FLAGS "${COVERAGE_${CMAKE_${LANG}_COMPILER_ID}_FLAGS}")
    if(FLAGS)
      list(APPEND COMPILE_OPTIONS "$<$<COMPILE_LANGUAGE:${LANG}>:${FLAGS}>")
      list(APPEND LINK_OPTIONS "$<$<LINK_LANGUAGE:${LANG}>:${FLAGS}>")
    endif()
  endforeach()

  if(NOT COMPILE_OPTIONS)
    message(AUTHOR_WARNING "Coverage disabled for target ${TNAME} because no "
      "coverage flags are available for any enabled language.")
    return()
  endif()

  # enable coverage for target
  target_compile_options(${TNAME} PRIVATE ${COMPILE_OPTIONS})
  target_link_options(${TNAME} PRIVATE ${LINK_OPTIONS})


  # Add gcov files generated by compiler to clean target.
  get_target_property(TSOURCES ${TNAME} SOURCES)
  set(CLEAN_FILES "")
  foreach(FILE IN LISTS TSOURCES)
    codecov_path_of_source("${FILE}" FILE)
    if(FILE STREQUAL "")
      continue()
    endif()

    codecov_lang_of_source("${FILE}" LANG)
    if(LANG)
      list(APPEND CLEAN_FILES
        "CMakeFiles/${TNAME}.dir/${FILE}.gcno"
        "CMakeFiles/${TNAME}.dir/${FILE}.gcda")
    endif()
  endforeach()

  set_property(DIRECTORY APPEND PROPERTY ADDITIONAL_CLEAN_FILES ${CLEAN_FILES})

  add_gcov_target(${TNAME})
endfunction()


# Add coverage support for target ${TNAME} and register target for coverage
# evaluation.
function(add_coverage TNAME)
  foreach(T IN LISTS ARGV)
    add_coverage_target(${T})
  endforeach()
endfunction()


# Include modules for parsing the collected data and output it in a readable
# format (like gcov).
find_package(Gcov)
