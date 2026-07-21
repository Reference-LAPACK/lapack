include_guard(GLOBAL)

set(EXTENDED_API_GENERATOR
  "${CMAKE_CURRENT_LIST_DIR}/Generate64BitSuffixedSource.cmake" CACHE INTERNAL
  "Script that generates 64-bit suffixed sources for the extended API")

if(BUILD_INDEX64_EXT_API)
  add_custom_target(64bit_codegen ALL
    COMMENT "Generating 64-bit suffixed sources for extended API")
endif()

# Generate 64-bit suffixed sources for the extended API. The generation happens
# at build time in the Generate64BitSuffixedSource.cmake script.
function(generate_64bit_suffixed_sources target source_list generated_sources)
  set(options NO_STRING_REPLACEMENTS)
  cmake_parse_arguments(PARSE_ARGV 3 extended_api "${options}" "" "")

  get_filename_component(destination "${target}_64_sources" ABSOLUTE BASE_DIR "${CMAKE_CURRENT_BINARY_DIR}")
  get_property(_generated_64bit_source_files GLOBAL PROPERTY EXTENDED_API_GENERATED_SOURCE_FILES)
  set(new_generated_source_files)
  set(generated_source_files)

  foreach(_source IN LISTS ${source_list})
    get_filename_component(source_abs "${_source}" ABSOLUTE BASE_DIR "${CMAKE_CURRENT_SOURCE_DIR}")
    get_filename_component(source_name "${_source}" NAME_WLE)
    get_filename_component(source_ext "${_source}" EXT)
    set(output_file "${destination}/${source_name}_64${source_ext}")

    set(_fortran_extensions ".f" ".F" ".f90" ".F90")
    if(NOT source_ext IN_LIST _fortran_extensions)
      message(WARNING "Skipping non-Fortran source '${_source}' for target '${target}'")
      continue()
    endif()

    # Make sure we only have one custom command generating a given output file
    if(NOT "${output_file}" IN_LIST _generated_64bit_source_files)
      set(generator_args
        "-DINPUT_FILE=${source_abs}"
        "-DOUTPUT_FILE=${output_file}")
      if(extended_api_NO_STRING_REPLACEMENTS)
        list(APPEND generator_args "-DREPLACE_IN_STRINGS=OFF")
      endif()

      add_custom_command(
        OUTPUT "${output_file}"
        COMMAND
          "${CMAKE_COMMAND}"
          ${generator_args}
          -P "${EXTENDED_API_GENERATOR}"
        DEPENDS
          "${source_abs}"
          "${EXTENDED_API_GENERATOR}"
        COMMENT "Generating 64-bit extended API source for ${_source}"
        VERBATIM)

      list(APPEND new_generated_source_files "${output_file}")
    endif()

    list(APPEND generated_source_files "${output_file}")
  endforeach()

  # Make sure each generated source file is only part of one target to
  # avoid multiple targets trying to generate the same file
  if(new_generated_source_files)
    add_custom_target("${target}_64bit_codegen" ALL
      DEPENDS ${new_generated_source_files}
      COMMENT "Generating 64-bit suffixed sources for target ${target}")
    add_dependencies(64bit_codegen "${target}_64bit_codegen")

    set_property(GLOBAL APPEND PROPERTY EXTENDED_API_GENERATED_SOURCE_FILES ${new_generated_source_files})
  endif()

  set(${generated_sources} ${generated_source_files} PARENT_SCOPE)
endfunction()
