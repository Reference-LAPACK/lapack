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
  get_filename_component(destination "${target}_64_sources" ABSOLUTE BASE_DIR "${CMAKE_CURRENT_BINARY_DIR}")
  set(generated_source_files)

  foreach(source IN LISTS ${source_list})
    get_filename_component(source_abs "${source}" ABSOLUTE BASE_DIR "${CMAKE_CURRENT_SOURCE_DIR}")
    get_filename_component(source_name "${source}" NAME_WLE)
    get_filename_component(source_ext "${source}" EXT)
    set(output_file "${destination}/${source_name}_64${source_ext}")

    set(_fortran_extensions ".f" ".F" ".f90" ".F90")
    if(NOT source_ext IN_LIST _fortran_extensions)
      message(WARNING "Skipping non-Fortran source '${source}' for target '${target}'")
      continue()
    endif()

    message("${EXTENDED_API_GENERATOR}")

    add_custom_command(
      OUTPUT "${output_file}"
      COMMAND
        "${CMAKE_COMMAND}"
        "-DINPUT_FILE=${source_abs}"
        "-DOUTPUT_FILE=${output_file}"
        -P "${EXTENDED_API_GENERATOR}"
      DEPENDS
        "${source_abs}"
        "${EXTENDED_API_GENERATOR}"
      COMMENT "Generating 64-bit extended API source for ${source}"
      VERBATIM)

    list(APPEND generated_source_files "${output_file}")
  endforeach()

  if(generated_source_files)
    add_custom_target("${target}_64bit_codegen" ALL
      DEPENDS ${generated_source_files}
      COMMENT "Generating 64-bit suffixed sources for target ${target}")
    add_dependencies(64bit_codegen "${target}_64bit_codegen")
    set_source_files_properties(${generated_source_files} PROPERTIES GENERATED TRUE)
  endif()

  set(${generated_sources} ${generated_source_files} PARENT_SCOPE)
endfunction()
