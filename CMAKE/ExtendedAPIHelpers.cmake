include_guard(GLOBAL)

set(EXTENDED_API_GENERATOR
  "${CMAKE_CURRENT_LIST_DIR}/GenerateSuffixedSource.cmake" CACHE INTERNAL
  "Script that generates suffixed sources for extended APIs")

if(BUILD_INDEX64_EXT_API)
  add_custom_target(64bit_codegen ALL
    COMMENT "Generating 64-bit suffixed sources for extended API")
endif()

# Generate suffixed sources for an extended API. The generation happens at
# build time in the GenerateSuffixedSource.cmake script. Arguments:
#   SUFFIX <suffix>            -- suffix appended to symbols and file names
#                                 (required)
#   SYMBOL_ALLOWLIST <sym>...  -- only rename the listed symbols
#                                 (case-insensitive); by default every symbol
#                                 defined in the sources is renamed
#   NO_STRING_REPLACEMENTS     -- do not rename symbols inside string literals
function(generate_suffixed_sources target source_list generated_sources)
  set(options NO_STRING_REPLACEMENTS)
  set(oneValueArgs SUFFIX)
  set(multiValueArgs SYMBOL_ALLOWLIST)
  cmake_parse_arguments(PARSE_ARGV 3 extended_api "${options}" "${oneValueArgs}" "${multiValueArgs}")

  if(NOT DEFINED extended_api_SUFFIX)
    message(FATAL_ERROR "generate_suffixed_sources: SUFFIX is required")
  endif()

  get_filename_component(destination "${target}${extended_api_SUFFIX}_sources" ABSOLUTE BASE_DIR "${CMAKE_CURRENT_BINARY_DIR}")
  get_property(_generated_suffixed_source_files GLOBAL PROPERTY EXTENDED_API_GENERATED_SOURCE_FILES)
  set(new_generated_source_files)
  set(generated_source_files)
  set(_pending_codegen_inputs)
  set(_pending_codegen_outputs)

  foreach(_source IN LISTS ${source_list})
    get_filename_component(source_abs "${_source}" ABSOLUTE BASE_DIR "${CMAKE_CURRENT_SOURCE_DIR}")
    get_filename_component(source_name "${_source}" NAME_WLE)
    get_filename_component(source_ext "${_source}" EXT)
    set(output_file "${destination}/${source_name}${extended_api_SUFFIX}${source_ext}")

    set(_fortran_extensions ".f" ".F" ".f90" ".F90")
    if(NOT source_ext IN_LIST _fortran_extensions)
      message(WARNING "Skipping non-Fortran source '${_source}' for target '${target}'")
      continue()
    endif()

    # Make sure we only have one custom command generating a given output file.
    # Check both the global registry (outputs queued by earlier calls) and
    # the outputs already queued by this call: source lists may repeat a
    # file (e.g. dchklq.f appears twice in DLINTST), which would otherwise
    # put the same OUTPUT in two different unity batches -> Ninja
    # "multiple rules generate ..." error (identical per-file rules are
    # tolerated, differing batch rules are not).
    if(NOT "${output_file}" IN_LIST _generated_suffixed_source_files AND NOT "${output_file}" IN_LIST new_generated_source_files)
      list(APPEND _pending_codegen_inputs "${source_abs}")
      list(APPEND _pending_codegen_outputs "${output_file}")
      list(APPEND new_generated_source_files "${output_file}")
    endif()

    list(APPEND generated_source_files "${output_file}")
  endforeach()

  # With CMAKE_UNITY_BUILD the many tiny per-file codegen rules become the
  # bottleneck (one cmake process per file), so group generation into
  # unity-like batches of CMAKE_UNITY_BUILD_BATCH_SIZE files sharing a
  # single cmake invocation.
  set(_use_unity_batches FALSE)
  if(CMAKE_UNITY_BUILD)
    set(_use_unity_batches TRUE)
  elseif(TARGET "${target}")
    get_target_property(_target_unity_build "${target}" UNITY_BUILD)
    if(_target_unity_build)
      set(_use_unity_batches TRUE)
    endif()
  endif()
  if(_use_unity_batches)
    if(DEFINED CMAKE_UNITY_BUILD_BATCH_SIZE AND CMAKE_UNITY_BUILD_BATCH_SIZE GREATER 0)
      set(_codegen_batch_size "${CMAKE_UNITY_BUILD_BATCH_SIZE}")
    else()
      set(_codegen_batch_size 8)
    endif()
  endif()

  # Options shared by every file of this call (hence by every batch).
  set(_common_generator_args "-DSUFFIX=${extended_api_SUFFIX}")
  if(extended_api_NO_STRING_REPLACEMENTS)
    list(APPEND _common_generator_args "-DREPLACE_IN_STRINGS=OFF")
  endif()
  if(DEFINED extended_api_SYMBOL_ALLOWLIST)
    # Join with $<SEMICOLON> so that the allowlist stays one -D argument
    # on the command line but still reaches the script as a CMake list.
    # An allowlist change re-runs the generation (the command changes).
    list(JOIN extended_api_SYMBOL_ALLOWLIST "$<SEMICOLON>" _symbol_allowlist)
    list(APPEND _common_generator_args "-DSYMBOL_ALLOWLIST=${_symbol_allowlist}")
  endif()

  if(_use_unity_batches AND _pending_codegen_inputs)
    list(LENGTH _pending_codegen_inputs _pending_count)
    set(_batch_index 0)
    set(_batch_offset 0)
    while(_batch_offset LESS _pending_count)
      set(_batch_inputs)
      set(_batch_outputs)
      set(_batch_count 0)
      while(_batch_count LESS _codegen_batch_size)
        math(EXPR _pending_index "${_batch_offset} + ${_batch_count}")
        if(NOT _pending_index LESS _pending_count)
          break()
        endif()
        list(GET _pending_codegen_inputs ${_pending_index} _batch_input)
        list(GET _pending_codegen_outputs ${_pending_index} _batch_output)
        list(APPEND _batch_inputs "${_batch_input}")
        list(APPEND _batch_outputs "${_batch_output}")
        math(EXPR _batch_count "${_batch_count} + 1")
      endwhile()
      list(JOIN _batch_inputs "$<SEMICOLON>" _batch_inputs_joined)
      list(JOIN _batch_outputs "$<SEMICOLON>" _batch_outputs_joined)
      list(LENGTH _batch_outputs _batch_len)
      add_custom_command(
        OUTPUT ${_batch_outputs}
        COMMAND
          "${CMAKE_COMMAND}"
          "-DINPUT_FILES=${_batch_inputs_joined}"
          "-DOUTPUT_FILES=${_batch_outputs_joined}"
          ${_common_generator_args}
          -P "${EXTENDED_API_GENERATOR}"
        DEPENDS
          ${_batch_inputs}
          "${EXTENDED_API_GENERATOR}"
        COMMENT "Generating ${extended_api_SUFFIX} extended API Unity/unity_${_batch_index} (${_batch_len} sources for target ${target})"
        VERBATIM)
      math(EXPR _batch_offset "${_batch_offset} + ${_batch_count}")
      math(EXPR _batch_index "${_batch_index} + 1")
    endwhile()
  else()
    list(LENGTH _pending_codegen_inputs _pending_count)
    if(_pending_count GREATER 0)
      math(EXPR _last_pending_index "${_pending_count} - 1")
      foreach(_pending_index RANGE 0 ${_last_pending_index})
        list(GET _pending_codegen_inputs ${_pending_index} _single_input)
        list(GET _pending_codegen_outputs ${_pending_index} _single_output)
        get_filename_component(_single_display "${_single_input}" NAME)
        add_custom_command(
          OUTPUT "${_single_output}"
          COMMAND
            "${CMAKE_COMMAND}"
            "-DINPUT_FILE=${_single_input}"
            "-DOUTPUT_FILE=${_single_output}"
            ${_common_generator_args}
            -P "${EXTENDED_API_GENERATOR}"
          DEPENDS
            "${_single_input}"
            "${EXTENDED_API_GENERATOR}"
          COMMENT "Generating ${extended_api_SUFFIX} extended API source for ${_single_display}"
          VERBATIM)
      endforeach()
    endif()
  endif()

  # Make sure each generated source file is only part of one target to
  # avoid multiple targets trying to generate the same file
  if(new_generated_source_files)
    add_custom_target("${target}${extended_api_SUFFIX}_codegen" ALL
      DEPENDS ${new_generated_source_files}
      COMMENT "Generating ${extended_api_SUFFIX} suffixed sources for target ${target}")
    if(extended_api_SUFFIX STREQUAL "_64" AND TARGET 64bit_codegen)
      add_dependencies(64bit_codegen "${target}${extended_api_SUFFIX}_codegen")
    endif()

    set_property(GLOBAL APPEND PROPERTY EXTENDED_API_GENERATED_SOURCE_FILES ${new_generated_source_files})
  endif()

  set(${generated_sources} ${generated_source_files} PARENT_SCOPE)
endfunction()

# Generate 64-bit suffixed sources for the extended API. Kept as a thin
# wrapper around generate_suffixed_sources for the Index-64 extended API.
function(generate_64bit_suffixed_sources target source_list generated_sources)
  set(options NO_STRING_REPLACEMENTS)
  cmake_parse_arguments(PARSE_ARGV 3 extended_api "${options}" "" "")

  set(_forwarded_options)
  if(extended_api_NO_STRING_REPLACEMENTS)
    list(APPEND _forwarded_options NO_STRING_REPLACEMENTS)
  endif()

  generate_suffixed_sources("${target}" "${source_list}" _generated_source_files
    SUFFIX "_64" ${_forwarded_options})

  set(${generated_sources} ${_generated_source_files} PARENT_SCOPE)
endfunction()
