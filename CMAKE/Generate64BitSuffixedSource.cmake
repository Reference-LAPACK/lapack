if(NOT DEFINED INPUT_FILE)
  message(FATAL_ERROR "INPUT_FILE must be set")
endif()

if(NOT DEFINED OUTPUT_FILE)
  message(FATAL_ERROR "OUTPUT_FILE must be set")
endif()

# Check whether the input file is fixed or free form based on its extension.
get_filename_component(input_extension "${INPUT_FILE}" LAST_EXT)
string(TOLOWER "${input_extension}" input_extension_lower)
if(input_extension_lower MATCHES "^\\.f(90|95|03|08)$")
  set(_rewrite_source_is_free_form TRUE)
else()
  set(_rewrite_source_is_free_form FALSE)
endif()

# Check if line is a comment line based on the form of the Fortran source.
function(_is_fortran_comment_line line result)
  if(_rewrite_source_is_free_form)
    if("${line}" MATCHES "^[ \t]*!")
      set(${result} TRUE PARENT_SCOPE)
    else()
      set(${result} FALSE PARENT_SCOPE)
    endif()
  elseif("${line}" MATCHES "^[cC\\*!]" OR "${line}" MATCHES "^[ \t]*!")
    set(${result} TRUE PARENT_SCOPE)
  else()
    set(${result} FALSE PARENT_SCOPE)
  endif()
endfunction()

# Check if the line is a fixed-form continuation line.
function(_is_fixed_form_continuation_line line result)
  if(_rewrite_source_is_free_form)
    set(${result} FALSE PARENT_SCOPE)
    return()
  endif()

  string(LENGTH "${line}" line_length)
  if(line_length GREATER 5)
    string(SUBSTRING "${line}" 5 1 continuation_char)
    if(NOT continuation_char STREQUAL " " AND NOT continuation_char STREQUAL "0")
      set(${result} TRUE PARENT_SCOPE)
      return()
    endif()
  endif()

  set(${result} FALSE PARENT_SCOPE)
endfunction()

# Check if the line ends with a free-form continuation character.
function(_has_trailing_free_form_continuation line result)
  if("${line}" MATCHES "&[ \t]*$")
    set(${result} TRUE PARENT_SCOPE)
  else()
    set(${result} FALSE PARENT_SCOPE)
  endif()
endfunction()

# Remove leading and trailing continuation characters and whitespace
function(_normalize_fortran_statement_fragment line result)
  if(_rewrite_source_is_free_form)
    string(REGEX REPLACE "^[ \t]*&[ \t]*" "" fragment "${line}")
    string(REGEX REPLACE "[ \t]*&[ \t]*$" "" fragment "${fragment}")
  else()
    string(LENGTH "${line}" line_length)
    if(line_length GREATER 5)
      string(SUBSTRING "${line}" 6 -1 fragment)
    endif()
  endif()

  set(${result} "${fragment}" PARENT_SCOPE)
endfunction()

# Create one logical Fortran statement by appending all continuation lines.
function(_collect_logical_fortran_line line_list_name start_index logical_line next_index)
  list(GET ${line_list_name} ${start_index} current_line)
  _normalize_fortran_statement_fragment("${current_line}" merged_line)

  list(LENGTH ${line_list_name} raw_line_count)
  set(last_line "${current_line}")
  math(EXPR current_index "${start_index} + 1")

  while(current_index LESS raw_line_count)
    list(GET ${line_list_name} ${current_index} continuation_line)
    _is_fortran_comment_line("${continuation_line}" continuation_is_comment)
    if(continuation_is_comment)
      break()
    endif()

    # Check for continuation based on the form of the Fortran source
    if(_rewrite_source_is_free_form)
      if("${last_line}" MATCHES "&[ \t]*$" OR "${continuation_line}" MATCHES "^[ \t]*&")
        _normalize_fortran_statement_fragment("${continuation_line}" continuation_fragment)
      else()
        break()
      endif()
    else()
      _is_fixed_form_continuation_line("${continuation_line}" is_fixed_form_continuation)
      if(is_fixed_form_continuation)
        _normalize_fortran_statement_fragment("${continuation_line}" continuation_fragment)
      else()
        break()
      endif()
    endif()

    string(APPEND merged_line " ${continuation_fragment}")
    set(last_line "${continuation_line}")
    math(EXPR current_index "${current_index} + 1")
  endwhile()

  set(${logical_line} "${merged_line}" PARENT_SCOPE)
  set(${next_index} ${current_index} PARENT_SCOPE)
endfunction()

# Extract function/subroutine names from one logical Fortran statement.
function(_extract_symbols statement_text result)
  set(symbols)

  if("${statement_text}" MATCHES "(external|EXTERNAL)")
    string(REGEX REPLACE
      "^[a-zA-Z0-9_ *]*(external|EXTERNAL)[ ]*" ""
      external_names "${statement_text}")
    string(REPLACE "," ";" external_names "${external_names}")
    foreach(external_name IN LISTS external_names)
      string(STRIP "${external_name}" external_name)
      if(NOT external_name STREQUAL "")
        list(APPEND symbols "${external_name}")
      endif()
    endforeach()
  elseif("${statement_text}" MATCHES "(subroutine|SUBROUTINE|function|FUNCTION)")
    string(REGEX REPLACE
      "^[a-zA-Z0-9_ *]*(subroutine|SUBROUTINE|function|FUNCTION)[ ]*" ""
      symbol_name "${statement_text}")
    string(REGEX REPLACE "[(].*$" "" symbol_name "${symbol_name}")
    string(STRIP "${symbol_name}" symbol_name)
    if(NOT symbol_name STREQUAL "")
      list(APPEND symbols "${symbol_name}")
    endif()
  endif()

  set(${result} ${symbols} PARENT_SCOPE)
endfunction()

get_filename_component(output_dir "${OUTPUT_FILE}" DIRECTORY)
file(MAKE_DIRECTORY "${output_dir}")

file(READ "${INPUT_FILE}" source_content)
set(rewritten_content "${source_content}")

# Extract symbol names from the source file
file(STRINGS "${INPUT_FILE}" raw_lines)
set(symbol_names)
list(LENGTH raw_lines raw_line_count)
set(line_index 0)
while(line_index LESS raw_line_count)
  list(GET raw_lines ${line_index} current_line)
  _is_fortran_comment_line("${current_line}" is_comment_line)
  if(is_comment_line)
    math(EXPR line_index "${line_index} + 1")
    continue()
  endif()

  _collect_logical_fortran_line(raw_lines ${line_index} logical_line next_index)
  _extract_symbols("${logical_line}" symbols)
  if(symbols)
    list(APPEND symbol_names ${symbols})
  endif()
  set(line_index ${next_index})
endwhile()
list(REMOVE_DUPLICATES symbol_names)

# Replace symbol names with their suffixed versions in the source content
foreach(symbol_name IN LISTS symbol_names)
  set(symbol_variants "${symbol_name}")
  string(TOLOWER "${symbol_name}" symbol_lower)
  string(TOUPPER "${symbol_name}" symbol_upper)
  foreach(symbol_variant "${symbol_lower}" "${symbol_upper}")
    string(
      REGEX REPLACE
      "(^|[^A-Za-z0-9_])${symbol_variant}([^A-Za-z0-9_]|$)"
      "\\1${symbol_variant}_64\\2"
      rewritten_content
      "${rewritten_content}"
    )
  endforeach()
endforeach()

if(EXISTS "${OUTPUT_FILE}")
  file(READ "${OUTPUT_FILE}" existing_output)
endif()

if(NOT DEFINED existing_output OR NOT existing_output STREQUAL rewritten_content)
  file(WRITE "${OUTPUT_FILE}" "${rewritten_content}")
endif()
