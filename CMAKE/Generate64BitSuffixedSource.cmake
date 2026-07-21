if(NOT DEFINED INPUT_FILE)
  message(FATAL_ERROR "INPUT_FILE must be set")
endif()

if(NOT DEFINED OUTPUT_FILE)
  message(FATAL_ERROR "OUTPUT_FILE must be set")
endif()

if(NOT DEFINED REPLACE_IN_STRINGS)
  set(REPLACE_IN_STRINGS ON)
endif()

# Check whether the input file is fixed or free form based on its extension.
get_filename_component(input_extension "${INPUT_FILE}" LAST_EXT)
string(TOLOWER "${input_extension}" input_extension_lower)
if(input_extension_lower MATCHES "^\\.f(90|95|03|08)$")
  set(_source_is_free_form TRUE)
else()
  set(_source_is_free_form FALSE)
endif()

# Check if line is a comment line based on the form of the Fortran source.
function(_is_fortran_comment_line line result)
  if(_source_is_free_form)
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
  if(_source_is_free_form)
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
  if(_source_is_free_form)
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

# Pop one physical source line from a newline-delimited string. Do not use a
# CMake list here: source comments may contain unmatched brackets, which changes
# how list semicolons are interpreted.
macro(_pop_source_line remaining_source line result)
  if("${${remaining_source}}" STREQUAL "")
    set(${line} "")
    set(${result} FALSE)
  else()
    string(FIND "${${remaining_source}}" "\n" _source_line_newline_index)
    if(_source_line_newline_index EQUAL -1)
      set(${line} "${${remaining_source}}")
      set(${remaining_source} "")
    else()
      string(SUBSTRING "${${remaining_source}}" 0
        ${_source_line_newline_index} ${line})
      math(EXPR _source_line_next_index
        "${_source_line_newline_index} + 1")
      string(SUBSTRING "${${remaining_source}}"
        ${_source_line_next_index} -1 ${remaining_source})
    endif()
    set(${result} TRUE)
  endif()
endmacro()

# Peek at the next physical line without consuming it.
macro(_peek_source_line remaining_source line result)
  if("${${remaining_source}}" STREQUAL "")
    set(${line} "")
    set(${result} FALSE)
  else()
    string(FIND "${${remaining_source}}" "\n" _source_line_newline_index)
    if(_source_line_newline_index EQUAL -1)
      set(${line} "${${remaining_source}}")
    else()
      string(SUBSTRING "${${remaining_source}}" 0
        ${_source_line_newline_index} ${line})
    endif()
    set(${result} TRUE)
  endif()
endmacro()

# Extract function/subroutine names from one logical Fortran statement.
function(_extract_symbols statement_text result)
  set(symbols)

  if("${statement_text}" MATCHES "(external|EXTERNAL)")
    string(REGEX REPLACE
      "^.*(external|EXTERNAL)[ \t]*(::)?[ \t]*" ""
      external_names "${statement_text}")
    string(REGEX REPLACE "^[ \t]*::[ \t]*" "" external_names "${external_names}")
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

# Find a safe fixed-form split point at or before the 72-column limit.
function(_find_fixed_form_split line result)
  set(split_pos -1)
  set(space_pos -1)
  set(index 71)

  while(index GREATER 6)
    string(SUBSTRING "${line}" ${index} 1 current_char)
    if(current_char STREQUAL ",")
      math(EXPR split_pos "${index} + 1")
      break()
    elseif(space_pos LESS 0 AND current_char MATCHES "[ \t]")
      set(space_pos ${index})
    endif()
    math(EXPR index "${index} - 1")
  endwhile()

  if(split_pos LESS 0 AND NOT space_pos LESS 0)
    set(split_pos ${space_pos})
  endif()

  set(${result} ${split_pos} PARENT_SCOPE)
endfunction()

# Wrap one fixed-form physical line so generated sources remain valid with the
# standard 72-column statement field.
function(_wrap_fixed_form_line line result)
  _is_fortran_comment_line("${line}" is_comment_line)
  if(is_comment_line OR "${line}" MATCHES "^#")
    set(${result} "${line}" PARENT_SCOPE)
    return()
  endif()

  set(wrapped_line "")
  set(current_line "${line}")

  while(1)
    string(LENGTH "${current_line}" line_length)
    if(NOT line_length GREATER 72)
      if(wrapped_line STREQUAL "")
        set(wrapped_line "${current_line}")
      else()
        string(APPEND wrapped_line "\n${current_line}")
      endif()
      break()
    endif()

    _find_fixed_form_split("${current_line}" split_pos)
    if(NOT split_pos GREATER 6)
      if(wrapped_line STREQUAL "")
        set(wrapped_line "${current_line}")
      else()
        string(APPEND wrapped_line "\n${current_line}")
      endif()
      break()
    endif()

    string(SUBSTRING "${current_line}" 0 ${split_pos} line_head)
    string(SUBSTRING "${current_line}" ${split_pos} -1 line_tail)
    string(STRIP "${line_tail}" line_tail)

    if(wrapped_line STREQUAL "")
      set(wrapped_line "${line_head}")
    else()
      string(APPEND wrapped_line "\n${line_head}")
    endif()

    set(current_line "     $ ${line_tail}")
  endwhile()

  set(${result} "${wrapped_line}" PARENT_SCOPE)
endfunction()

# Wrap fixed-form Fortran source code to ensure it remains valid with
# the standard 72-column statement field.
function(_wrap_fixed_form_source source_text result)
  if(_source_is_free_form)
    set(${result} "${source_text}" PARENT_SCOPE)
    return()
  endif()

  string(REPLACE "\r\n" "\n" remaining_source "${source_text}")
  string(REPLACE "\r" "\n" remaining_source "${remaining_source}")
  set(wrapped_source "")

  while(1)
    string(FIND "${remaining_source}" "\n" newline_index)
    if(newline_index EQUAL -1)
      if(NOT remaining_source STREQUAL "")
        _wrap_fixed_form_line("${remaining_source}" wrapped_line)
        string(APPEND wrapped_source "${wrapped_line}")
      endif()
      break()
    endif()

    string(SUBSTRING "${remaining_source}" 0 ${newline_index} source_line)
    math(EXPR next_index "${newline_index} + 1")
    string(SUBSTRING "${remaining_source}" ${next_index} -1 remaining_source)

    _wrap_fixed_form_line("${source_line}" wrapped_line)
    string(APPEND wrapped_source "${wrapped_line}\n")
  endwhile()

  set(${result} "${wrapped_source}" PARENT_SCOPE)
endfunction()

# Protect string literals in the source content by replacing them with
# placeholders, and save the original literals in variables for later restoration.
function(_protect_fortran_string_literals input_text result count_result)
  set(output_text "")
  set(remaining_source "${input_text}")
  set(literal_count 0)
  set(string_literal_regex "'([^']|'')*'|\"([^\"]|\"\")*\"")

  while(NOT remaining_source STREQUAL "")
    string(FIND "${remaining_source}" "\n" newline_index)
    if(newline_index EQUAL -1)
      set(current_line "${remaining_source}")
      set(remaining_source "")
      set(has_newline FALSE)
    else()
      string(SUBSTRING "${remaining_source}" 0 ${newline_index} current_line)
      math(EXPR next_index "${newline_index} + 1")
      string(SUBSTRING "${remaining_source}" ${next_index} -1
        remaining_source)
      set(has_newline TRUE)
    endif()

    string(REGEX MATCHALL "${string_literal_regex}" string_literals
      "${current_line}")
    foreach(string_literal IN LISTS string_literals)
      set(placeholder "@@LAPACK_64_STRING_LITERAL_${literal_count}@@")
      string(REPLACE "${string_literal}" "${placeholder}" current_line
        "${current_line}")
      set(protected_string_literal_${literal_count} "${string_literal}"
        PARENT_SCOPE)
      math(EXPR literal_count "${literal_count} + 1")
    endforeach()

    string(APPEND output_text "${current_line}")
    if(has_newline)
      string(APPEND output_text "\n")
    endif()
  endwhile()

  set(${result} "${output_text}" PARENT_SCOPE)
  set(${count_result} "${literal_count}" PARENT_SCOPE)
endfunction()

# Restore string literals in the rewritten source content by replacing the
# placeholders with the original literals saved from the input content.
function(_restore_fortran_string_literals input_text literal_count result)
  set(output_text "${input_text}")
  if(literal_count GREATER 0)
    math(EXPR last_literal_index "${literal_count} - 1")
    foreach(index RANGE 0 ${last_literal_index})
      set(placeholder "@@LAPACK_64_STRING_LITERAL_${index}@@")
      string(REPLACE "${placeholder}" "${protected_string_literal_${index}}"
        output_text "${output_text}")
    endforeach()
  endif()

  set(${result} "${output_text}" PARENT_SCOPE)
endfunction()

get_filename_component(output_dir "${OUTPUT_FILE}" DIRECTORY)
file(MAKE_DIRECTORY "${output_dir}")

file(READ "${INPUT_FILE}" source_content)
set(rewritten_content "${source_content}")

# Extract symbol names from the source file
string(REPLACE "\r\n" "\n" symbol_scan_source "${source_content}")
string(REPLACE "\r" "\n" symbol_scan_source "${symbol_scan_source}")
set(symbol_names)

while(1)
  _pop_source_line(symbol_scan_source current_line has_line)
  if(NOT has_line)
    break()
  endif()

  _is_fortran_comment_line("${current_line}" is_comment_line)
  if(is_comment_line)
    continue()
  endif()

  _normalize_fortran_statement_fragment("${current_line}" logical_line)
  set(last_line "${current_line}")

  while(1)
    _peek_source_line(symbol_scan_source continuation_line has_line)
    if(NOT has_line)
      break()
    endif()

    _is_fortran_comment_line("${continuation_line}" continuation_is_comment)
    if(continuation_is_comment)
      break()
    endif()

    set(consume_continuation FALSE)
    if(_source_is_free_form)
      if("${last_line}" MATCHES "&[ \t]*$" OR
          "${continuation_line}" MATCHES "^[ \t]*&")
        set(consume_continuation TRUE)
      endif()
    else()
      _is_fixed_form_continuation_line("${continuation_line}"
        is_fixed_form_continuation)
      if(is_fixed_form_continuation)
        set(consume_continuation TRUE)
      endif()
    endif()

    if(NOT consume_continuation)
      break()
    endif()

    _pop_source_line(symbol_scan_source continuation_line has_line)
    _normalize_fortran_statement_fragment("${continuation_line}"
      continuation_fragment)
    string(APPEND logical_line " ${continuation_fragment}")
    set(last_line "${continuation_line}")
  endwhile()

  _extract_symbols("${logical_line}" symbols)
  if(symbols)
    list(APPEND symbol_names ${symbols})
  endif()
endwhile()
list(REMOVE_DUPLICATES symbol_names)
list(REMOVE_ITEM symbol_names ETIME etime ETIME_ etime_)

# If string literals should not be modified, protect them before performing replacements
if(NOT REPLACE_IN_STRINGS)
  _protect_fortran_string_literals(
    "${rewritten_content}" rewritten_content protected_string_literal_count)
endif()

# Replace symbol names with their suffixed versions in the source content
foreach(symbol_name IN LISTS symbol_names)
  string(TOLOWER "${symbol_name}" symbol_lower)
  string(TOUPPER "${symbol_name}" symbol_upper)
  foreach(symbol_variant "${symbol_name}" "${symbol_lower}" "${symbol_upper}")
    set(match_regex "(^|[^A-Za-z0-9_])${symbol_variant}([^A-Za-z0-9_]|$)")
    set(replacement "\\1${symbol_variant}_64\\2")
    string(REGEX REPLACE
      "${match_regex}" "${replacement}"
      rewritten_content "${rewritten_content}")
  endforeach()
endforeach()

# Restore string literals if they were protected
if(NOT REPLACE_IN_STRINGS)
  _restore_fortran_string_literals(
    "${rewritten_content}" "${protected_string_literal_count}" rewritten_content)
endif()

_wrap_fixed_form_source("${rewritten_content}" rewritten_content)

if(EXISTS "${OUTPUT_FILE}")
  file(READ "${OUTPUT_FILE}" existing_output)
endif()

if(NOT DEFINED existing_output OR NOT existing_output STREQUAL rewritten_content)
  file(WRITE "${OUTPUT_FILE}" "${rewritten_content}")
endif()
