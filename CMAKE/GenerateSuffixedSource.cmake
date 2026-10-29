# Use this to test GenerateSuffixedSource.cmake. Run from top of repo:
# bash CMAKE/GenerateSuffixedSource.test

# Empty list elements are load-bearing below (blank physical lines must
# survive the split/iterate/rejoin round-trip), so pin the new behavior
# explicitly: script mode does not inherit the project's policies.
cmake_policy(SET CMP0007 NEW)

# Matcher for fixed-form lines exceeding the 72-column statement field.
# CMake's regex engine has no {n} interval quantifiers, so spell the 73
# characters out.
set(_overlong_line_pat "\n")
foreach(_overlong_col RANGE 1 73)
  string(APPEND _overlong_line_pat "[^\n]")
endforeach()
string(APPEND _overlong_line_pat "[^\n]*")

# Single-file mode (INPUT_FILE/OUTPUT_FILE) and batched mode
# (INPUT_FILES/OUTPUT_FILES, semicolon-separated lists of equal length) share
# the same per-file processing below. Batching amortizes the cmake startup
# cost when CMAKE_UNITY_BUILD groups compilation into a few unity objects.
if(DEFINED INPUT_FILES OR DEFINED OUTPUT_FILES)
  if(NOT DEFINED INPUT_FILES)
    message(FATAL_ERROR "INPUT_FILES must be set in batch mode")
  endif()
  if(NOT DEFINED OUTPUT_FILES)
    message(FATAL_ERROR "OUTPUT_FILES must be set in batch mode")
  endif()
else()
  if(NOT DEFINED INPUT_FILE)
    message(FATAL_ERROR "INPUT_FILE must be set")
  endif()
  if(NOT DEFINED OUTPUT_FILE)
    message(FATAL_ERROR "OUTPUT_FILE must be set")
  endif()
  set(INPUT_FILES "${INPUT_FILE}")
  set(OUTPUT_FILES "${OUTPUT_FILE}")
  unset(INPUT_FILE)
  unset(OUTPUT_FILE)
endif()

list(LENGTH INPUT_FILES _num_input_files)
list(LENGTH OUTPUT_FILES _num_output_files)
if(NOT _num_input_files EQUAL _num_output_files)
  message(FATAL_ERROR "INPUT_FILES (${_num_input_files}) and OUTPUT_FILES (${_num_output_files}) must have the same length")
endif()

if(NOT DEFINED REPLACE_IN_STRINGS)
  set(REPLACE_IN_STRINGS ON)
endif()

if(NOT DEFINED SUFFIX)
  message(FATAL_ERROR "SUFFIX must be set")
endif()

# Optional allowlist: if given, only the listed symbols (case-insensitive)
# are suffixed. Used to route a subset of the routines to alternative
# implementations (e.g. LAPACKE test wrappers) while all other references
# keep resolving to their default implementation.
if(DEFINED SYMBOL_ALLOWLIST)
  set(_symbol_allowlist)
  foreach(_allowlist_entry IN LISTS SYMBOL_ALLOWLIST)
    string(TOLOWER "${_allowlist_entry}" _allowlist_entry)
    list(APPEND _symbol_allowlist "${_allowlist_entry}")
  endforeach()
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

# Split a newline-delimited string into physical lines for single-pass
# foreach iteration (much faster than popping lines with repeated
# full-string copies). Three hazards: file contents (comments, strings)
# may contain `;`, which would corrupt list iteration, and unmatched `[`
# or `]` (common in comments, e.g. `(VL,VU]`), which also corrupt how list
# semicolons are interpreted -- so all three are masked behind placeholders
# while scanning; and foreach skips empty elements, so blank lines are
# shielded behind a further placeholder that callers map back to "". The
# placeholders contain no regex-active characters and never appear in real
# sources. Callers that need the exact line text restore all placeholders.
macro(_split_lines text lines_var)
  string(REPLACE ";" "@@LAPACK_GEN_SEMICOLON@@" _split_masked "${${text}}")
  string(REPLACE "[" "@@LAPACK_GEN_LBRACKET@@" _split_masked "${_split_masked}")
  string(REPLACE "]" "@@LAPACK_GEN_RBRACKET@@" _split_masked "${_split_masked}")
  if(_split_masked MATCHES "^\n")
    string(REGEX REPLACE "^\n" "@@LAPACK_GEN_EMPTY_LINE@@\n" _split_masked "${_split_masked}")
  endif()
  if(_split_masked MATCHES "\n$")
    string(APPEND _split_masked "@@LAPACK_GEN_EMPTY_LINE@@")
  endif()
  while(_split_masked MATCHES "\n\n")
    string(REPLACE "\n\n" "\n@@LAPACK_GEN_EMPTY_LINE@@\n" _split_masked "${_split_masked}")
  endwhile()
  string(REPLACE "\n" ";" ${lines_var} "${_split_masked}")
endmacro()

# Map one split line back to its exact source text.
macro(_restore_line line)
  if("${${line}}" STREQUAL "@@LAPACK_GEN_EMPTY_LINE@@")
    set(${line} "")
  else()
    string(REPLACE "@@LAPACK_GEN_SEMICOLON@@" ";" ${line} "${${line}}")
    string(REPLACE "@@LAPACK_GEN_LBRACKET@@" "[" ${line} "${${line}}")
    string(REPLACE "@@LAPACK_GEN_RBRACKET@@" "]" ${line} "${${line}}")
  endif()
endmacro()

# Extract function/subroutine names from one logical Fortran statement.
function(_extract_symbols statement_text result)
  set(symbols)

  if("${statement_text}" MATCHES "(external|EXTERNAL)")
    string(REGEX REPLACE
      "^.*${CMAKE_MATCH_1}[ \t]*(::)?[ \t]*" ""
      external_names "${statement_text}")
    string(REGEX REPLACE "^[ \t]*::[ \t]*" "" external_names "${external_names}")
    string(REPLACE "," ";" external_names "${external_names}")
    foreach(external_name IN LISTS external_names)
      string(STRIP "${external_name}" external_name)
      if(NOT external_name STREQUAL "")
        list(APPEND symbols "${external_name}")
      endif()
    endforeach()
  elseif("${statement_text}" MATCHES "^ *(use|USE|module|MODULE)  *((xerbla|XERBLA)[a-zA-Z0-9_]*)")
    set(symbol_name "${CMAKE_MATCH_2}")
    list(APPEND symbols "${symbol_name}")
  elseif("${statement_text}" MATCHES "(subroutine|SUBROUTINE|function|FUNCTION)")
    string(REGEX REPLACE
      "^[a-zA-Z0-9_ *]*${CMAKE_MATCH_1}[ ]*" ""
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

  string(REPLACE "\r\n" "\n" _wrap_text "${source_text}")
  string(REPLACE "\r" "\n" _wrap_text "${_wrap_text}")
  # Fast path: fixed-form files that fit in 72 columns need no wrapping.
  if(NOT _wrap_text MATCHES "${_overlong_line_pat}")
    set(${result} "${_wrap_text}" PARENT_SCOPE)
    return()
  endif()

  # Wrap only the over-long lines. Collecting the distinct ones in a single
  # scan is far cheaper than walking every line: `;` would corrupt list
  # iteration, so match on masked text (no `;` left, elements stay intact),
  # then map each candidate back to exact text. Replacements are anchored
  # on the leading newline, so one line can never match inside another;
  # each distinct line is therefore wrapped exactly once.
  string(REPLACE ";" "@@LAPACK_GEN_SEMICOLON@@" _wrap_masked "\n${_wrap_text}")
  string(REPLACE "[" "@@LAPACK_GEN_LBRACKET@@" _wrap_masked "${_wrap_masked}")
  string(REPLACE "]" "@@LAPACK_GEN_RBRACKET@@" _wrap_masked "${_wrap_masked}")
  string(REGEX MATCHALL "${_overlong_line_pat}" _wrap_candidates "${_wrap_masked}")
  list(REMOVE_DUPLICATES _wrap_candidates)
  set(_wrap_work "\n${_wrap_text}")
  foreach(_wrap_candidate IN LISTS _wrap_candidates)
    string(REGEX REPLACE "^\n" "" _wrap_body "${_wrap_candidate}")
    string(REPLACE "@@LAPACK_GEN_SEMICOLON@@" ";" _wrap_body "${_wrap_body}")
    string(REPLACE "@@LAPACK_GEN_LBRACKET@@" "[" _wrap_body "${_wrap_body}")
    string(REPLACE "@@LAPACK_GEN_RBRACKET@@" "]" _wrap_body "${_wrap_body}")
    string(LENGTH "${_wrap_body}" _wrap_body_len)
    if(NOT _wrap_body_len GREATER 72)
      continue()
    endif()
    _wrap_fixed_form_line("${_wrap_body}" _wrap_body_wrapped)
    if(_wrap_body_wrapped STREQUAL _wrap_body)
      continue()
    endif()
    string(REPLACE "\n${_wrap_body}" "\n${_wrap_body_wrapped}" _wrap_work "${_wrap_work}")
  endforeach()
  string(SUBSTRING "${_wrap_work}" 1 -1 _wrap_text)

  set(${result} "${_wrap_text}" PARENT_SCOPE)
endfunction()

# Protect string literals in the source content by replacing them with
# placeholders, and save the original literals in variables for later restoration.
function(_protect_fortran_string_literals input_text result count_result)
  set(output_text "")
  set(literal_count 0)
  set(string_literal_regex "'([^']|'')*'|\"([^\"]|\"\")*\"")

  _split_lines(input_text _protect_lines)
  set(_protect_first TRUE)
  foreach(_protect_masked IN LISTS _protect_lines)
    if(_protect_first)
      set(_protect_first FALSE)
    else()
      string(APPEND output_text "\n")
    endif()
    # Match and replace on the masked line: literals may contain `;` or
    # unmatched brackets, which would corrupt list iteration. The saved
    # copy is restored to exact text, since the restore phase works on
    # the real (unmasked) rewritten content.
    set(current_line "${_protect_masked}")
    string(REGEX MATCHALL "${string_literal_regex}" string_literals
      "${current_line}")
    foreach(string_literal IN LISTS string_literals)
      set(placeholder "@@LAPACK_SUFFIX_STRING_LITERAL_${literal_count}@@")
      string(REPLACE "${string_literal}" "${placeholder}" current_line
        "${current_line}")
      set(_protected_literal "${string_literal}")
      _restore_line(_protected_literal)
      set(protected_string_literal_${literal_count} "${_protected_literal}"
        PARENT_SCOPE)
      math(EXPR literal_count "${literal_count} + 1")
    endforeach()
    _restore_line(current_line)

    string(APPEND output_text "${current_line}")
  endforeach()

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
      set(placeholder "@@LAPACK_SUFFIX_STRING_LITERAL_${index}@@")
      string(REPLACE "${placeholder}" "${protected_string_literal_${index}}"
        output_text "${output_text}")
    endforeach()
  endif()

  set(${result} "${output_text}" PARENT_SCOPE)
endfunction()

math(EXPR _last_pair_index "${_num_input_files} - 1")
foreach(_pair_index RANGE 0 ${_last_pair_index})
  list(GET INPUT_FILES ${_pair_index} INPUT_FILE)
  list(GET OUTPUT_FILES ${_pair_index} OUTPUT_FILE)

  # Check whether this input file is fixed or free form based on extension.
  get_filename_component(input_extension "${INPUT_FILE}" LAST_EXT)
  string(TOLOWER "${input_extension}" input_extension_lower)
  if(input_extension_lower MATCHES "^\\.f(90|95|03|08)$")
    set(_source_is_free_form TRUE)
  else()
    set(_source_is_free_form FALSE)
  endif()

  unset(existing_output)
  unset(protected_string_literal_count)

  get_filename_component(output_dir "${OUTPUT_FILE}" DIRECTORY)
  file(MAKE_DIRECTORY "${output_dir}")

  file(READ "${INPUT_FILE}" source_content)
set(rewritten_content "${source_content}")

# Extract symbol names from the source file
string(REPLACE "\r\n" "\n" symbol_scan_source "${source_content}")
string(REPLACE "\r" "\n" symbol_scan_source "${symbol_scan_source}")
set(symbol_names)

_split_lines(symbol_scan_source _scan_lines)
set(_have_logical FALSE)
set(logical_line "")
set(last_line "")
foreach(_scan_masked IN LISTS _scan_lines)
  # Map blank-line shields back; `;` stays masked while scanning (the
  # declaration patterns never need the literal `;`).
  if("${_scan_masked}" STREQUAL "@@LAPACK_GEN_EMPTY_LINE@@")
    set(current_line "")
  else()
    set(current_line "${_scan_masked}")
  endif()

  _is_fortran_comment_line("${current_line}" is_comment_line)
  if(is_comment_line)
    # A comment ends the open logical line, exactly like the lookahead
    # version breaking its continuation join on a comment line.
    if(_have_logical)
      if("${logical_line}" MATCHES "(external|EXTERNAL|subroutine|SUBROUTINE|function|FUNCTION|use|USE|module|MODULE)")
        _extract_symbols("${logical_line}" symbols)
      else()
        set(symbols)
      endif()
      if(symbols)
        list(APPEND symbol_names ${symbols})
      endif()
      set(_have_logical FALSE)
    endif()
    continue()
  endif()

  # Does this line continue the open logical line? Same test as the old
  # peek-at-next-line version, asked from the other side: the current
  # line joins when it looks like a continuation of the previous one.
  set(consume_continuation FALSE)
  if(_have_logical)
    if(_source_is_free_form)
      if("${last_line}" MATCHES "&[ \t]*$" OR
          "${current_line}" MATCHES "^[ \t]*&")
        set(consume_continuation TRUE)
      endif()
    else()
      _is_fixed_form_continuation_line("${current_line}"
        is_fixed_form_continuation)
      if(is_fixed_form_continuation)
        set(consume_continuation TRUE)
      endif()
    endif()
  endif()

  _normalize_fortran_statement_fragment("${current_line}" fragment)
  if(consume_continuation)
    string(APPEND logical_line " ${fragment}")
  else()
    if(_have_logical)
      # Cheap trigger: most lines cannot declare anything, so skip the
      # heavier _extract_symbols call unless a declaration keyword is
      # present (it applies the precise anchored patterns itself).
      if("${logical_line}" MATCHES "(external|EXTERNAL|subroutine|SUBROUTINE|function|FUNCTION|use|USE|module|MODULE)")
        _extract_symbols("${logical_line}" symbols)
      else()
        set(symbols)
      endif()
      if(symbols)
        list(APPEND symbol_names ${symbols})
      endif()
    endif()
    set(logical_line "${fragment}")
    set(_have_logical TRUE)
  endif()
  set(last_line "${current_line}")
endforeach()
if(_have_logical)
  if("${logical_line}" MATCHES "(external|EXTERNAL|subroutine|SUBROUTINE|function|FUNCTION|use|USE|module|MODULE)")
    _extract_symbols("${logical_line}" symbols)
  else()
    set(symbols)
  endif()
  if(symbols)
    list(APPEND symbol_names ${symbols})
  endif()
endif()
list(REMOVE_DUPLICATES symbol_names)
list(REMOVE_ITEM symbol_names ETIME etime ETIME_ etime_)

# Restrict the renaming to allowlisted symbols if an allowlist was given
if(DEFINED SYMBOL_ALLOWLIST)
  set(_filtered_symbol_names)
  foreach(symbol_name IN LISTS symbol_names)
    string(TOLOWER "${symbol_name}" _symbol_name_lower)
    list(FIND _symbol_allowlist "${_symbol_name_lower}" _symbol_allowlist_index)
    if(NOT _symbol_allowlist_index EQUAL -1)
      list(APPEND _filtered_symbol_names "${symbol_name}")
    endif()
  endforeach()
  set(symbol_names ${_filtered_symbol_names})
endif()

# If string literals should not be modified, protect them before performing replacements
if(NOT REPLACE_IN_STRINGS)
  _protect_fortran_string_literals(
    "${rewritten_content}" rewritten_content protected_string_literal_count)
endif()

# Replace symbol names with their suffixed versions in the source content.
# One alternation over the whole file replaces one full-file regex scan per
# symbol and case variant. Matches are non-overlapping within a scan, so a
# name adjacent to a renamed one (e.g. `X = FOO()+BAR()`) can be missed when
# the shared boundary character is consumed; repeat to quiescence (renaming
# only appends word characters, so blocked names become matchable and no
# pass can undo an earlier one).
if(symbol_names)
  set(_rename_alternatives)
  foreach(symbol_name IN LISTS symbol_names)
    string(TOLOWER "${symbol_name}" symbol_lower)
    string(TOUPPER "${symbol_name}" symbol_upper)
    foreach(symbol_variant "${symbol_name}" "${symbol_lower}" "${symbol_upper}")
      list(APPEND _rename_alternatives "${symbol_variant}")
    endforeach()
  endforeach()
  list(REMOVE_DUPLICATES _rename_alternatives)
  # Longest first so overlapping names prefer the longest match.
  set(_sorted_alternatives)
  foreach(_alternative IN LISTS _rename_alternatives)
    string(LENGTH "${_alternative}" _alternative_len)
    set(_inserted FALSE)
    set(_new_sorted)
    foreach(_current IN LISTS _sorted_alternatives)
      if(NOT _inserted)
        string(LENGTH "${_current}" _current_len)
        if(_alternative_len GREATER _current_len)
          list(APPEND _new_sorted "${_alternative}")
          set(_inserted TRUE)
        endif()
      endif()
      list(APPEND _new_sorted "${_current}")
    endforeach()
    if(NOT _inserted)
      list(APPEND _new_sorted "${_alternative}")
    endif()
    set(_sorted_alternatives "${_new_sorted}")
  endforeach()
  list(JOIN _sorted_alternatives "|" _rename_pattern)
  set(_rename_regex "(^|[^A-Za-z0-9_])(${_rename_pattern})([^A-Za-z0-9_]|$)")
  set(_rename_replacement "\\1\\2${SUFFIX}\\3")
  set(_rename_pass 0)
  while(_rename_pass LESS 10)
    set(_rename_before "${rewritten_content}")
    string(REGEX REPLACE "${_rename_regex}" "${_rename_replacement}"
      rewritten_content "${rewritten_content}")
    if(rewritten_content STREQUAL _rename_before)
      break()
    endif()
    math(EXPR _rename_pass "${_rename_pass} + 1")
  endwhile()
endif()

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
endforeach()
