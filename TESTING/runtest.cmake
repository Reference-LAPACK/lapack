set(ARGS )
if(DEFINED OUTPUT)
  set(ARGS OUTPUT_FILE "${OUTPUT}"  ERROR_FILE "${OUTPUT}.err")
endif()
if(DEFINED INPUT)
  list(APPEND ARGS INPUT_FILE "${INPUT}")
endif()
message("Running: ${TEST}")
message("ARGS= ${ARGS}")
execute_process(COMMAND "${TEST}"
  ${ARGS}
  RESULT_VARIABLE RET)
if(DEFINED OUTPUT)
  file(READ "${OUTPUT}" TEST_OUTPUT)
  file(READ "${OUTPUT}.err" TEST_ERROR)
  message("Test OUTPUT:\n${TEST_OUTPUT}")
  message("Test ERROR:\n${TEST_ERROR}")
endif()

# if the test does not return 0, then fail it
if(NOT RET STREQUAL "0")
  message(FATAL_ERROR "Test ${TEST} returned ${RET}")
endif()
message( "Test ${TEST} returned ${RET}")
