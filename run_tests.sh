#!/bin/env bash

# this script is responsible for running the selected precision scripts
POSITIONAL_ARGS=()
RUN_SINGLE=""
RUN_DOUBLE=""
RUN_COMPLEX=""
RUN_COMPLEX_16=""

while [[ $# -gt 0 ]]; do
  case $1 in
    -s|--single)
      RUN_SINGLE="1"
      shift
      ;;
    -d|--double)
      RUN_DOUBLE="1"
      shift
      ;;
    -c|--complex)
      RUN_COMPLEX="1"
      shift
      ;;
    -z|--doublecomplex)
      RUN_COMPLEX_16="1"
      shift
      ;;
    *)
      # ignore bad input
      shift
      ;;
  esac
done

TEST_LINE=""
if [ -n "$RUN_SINGLE" ]; then
  TEST_LINE="${TEST_LINE} single"
fi
if [ -n "$RUN_DOUBLE" ]; then
  TEST_LINE="${TEST_LINE} double"
fi
if [ -n "$RUN_COMPLEX" ]; then
  TEST_LINE="${TEST_LINE} complex"
fi
if [ -n "$RUN_COMPLEX_16" ]; then
  TEST_LINE="${TEST_LINE} complex16"
fi

# Now we recompile all necessary libraries
make clean

make -C SRC/ -j

make -C BLAS/ -j

make tmglib -j4

# Run our desired tests
make -C TESTING $TEST_LINE -j4

# Print the testing info to the terminal
./lapack_testing.py
