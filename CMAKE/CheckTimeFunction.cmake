# - Check if the Fortran function exists.
# CHECK_TIME_FUNCTION(FUNCTION VARIABLE TYPE)
# - macro which checks if the Fortran function exists
#  FUNCTION - the name of the Fortran function
#  VARIABLE - variable to store the result
#

macro(CHECK_TIME_FUNCTION FUNCTION VARIABLE)

  try_compile(RES
    ${PROJECT_BINARY_DIR}/INSTALL
    ${PROJECT_SOURCE_DIR}/INSTALL
    TIMING secondtst_${FUNCTION}
    CMAKE_FLAGS
      -DCMAKE_OSX_DEPLOYMENT_TARGET:STRING=${CMAKE_OSX_DEPLOYMENT_TARGET}
      -DCMAKE_Fortran_FLAGS:STRING=${CMAKE_Fortran_FLAGS}
      -DCMAKE_EXE_LINKER_FLAGS:STRING=${CMAKE_EXE_LINKER_FLAGS}
      -DCMAKE_VERBOSE_MAKEFILE=ON
    OUTPUT_VARIABLE OUTPUT)

  if(RES)
    set(${VARIABLE} ${FUNCTION} CACHE INTERNAL "Have Fortran function ${FUNCTION}")
    message(STATUS "Looking for Fortran ${FUNCTION} - found")
    file(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeOutput.log
      "Fortran ${FUNCTION} exists. ${OUTPUT} \n\n")
  else()
    message(STATUS "Looking for Fortran ${FUNCTION} - not found")
    file(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeError.log
      "Fortran ${FUNCTION} does not exist. \n ${OUTPUT} \n")
  endif()
endmacro()
