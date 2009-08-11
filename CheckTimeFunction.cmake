# - Check if the Fortran function exists.
# CHECK_TIME_FUNCTION(FUNCTION VARIABLE TYPE)
# - macro which checks if the Fortran function exists
#  FUNCTION - the name of the Fortran function
#  VARIABLE - variable to store the result
#

macro(CHECK_TIME_FUNCTION FUNCTION VARIABLE)

    execute_process(COMMAND ${CMAKE_Fortran_COMPILER} -o secondtst second_${FUNCTION}.f secondtst.f
    WORKING_DIRECTORY ${LAPACK_SOURCE_DIR}/INSTALL
    RESULT_VARIABLE RES
    OUTPUT_VARIABLE OUT
    ERROR_VARIABLE  ERR
    )

    if(${RES} EQUAL 0)
      set(${VARIABLE} ${FUNCTION} CACHE INTERNAL "Have Fortran function ${FUNCTION}")
      message(STATUS "Looking for Fortran ${FUNCTION} - found")
      file(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeOutput.log 
        "Fortran ${FUNCTION} exists. ${OUT} \n\n")
    else(${RES} EQUAL 0)
      message(STATUS "Looking for Fortran ${FUNCTION} - not found")
      file(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeError.log 
        "Fortran ${FUNCTION} does not exist. \n ${ERR} \n")
    endif(${RES} EQUAL 0)
endmacro(CHECK_TIME_FUNCTION)
