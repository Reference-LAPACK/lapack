cmake_minimum_required(VERSION 3.13)
project(TIMING Fortran)

# Add the CMake directory for custom CMake modules
set(CMAKE_MODULE_PATH "${TIMING_SOURCE_DIR}/../CMAKE" ${CMAKE_MODULE_PATH})

# Check for any necessary platform specific compiler flags
include(CheckLAPACKCompilerFlags)
CheckLAPACKCompilerFlags()

add_executable(secondtst_NONE second_NONE.f secondtst.f)
add_executable(secondtst_EXT_ETIME second_EXT_ETIME.f secondtst.f)
add_executable(secondtst_EXT_ETIME_ second_EXT_ETIME_.f secondtst.f)
add_executable(secondtst_INT_ETIME second_INT_ETIME.f secondtst.f)
add_executable(secondtst_INT_CPU_TIME second_INT_CPU_TIME.f secondtst.f)
add_executable(testieee tstiee.f)
add_executable(testversion ilaver.f LAPACK_version.f)

if( NOT TARGET test_zcomplexabs )
  add_executable( test_zcomplexabs test_zcomplexabs.f )
endif()
if( NOT TARGET test_zcomplexdiv )
  add_executable( test_zcomplexdiv test_zcomplexdiv.f )
endif()
if( NOT TARGET test_zcomplexmult )
  add_executable( test_zcomplexmult test_zcomplexmult.f )
endif()
if( NOT TARGET test_zminMax )
  add_executable( test_zminMax test_zminMax.f )
endif()
