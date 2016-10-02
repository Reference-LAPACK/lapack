find_package(Boost 1.48.0 QUIET COMPONENTS unit_test_framework timer)

if(NOT Boost_FOUND)
	message("Boost libraries not found. C++ tests will not be build.")
	return()
endif(NOT Boost_FOUND)



set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wextra -Wall -std=c++14 -pedantic")

add_definitions(-DBOOST_TEST_DYN_LINK)
include_directories(${CMAKE_CURRENT_SOURCE_DIR} ${Boost_INCLUDE_DIR})


if(BUILD_DOUBLE)
	add_executable(cpp_tests_double
		main.cpp
		dggqrcs.cpp)
	target_link_libraries(cpp_tests_double
		tmglib
		${LAPACK_LIBRARIES}
		${Boost_LIBRARIES})


	add_test(NAME cpp_tests_double COMMAND cpp_tests_double)


	add_executable(ggsvd_vs_ggqrcs ggsvd_vs_ggqrcs.cpp)
	target_link_libraries(
		ggsvd_vs_ggqrcs ${LAPACK_LIBRARIES} ${Boost_LIBRARIES})
endif()
