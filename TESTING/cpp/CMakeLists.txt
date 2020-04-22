find_package(Boost 1.48.0 QUIET COMPONENTS unit_test_framework timer)

if(NOT Boost_FOUND)
	message("Boost libraries not found. C++ tests will not be build.")
	return()
endif(NOT Boost_FOUND)



# use C++14 because of a 'auto' parameter in a lambda
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wextra -Wall -std=c++14 -pedantic")

add_definitions(-DBOOST_TEST_DYN_LINK)
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${Boost_INCLUDE_DIR})
include_directories(${CMAKE_SOURCE_DIR}/LAPACKE/include)


if(BUILD_SINGLE)
	add_executable(cpp_tests_single
		main.cpp
		sggqrcs.cpp)
	target_link_libraries(cpp_tests_single
		tmglib
		${LAPACK_LIBRARIES}
		${Boost_LIBRARIES})


	add_test(NAME cpp_tests_single COMMAND cpp_tests_single)
endif()


if(BUILD_DOUBLE)
	add_executable(cpp_tests_double
		main.cpp
		dggqrcs.cpp)
	target_link_libraries(cpp_tests_double
		tmglib
		${LAPACK_LIBRARIES}
		${Boost_LIBRARIES})


	add_test(NAME cpp_tests_double COMMAND cpp_tests_double)


	add_executable(dggsvd_vs_dggqrcs ggsvd_vs_ggqrcs.cpp)
	target_link_libraries(
		dggsvd_vs_dggqrcs ${LAPACK_LIBRARIES} ${Boost_LIBRARIES})
endif()


if(BUILD_COMPLEX)
	add_executable(cpp_tests_complex
		main.cpp
		cggqrcs.cpp)
	target_link_libraries(cpp_tests_complex
		tmglib
		${LAPACK_LIBRARIES}
		${Boost_LIBRARIES})

	add_test(NAME cpp_tests_complex COMMAND cpp_tests_complex)
endif()


if(BUILD_COMPLEX16)
	add_executable(cpp_tests_complex16
		main.cpp
		zggqrcs.cpp)
	target_link_libraries(cpp_tests_complex16
		tmglib
		${LAPACK_LIBRARIES}
		${Boost_LIBRARIES})

	add_test(NAME cpp_tests_complex16 COMMAND cpp_tests_complex16)
endif()
