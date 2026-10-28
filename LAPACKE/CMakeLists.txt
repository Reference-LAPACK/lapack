add_library(${LAPACKELIB}_obj OBJECT)
target_sources(${LAPACKELIB} PRIVATE $<TARGET_OBJECTS:${LAPACKELIB}_obj>)
if(BUILD_INDEX64_EXT_API)
  add_library(${LAPACKELIB}_64_obj OBJECT)
  target_compile_definitions(${LAPACKELIB}_64_obj PRIVATE
    LAPACK_ILP64 LAPACKE_API64 WeirdNEC CBLAS_API64)
  target_sources(${LAPACKELIB} PRIVATE $<TARGET_OBJECTS:${LAPACKELIB}_64_obj>)
endif()
add_subdirectory(src)
add_subdirectory(utils)
lapack_configure_target(${LAPACKELIB}_obj)
if(BUILD_INDEX64_EXT_API)
  lapack_configure_target(${LAPACKELIB}_64_obj)
endif()
if(BUILD_TESTING)
  if(LAPACKE_BUILD_DOUBLE)
    add_subdirectory(example)
  endif()
  add_subdirectory(testing)
endif()

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/lapacke.pc.in ${CMAKE_CURRENT_BINARY_DIR}/${LAPACKELIB}.pc @ONLY)
install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/${LAPACKELIB}.pc
  DESTINATION ${PKG_CONFIG_DIR}
  COMPONENT Development
  )
