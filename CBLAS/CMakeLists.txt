if(WIN32)
  set(HAS_ATTRIBUTE_WEAK_SUPPORT FALSE)
else()
  include(CheckCSourceCompiles)
  check_c_source_compiles("int __attribute__((weak)) main(void) { return 0; }"
    HAS_ATTRIBUTE_WEAK_SUPPORT)
endif()
add_subdirectory(src)
if(BUILD_TESTING)
  add_subdirectory(testing)
  if(BUILD_DOUBLE)
    add_subdirectory(examples)
  endif()
endif()

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cblas.pc.in
  ${CMAKE_CURRENT_BINARY_DIR}/${CBLASLIB}.pc @ONLY)
install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/${CBLASLIB}.pc
  DESTINATION ${PKG_CONFIG_DIR}
  COMPONENT Development
  )
