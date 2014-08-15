cmake_minimum_required(VERSION 2.8.10)
project(CBLAS C)
enable_language(Fortran)

include_directories( include )
add_subdirectory(src)

