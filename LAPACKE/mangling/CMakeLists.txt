cmake_minimum_required(VERSION 3.18)
project(MANGLING C Fortran)

add_executable(xintface Fintface.f Cintface.c)
# Single test tool mixing Fortran and C; keep it out of unity batches.
set_target_properties(xintface PROPERTIES UNITY_BUILD OFF)
