cmake_minimum_required(VERSION 3.18)
project(MANGLING C Fortran)

add_executable(xintface Fintface.f Cintface.c)
