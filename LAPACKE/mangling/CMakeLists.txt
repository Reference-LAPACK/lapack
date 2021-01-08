cmake_minimum_required(VERSION 3.2)
project(MANGLING C Fortran)

add_executable(xintface Fintface.f Cintface.c)
