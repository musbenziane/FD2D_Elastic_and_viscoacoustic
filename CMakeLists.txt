cmake_minimum_required(VERSION 3.19)
project(FD2D_ELASTIC Fortran)


enable_language(Fortran)
add_executable(FD2D_ELASTIC main.f90 readmodelfiles.f90 ricker.f90)
