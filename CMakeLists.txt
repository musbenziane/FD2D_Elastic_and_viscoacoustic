cmake_minimum_required(VERSION 3.19)
project(FD2D_ELASTIC Fortran)
if(CMAKE_Fortran_COMPILER_ID MATCHES GNU)
    set(CMAKE_Fortran_FLAGS         "-fopenmp")
    set(CMAKE_Fortran_FLAGS         "${CMAKE_Fortran_FLAGS} -Wall")
    set(CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g3")
    set(CMAKE_Fortran_FLAGS_RELEASE "-Ofast -march=native")
endif()


enable_language(Fortran)
add_executable(FD2D_ELASTIC main.f90 readmodelfiles.f90 ricker.f90 read_acqui.f90)
