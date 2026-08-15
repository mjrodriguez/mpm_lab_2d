cmake_minimum_required(VERSION 3.10)
project(mpm_lab_2d CXX)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pedantic")
set(CMAKE_CXX_FLAGS_RELEASE "-O2")

add_executable(mpm
    main.cpp
    source/ELASTOPLASTIC.cpp
    source/FILE_IO.cpp
    source/GRID.cpp
    source/INTERPOLATION.cpp
    source/PARTICLES.cpp
    source/SIMULATION_PARAMETERS.cpp
    source/TOOLS.cpp
)
