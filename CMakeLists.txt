cmake_minimum_required(VERSION 3.14)
project(mpm_lab_2d CXX)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pedantic")
set(CMAKE_CXX_FLAGS_RELEASE "-O2")

# Simulation logic lives in a library so both the mpm executable and the
# test suite can link against it.
add_library(mpm_lib
    source/ELASTOPLASTIC.cpp
    source/FILE_IO.cpp
    source/GRID.cpp
    source/INTERPOLATION.cpp
    source/PARTICLES.cpp
    source/SIMULATION_PARAMETERS.cpp
    source/TOOLS.cpp
)
target_include_directories(mpm_lib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(mpm main.cpp)
target_link_libraries(mpm PRIVATE mpm_lib)

enable_testing()
add_subdirectory(tests)
