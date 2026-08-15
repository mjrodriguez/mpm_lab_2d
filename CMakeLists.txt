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

# Eigen is header-only. Prefer a system install (e.g. `brew install eigen` /
# `apt install libeigen3-dev`); if it isn't found, fetch a pinned version
# instead of vendoring the source tree in this repo.
find_package(Eigen3 3.3 QUIET NO_MODULE)
if(NOT Eigen3_FOUND)
    message(STATUS "Eigen3 not found on system; fetching Eigen 3.4.0")
    include(FetchContent)
    set(EIGEN_BUILD_DOC OFF CACHE BOOL "" FORCE)
    set(EIGEN_BUILD_TESTING OFF CACHE BOOL "" FORCE)
    set(BUILD_TESTING OFF CACHE BOOL "" FORCE)
    FetchContent_Declare(
        Eigen3
        GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
        GIT_TAG 3.4.0
        GIT_SHALLOW TRUE
    )
    FetchContent_MakeAvailable(Eigen3)
endif()

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

target_link_libraries(mpm PRIVATE Eigen3::Eigen)
