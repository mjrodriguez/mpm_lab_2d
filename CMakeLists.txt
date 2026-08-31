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
    message(STATUS "Eigen3 not found on system; fetching Eigen 5.0.1")
    include(FetchContent)
    set(EIGEN_BUILD_DOC OFF CACHE BOOL "" FORCE)
    set(EIGEN_BUILD_TESTING OFF CACHE BOOL "" FORCE)
    set(BUILD_TESTING OFF CACHE BOOL "" FORCE)
    FetchContent_Declare(
        Eigen3
        GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
        GIT_TAG 5.0.1
        GIT_SHALLOW TRUE
    )
    FetchContent_MakeAvailable(Eigen3)
endif()

# nlohmann/json for reading run configs. Same find-then-fetch pattern as Eigen.
find_package(nlohmann_json 3.11 QUIET)
if(NOT nlohmann_json_FOUND)
    message(STATUS "nlohmann_json not found on system; fetching 3.12.0")
    include(FetchContent)
    FetchContent_Declare(
        nlohmann_json
        GIT_REPOSITORY https://github.com/nlohmann/json.git
        GIT_TAG v3.12.0
        GIT_SHALLOW TRUE
    )
    FetchContent_MakeAvailable(nlohmann_json)
endif()

# Simulation logic lives in a library so both the mpm executable and the
# test suite can link against it.
add_library(mpm_lib
    source/CONFIG.cpp
    source/ELASTOPLASTIC.cpp
    source/FILE_IO.cpp
    source/GRID.cpp
    source/INTERPOLATION.cpp
    source/PARTICLES.cpp
    source/SIMULATION_PARAMETERS.cpp
    source/TOOLS.cpp
)
target_include_directories(mpm_lib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(mpm_lib PUBLIC Eigen3::Eigen nlohmann_json::nlohmann_json)

add_executable(mpm main.cpp)
target_link_libraries(mpm PRIVATE mpm_lib)

enable_testing()
add_subdirectory(tests)
