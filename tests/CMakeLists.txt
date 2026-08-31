include(FetchContent)
FetchContent_Declare(
    Catch2
    GIT_REPOSITORY https://github.com/catchorg/Catch2.git
    GIT_TAG v3.7.1
    GIT_SHALLOW TRUE
)
FetchContent_MakeAvailable(Catch2)

add_executable(unit_tests
    sample_test.cpp
    tools_test.cpp
    interpolation_test.cpp
    grid_test.cpp
    elastoplastic_test.cpp
    particles_test.cpp
    levelset_test.cpp
    single_timestep_test.cpp
    file_io_test.cpp
    config_test.cpp
)
target_link_libraries(unit_tests PRIVATE mpm_lib Catch2::Catch2WithMain)

list(APPEND CMAKE_MODULE_PATH ${catch2_SOURCE_DIR}/extras)
include(Catch)
catch_discover_tests(unit_tests)
