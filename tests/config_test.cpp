#include <cstdio>
#include <fstream>
#include <stdexcept>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "include/CONFIG.h"

namespace {

// Writes `content` to a fixed temp path and returns it. Tests remove the
// file themselves once done with it.
string WriteTempConfig(const string& filename, const string& content) {
    string path = "/tmp/" + filename;
    ofstream file(path);
    file << content;
    file.close();
    return path;
}

}  // namespace

TEST_CASE("LoadConfig reads a valid config file", "[config]") {
    string path = WriteTempConfig("mpm_test_valid.json", R"({
        "output_directory": "./output",
        "simulation_name": "test_run",
        "grid": "GRID_HAT_OBSTACLE",
        "particles": "RECTANGLE_FREEFALL_HAT_OBSTACLE",
        "parameters": "DEFAULT_PARAMETERS"
    })");

    RUN_CONFIG config = LoadConfig(path);

    REQUIRE(config.outputDirectory == "./output");
    REQUIRE(config.simulationName == "test_run");
    REQUIRE(config.gridType == "GRID_HAT_OBSTACLE");
    REQUIRE(config.particlesType == "RECTANGLE_FREEFALL_HAT_OBSTACLE");
    REQUIRE(config.parametersType == "DEFAULT_PARAMETERS");

    remove(path.c_str());
}

TEST_CASE("LoadConfig throws when the file doesn't exist", "[config]") {
    REQUIRE_THROWS_WITH(
        LoadConfig("/tmp/mpm_test_definitely_does_not_exist.json"),
        Catch::Matchers::ContainsSubstring("Could not open config file"));
}

TEST_CASE("LoadConfig throws on invalid JSON", "[config]") {
    string path = WriteTempConfig("mpm_test_invalid.json", "{ this is not valid json");

    REQUIRE_THROWS_AS(LoadConfig(path), runtime_error);

    remove(path.c_str());
}

TEST_CASE("LoadConfig throws when a required key is missing", "[config]") {
    string path = WriteTempConfig("mpm_test_missing_key.json", R"({
        "output_directory": "./output",
        "simulation_name": "test_run",
        "grid": "GRID_HAT_OBSTACLE",
        "particles": "RECTANGLE_FREEFALL_HAT_OBSTACLE"
    })");  // "parameters" is missing

    REQUIRE_THROWS_WITH(LoadConfig(path), Catch::Matchers::ContainsSubstring("parameters"));

    remove(path.c_str());
}

TEST_CASE("LoadConfig throws when a required key isn't a string", "[config]") {
    string path = WriteTempConfig("mpm_test_wrong_type.json", R"({
        "output_directory": "./output",
        "simulation_name": "test_run",
        "grid": 42,
        "particles": "RECTANGLE_FREEFALL_HAT_OBSTACLE",
        "parameters": "DEFAULT_PARAMETERS"
    })");

    REQUIRE_THROWS_WITH(LoadConfig(path), Catch::Matchers::ContainsSubstring("grid"));

    remove(path.c_str());
}
