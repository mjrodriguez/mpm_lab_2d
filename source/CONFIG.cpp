#include "../include/CONFIG.h"

#include <fstream>
#include <stdexcept>

#include <nlohmann/json.hpp>

using json = nlohmann::json;

namespace {

std::string RequireString(const json& j, const std::string& key, const std::string& configPath) {
    if (!j.contains(key)) {
        throw std::runtime_error("Config \"" + configPath + "\" is missing required key \"" + key + "\"");
    }
    if (!j.at(key).is_string()) {
        throw std::runtime_error("Config \"" + configPath + "\" key \"" + key + "\" must be a string");
    }
    return j.at(key).get<std::string>();
}

}  

RUN_CONFIG LoadConfig(const std::string& path) {
    std::ifstream file(path);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open config file \"" + path + "\"");
    }

    json j;
    try {
        file >> j;
    } catch (const json::parse_error& e) {
        throw std::runtime_error("Failed to parse \"" + path + "\" as JSON: " + e.what());
    }

    RUN_CONFIG config;
    config.outputDirectory = RequireString(j, "output_directory", path);
    config.simulationName = RequireString(j, "simulation_name", path);
    config.gridType = RequireString(j, "grid", path);
    config.particlesType = RequireString(j, "particles", path);
    config.parametersType = RequireString(j, "parameters", path);
    return config;
}

std::unique_ptr<GRID> MakeGrid(const std::string& name) {
    if (name == "GRID_NO_OBSTACLE") return std::make_unique<GRID_NO_OBSTACLE>();
    if (name == "GRID_CYLINDER_OBSTACLE") return std::make_unique<GRID_CYLINDER_OBSTACLE>();
    if (name == "GRID_HAT_OBSTACLE") return std::make_unique<GRID_HAT_OBSTACLE>();

    throw std::runtime_error(
        "Unknown grid type \"" + name + "\". Valid options: "
        "GRID_NO_OBSTACLE, GRID_CYLINDER_OBSTACLE, GRID_HAT_OBSTACLE");
}

std::unique_ptr<PARTICLES> MakeParticles(const std::string& name) {
    if (name == "CUBE_TO_CUBE_COLLISION") return std::make_unique<CUBE_TO_CUBE_COLLISION>();
    if (name == "RECTANGLE_FREEFALL_CYLINDER") return std::make_unique<RECTANGLE_FREEFALL_CYLINDER>();
    if (name == "RECTANGLE_FREEFALL_HAT_OBSTACLE") return std::make_unique<RECTANGLE_FREEFALL_HAT_OBSTACLE>();
    if (name == "CUBE_CRASH_PADDED_GROUND") return std::make_unique<CUBE_CRASH_PADDED_GROUND>();
    if (name == "CUBE_CRASH_WALL") return std::make_unique<CUBE_CRASH_WALL>();
    if (name == "CUBE_TO_CUBE_FREEFALL") return std::make_unique<CUBE_TO_CUBE_FREEFALL>();

    // NOTE: CUBE_FREEFALL is intentionally excluded. It's declared in
    // PARTICLES.h but has no method definitions anywhere in PARTICLES.cpp --
    // referencing `new CUBE_FREEFALL()` anywhere in the compiled binary
    // fails to link ("missing vtable"), which would break the whole build,
    // not just runs that select it. See scratch/REVIEW.md.
    throw std::runtime_error(
        "Unknown particles type \"" + name + "\". Valid options: "
        "CUBE_TO_CUBE_COLLISION, RECTANGLE_FREEFALL_CYLINDER, RECTANGLE_FREEFALL_HAT_OBSTACLE, "
        "CUBE_CRASH_PADDED_GROUND, CUBE_CRASH_WALL, CUBE_TO_CUBE_FREEFALL");
}

std::unique_ptr<SIMULATION_PARAMETERS> MakeSimulationParameters(const std::string& name) {
    if (name == "DEFAULT_PARAMETERS") return std::make_unique<DEFAULT_PARAMETERS>();
    if (name == "DEFAULT_IMPLICIT_PARAMETERS") return std::make_unique<DEFAULT_IMPLICIT_PARAMETERS>();
    if (name == "LOWER_YOUNGS_MODULUS") return std::make_unique<LOWER_YOUNGS_MODULUS>();
    if (name == "LOWER_YOUNGS_MODULUS_IMPLICIT") return std::make_unique<LOWER_YOUNGS_MODULUS_IMPLICIT>();
    if (name == "LOWER_CRITICAL_COMPRESSION_STRETCH_PARAMETERS")
        return std::make_unique<LOWER_CRITICAL_COMPRESSION_STRETCH_PARAMETERS>();
    if (name == "LOWER_HARDENING") return std::make_unique<LOWER_HARDENING>();
    if (name == "LOWER_CRITICAL_COMPRESSION_PARAMETERS") return std::make_unique<LOWER_CRITICAL_COMPRESSION_PARAMETERS>();
    if (name == "LOWER_CRITICAL_STRETCH_PARAMETERS") return std::make_unique<LOWER_CRITICAL_STRETCH_PARAMETERS>();
    if (name == "HYPERELASTICITY") return std::make_unique<HYPERELASTICITY>();
    if (name == "SCALED_DEFAULT_PARAMETERS") return std::make_unique<SCALED_DEFAULT_PARAMETERS>();
    if (name == "SCALED_HYPERELASTICITY") return std::make_unique<SCALED_HYPERELASTICITY>();

    throw std::runtime_error(
        "Unknown parameters type \"" + name + "\". Valid options: "
        "DEFAULT_PARAMETERS, DEFAULT_IMPLICIT_PARAMETERS, LOWER_YOUNGS_MODULUS, "
        "LOWER_YOUNGS_MODULUS_IMPLICIT, LOWER_CRITICAL_COMPRESSION_STRETCH_PARAMETERS, "
        "LOWER_HARDENING, LOWER_CRITICAL_COMPRESSION_PARAMETERS, LOWER_CRITICAL_STRETCH_PARAMETERS, "
        "HYPERELASTICITY, SCALED_DEFAULT_PARAMETERS, SCALED_HYPERELASTICITY");
}
