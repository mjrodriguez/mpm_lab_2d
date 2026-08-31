#ifndef MPM_LAB_2D_CONFIG_H
#define MPM_LAB_2D_CONFIG_H

#include <memory>
#include <string>

#include "GRID.h"
#include "PARTICLES.h"
#include "SIMULATION_PARAMETERS.h"

// Everything needed to select and name a simulation run, loaded from a JSON
// config file (see LoadConfig). Replaces the old approach of hardcoding an
// output path and commenting/uncommenting class declarations in main.cpp.
struct RUN_CONFIG {
    std::string outputDirectory;
    std::string simulationName;
    std::string gridType;
    std::string particlesType;
    std::string parametersType;
};

// Reads and validates a JSON run config from `path`. Throws
// std::runtime_error with a clear, user-facing message if the file can't be
// opened, isn't valid JSON, or is missing a required key.
RUN_CONFIG LoadConfig(const std::string& path);

// Factories mapping a config's type name (e.g. "GRID_HAT_OBSTACLE") to a
// concrete instance. Each throws std::runtime_error listing the valid names
// if `name` isn't recognized.
std::unique_ptr<GRID> MakeGrid(const std::string& name);
std::unique_ptr<PARTICLES> MakeParticles(const std::string& name);
std::unique_ptr<SIMULATION_PARAMETERS> MakeSimulationParameters(const std::string& name);

#endif //MPM_LAB_2D_CONFIG_H
