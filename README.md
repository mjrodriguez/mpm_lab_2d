# mpm_lab_2d
[![CI](https://github.com/mjrodriguez/mpm_lab_2d/actions/workflows/ci.yml/badge.svg)](https://github.com/mjrodriguez/mpm_lab_2d/actions/workflows/ci.yml)

A 2D implementation of the Material Point Method (MPM) constitutive model and time integration scheme described in:

> Alexey Stomakhin, Craig Schroeder, Lawrence Chai, Joseph Teran, and Andrew Selle. 2013. **A material point method for snow simulation**. *ACM Trans. Graph.* 32, 4, Article 102 (July 2013), 12 pages. https://doi.org/10.1145/2461912.2461948

This code uses the paper's fixed-corotated hyperelastic model coupled with plastic deformation determined by the principal stretches, its cubic B-spline interpolation kernel, and its explicit time integration scheme. It is implemented in two dimensions only — the paper's method is 3D — and does not include the rendering pipeline. Many of the details can be found in my MS Thesis report, feel free to contact me if you'd like a copy.


# Compiling:
This project uses **CMake** and depends on the **Eigen** linear algebra library and **nlohmann/json**. You don't need to install anything up front: CMake will use system-installed versions if it finds them (Eigen `3.3+`, nlohmann/json `3.11+`), and otherwise will automatically download pinned versions at configure time.

```
mkdir build && cd build
cmake ..
cmake --build .
```

The resulting `mpm` executable will be in the `build` directory.

To skip the automatic download and speed up configuring, install these with your package manager first, e.g. `brew install eigen nlohmann-json` (macOS) or `sudo apt install libeigen3-dev nlohmann-json3-dev` (Debian/Ubuntu).

# Running:
`mpm` takes a single argument: the path to a JSON config file describing what to run.

```
./mpm ../configs/default.json
```

A config file specifies the output directory and which scenario/parameter classes to run:

```json
{
  "output_directory": "./output",
  "simulation_name": "default_run",
  "grid": "GRID_HAT_OBSTACLE",
  "particles": "RECTANGLE_FREEFALL_HAT_OBSTACLE",
  "parameters": "DEFAULT_PARAMETERS"
}
```

- `grid`: `GRID_NO_OBSTACLE`, `GRID_CYLINDER_OBSTACLE`, `GRID_HAT_OBSTACLE`
- `particles`: `CUBE_TO_CUBE_COLLISION`, `RECTANGLE_FREEFALL_CYLINDER`, `RECTANGLE_FREEFALL_HAT_OBSTACLE`, `CUBE_CRASH_PADDED_GROUND`, `CUBE_CRASH_WALL`, `CUBE_TO_CUBE_FREEFALL`
- `parameters`: `DEFAULT_PARAMETERS`, `DEFAULT_IMPLICIT_PARAMETERS`, `LOWER_YOUNGS_MODULUS`, `LOWER_YOUNGS_MODULUS_IMPLICIT`, `LOWER_CRITICAL_COMPRESSION_STRETCH_PARAMETERS`, `LOWER_HARDENING`, `LOWER_CRITICAL_COMPRESSION_PARAMETERS`, `LOWER_CRITICAL_STRETCH_PARAMETERS`, `HYPERELASTICITY`, `SCALED_DEFAULT_PARAMETERS`, `SCALED_HYPERELASTICITY`

Passing an unrecognized name for any of these prints the valid options and exits. `output_directory` must already exist (it isn't created automatically) and is resolved relative to wherever `mpm` is run from, not relative to the config file — that's why `configs/default.json` uses `../output` (it's meant to be run from `build/`, matching the example above).

# Testing:
Unit tests use **Catch2**, fetched automatically by CMake. From the `build` directory:

```
ctest
```
