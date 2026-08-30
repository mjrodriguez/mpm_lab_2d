# mpm_lab_2d
[![CI](https://github.com/mjrodriguez/mpm_lab_2d/actions/workflows/ci.yml/badge.svg)](https://github.com/mjrodriguez/mpm_lab_2d/actions/workflows/ci.yml)

A 2D implementation of the Material Point Method (MPM) constitutive model and time integration scheme described in:

> Alexey Stomakhin, Craig Schroeder, Lawrence Chai, Joseph Teran, and Andrew Selle. 2013. **A material point method for snow simulation**. *ACM Trans. Graph.* 32, 4, Article 102 (July 2013), 12 pages. https://doi.org/10.1145/2461912.2461948

This code uses the paper's fixed-corotated hyperelastic model coupled with plastic deformation determined by the principal stretches, its cubic B-spline interpolation kernel, and its explicit time integration scheme. It is implemented in two dimensions only — the paper's method is 3D — and does not include the rendering pipeline. Many of the details can be found in my MS Thesis report, feel free to contact me if you'd like a copy.


# Compiling:
This project uses **CMake** and depends on the **Eigen** linear algebra library. You don't need to install anything up front: CMake will use a system-installed Eigen (`3.3+`) if it finds one, and otherwise will automatically download a pinned version (`5.0.1`) at configure time.

```
mkdir build && cd build
cmake ..
cmake --build .
```

The resulting `mpm` executable will be in the `build` directory.

To skip the automatic download and speed up configuring, install Eigen with your package manager first, e.g. `brew install eigen` (macOS) or `sudo apt install libeigen3-dev` (Debian/Ubuntu).

# Testing:
Unit tests use **Catch2**, fetched automatically by CMake. From the `build` directory:

```
ctest
```
