# mpm_lab_2d
Material point method in 2D. Based off the paper from Stomakhin et al.

This code uses a hyperelastic model coupled with plastic deformation determined by the principal stretches. The code currently only includes an explicit time integration scheme and uses a cubic B-Spline. The code also is implemented in two-dimensions but an extension to three-dimensions is in process. Many of the details can be found in my MS Thesis report, feel free to contact me if you'd like a copy.


# Compiling:
This project uses **CMake**. Since the **Eigen Linear Algebra** package is just header files, no need to install or compile it separately.

```
mkdir build && cd build
cmake ..
cmake --build .
```

The resulting `mpm` executable will be in the `build` directory.

# Testing:
Unit tests use **Catch2**, fetched automatically by CMake. From the `build` directory:

```
ctest
```
