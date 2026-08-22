# mpm_lab_2d
A 2D implementation of the Material Point Method (MPM) constitutive model and time integration scheme described in:

> Alexey Stomakhin, Craig Schroeder, Lawrence Chai, Joseph Teran, and Andrew Selle. 2013. **A material point method for snow simulation**. *ACM Trans. Graph.* 32, 4, Article 102 (July 2013), 12 pages. https://doi.org/10.1145/2461912.2461948

This code uses the paper's fixed-corotated hyperelastic model coupled with plastic deformation determined by the principal stretches, its cubic B-spline interpolation kernel, and its explicit time integration scheme. It is implemented in two dimensions only — the paper's method is 3D — and does not include the rendering pipeline. Many of the details can be found in my MS Thesis report, feel free to contact me if you'd like a copy.


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
