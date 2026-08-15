# mpm_lab_2d
Material point method in 2D. Based off the paper from Stomakhin et al.

This code uses a hyperelastic model coupled with plastic deformation determined by the principal stretches. The code currently only includes an explicit time integration scheme and uses a cubic B-Spline. The code also is implemented in two-dimensions but an extension to three-dimensions is in process. Many of the details can be found in my MS Thesis report, feel free to contact me if you'd like a copy.


# Compiling:
This project uses **CMake** and depends on the **Eigen** linear algebra library. You don't need to install anything up front: CMake will use a system-installed Eigen (`3.3+`) if it finds one, and otherwise will automatically download a pinned version (`3.4.0`) at configure time.

```
mkdir build && cd build
cmake ..
cmake --build .
```

The resulting `mpm` executable will be in the `build` directory.

To skip the automatic download and speed up configuring, install Eigen with your package manager first, e.g. `brew install eigen` (macOS) or `sudo apt install libeigen3-dev` (Debian/Ubuntu).
