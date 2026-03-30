This repo contains a few of the most common methods for solving the n-body problem.

Currently implmented:
 - Exact analytical solution (only for 2 bodies)
 - Euler
 - Symplectic Euler
 - LeapFrog
 - RK2
 - RK4

<!-- ![haha](.github/media/4-body-cmp.mp4) -->

## Build

To build the simulation, run `./run` in the project directory. This should generate text files with the simulation results in place. To change the methods used in the simulation, edit `src/main.cpp`.

## Visualization

mp4 files of the simulations can be generated with [manim](https://www.manim.community/). It can be a pain to install locally, you best bet is to run `pip install manim` and follow the errors.