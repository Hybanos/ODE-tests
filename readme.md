N-body problem solver, using a few of the most common numerical integrators.

Currently implmented:
 - Euler
 - Symplectic Euler
 - LeapFrog
 - RK2
 - RK4
 - RK4(5)
 - DOP8(5,3)
 - Reference DOP8(5,3) implementation in Fortran [link](https://www.unige.ch/~hairer/prog/nonstiff/dop853.f)

https://github.com/user-attachments/assets/c05d0783-b1b9-4470-976a-01f7c282c02a

## Build

To build the simulation, run `./run` in the project directory. This should generate text files with the simulation results in place. To change the methods used in the simulation, edit `src/main.cpp`.

## Visualization

mp4 files of the simulations can be generated with [manim](https://www.manim.community/). It can be a pain to install locally, you best bet is to run `pip install manim` and follow the errors.