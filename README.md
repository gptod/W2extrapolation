# Extrapolation in the Wasserstein space

Computation of the extrapolation according to the quadratic Wasserstein distance between two probability measures discretized as point clouds. The code is an implementation of the algorithm cointained in:

[1] Gallouët, T. O., Natale, A., & Todeschi, G. (2025). [_Metric extrapolation in the Wasserstein space_](https://arxiv.org/abs/2407.10516). Calculus of Variations and Partial Differential Equations, 64(5), 147.

The algorithm is implemented in fortran and wrapped into a python module using f2py, for both performance and simplicity of use.

## Contents of the repository

The repository contains:
1. The fortran source code containing the implementation of the algorithm and the necessary subroutines: global.f90 and w2extrapolation.f90
2. The python interface function: extra_py.py
3. The python code to reproduce the main test from [1]: test_extra_n.py and plot_extra_n.py
4. The data to reproduce the abovementioned test: measures (with a smaller version for test measures_small)

## Use
To compile the fortran code and build the python wrapper, use the following command (requires installation of [f2py](https://numpy.org/doc/stable/f2py/))

python3 -m numpy.f2py --fcompiler=* -c global.f90 w2extrapolation.f90 -m w2extrapolation

where * is the fortran compiler to be used (e.g. gfortran).

The source global.f90 is only used to define the global precision for w2extrapolation.f90, single (for faster computation) or double. This can be easily changed but require then to recompile the fortran source code.

To redo the simulation from figure (2) in [1], run

python3 test_extra_n.py

to produce the extrapolation (file 'extrapolations'). Data and outputs are saved/loaded in binary format using the [pickle](https://docs.python.org/3/library/pickle.html) module. Run

python3 plot_extra_n.py

to plot the results. Create an adhoc folder where to save the images (e.g. 'images_extra_n') and specify the path in the code. 






