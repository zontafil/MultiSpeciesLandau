# Multi specie structure preserving discrete Landau operator

## Install

This project requires the Eigen library (https://eigen.tuxfamily.org/). It is supposed to be downloaded and installed locally in `./Eigen/` folder.

Alternatively, Eigen can be installed globally, and `coulomb_utils.h` needs to be changed.

A better way to integrate Eigen is needed.

## Configure

Copy an example simulation from `examples/` to `main.cpp`

## Compile

`CUDA=1 make`

### Logging levels and debug symbols

Set the ENV variables `VERBOSE=1` or `SILLY=1` to have more verbose debug logs.

Set `DEBUG=1` to enable the debug symbols and disable memory optimizations.

## Run

`./coulomb`

The output files are put in `out/data`. For every n-th timestep the distribution functions are plotted along with some other quantities, i.e. Energy, temperature, momentum.

## Scripts

`python plot_distribution.py`: Plot from the output files. The plots are put in `out/`

`saveSim.sh [LABEL]`: Save the configuration, plots and out data to `examples/`, `plots/`

`deleteSim.sh example\{file}.cpp`: Delete a saved simulation






