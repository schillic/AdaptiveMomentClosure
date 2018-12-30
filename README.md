This is a tool for the identification of unknown rate parameters of chemical reaction networks.

The approach was presented in the paper [*Adaptive moment closure for parameter inference of biochemical reaction networks*](https://doi.org/10.1016/j.biosystems.2016.07.005) in BioSystems 2016.

# INSTALLATION

You need to have Matlab.

We use the following libraries:
- [StochDynTools](https://www.ece.ucsb.edu/~hespanha/software/stochdyntool.html)
- [SUNDIALS](http://computation.llnl.gov/casc/sundials/main.html) which provides the *CVODE* an optional but much faster solver of ordinary differential equations.

We deploy The libraries, but it may be necessary to install them anyway.

# RUN

## Initialization

Open Matlab.
Run `init` once per session.

## Running a model or benchmarks

There are some model scripts in the folder `examples/`.
Just run them with `runX`, where `X` is the model/folder name.

There are also some benchmark scripts in the folder `functions/benchmarks/`.
