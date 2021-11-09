This is a tool for the identification of unknown rate parameters of chemical reaction networks.

The approach was presented in the article [*Adaptive moment closure for parameter inference of biochemical reaction networks*](https://doi.org/10.1016/j.biosystems.2016.07.005) in BioSystems, 2016.
See [below](#citation) for how to cite this work.

# Installation

You need to have Matlab.

We use the following libraries:
- [StochDynTools](https://www.ece.ucsb.edu/~hespanha/software/stochdyntool.html)
- [SUNDIALS](http://computation.llnl.gov/casc/sundials/main.html) (which provides *CVODE*, an optional but much faster solver of ordinary differential equations)

We deploy the libraries, but it may be necessary to install them anyway.

# Run

## Initialization

Open Matlab.
Run `init` at the beginning of each session.

## Running a model or benchmarks

There are some model scripts in the folder `examples/`.
Just run them with `runX`, where `X` is the model/folder name.

There are also some benchmark scripts in the folder `functions/benchmarks/`.

# Citation

```bibtex
@article{SchillingBHPR16,
  author    = {Christian Schilling and
               Sergiy Bogomolov and
               Thomas A. Henzinger and
               Andreas Podelski and
               Jakob Ruess},
  title     = {Adaptive moment closure for parameter inference of biochemical reaction
               networks},
  journal   = {Biosyst.},
  volume    = {149},
  pages     = {15--25},
  year      = {2016},
  url       = {https://doi.org/10.1016/j.biosystems.2016.07.005},
  doi       = {10.1016/j.biosystems.2016.07.005}
}
```
