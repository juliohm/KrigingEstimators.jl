# KrigingEstimators.jl

[![][travis-img]][travis-url] [![][julia-pkg-img]][julia-pkg-url] [![][codecov-img]][codecov-url]

This package provides high-performance implementations of Kriging estimators introduced by
[Matheron 1971](https://books.google.com/books/about/The_Theory_of_Regionalized_Variables_and.html?id=TGhGAAAAYAAJ),
and Kriging-based solvers (e.g. sequential Gaussian simulation) for the [GeoStats.jl](https://github.com/juliohm/GeoStats.jl)
framework.

Currently, the following Kriging variants are implemented:

- Simple Kriging
- Ordinary Kriging
- Universal Kriging
- External Drift Kriging

Unlike most popular estimators in statistics, Kriging does **not** assume independent and identically distributed residuals.
No distributional assumptions are required in the derivation of Kriging estimators, which makes these methods quite useful
in problems with real world data.

## Installation

Get the latest stable release with Julia's package manager:

```julia
Pkg.add("KrigingEstimators")
```

## Usage

This package is part of the [GeoStats.jl](https://github.com/juliohm/GeoStats.jl) framework.

For a simple example of usage, please check [this notebook](docs/Usage.ipynb).

## Asking for help

If you have any questions, please [open an issue](https://github.com/juliohm/KrigingEstimators.jl/issues).

[travis-img]: https://travis-ci.org/juliohm/KrigingEstimators.jl.svg?branch=master
[travis-url]: https://travis-ci.org/juliohm/KrigingEstimators.jl

[julia-pkg-img]: http://pkg.julialang.org/badges/KrigingEstimators_0.6.svg
[julia-pkg-url]: http://pkg.julialang.org/?pkg=KrigingEstimators

[codecov-img]: https://codecov.io/gh/juliohm/KrigingEstimators.jl/branch/master/graph/badge.svg
[codecov-url]: https://codecov.io/gh/juliohm/KrigingEstimators.jl
