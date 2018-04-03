# KrigingEstimators.jl

[![][travis-img]][travis-url] [![][julia-pkg-img]][julia-pkg-url] [![][codecov-img]][codecov-url]

This package provides high-performance implementations of Kriging estimators introduced by
[Matheron 1971](https://books.google.com/books/about/The_Theory_of_Regionalized_Variables_and.html?id=TGhGAAAAYAAJ).
As the most general form of estimation with covariance/variogram models, Kriging does **not** require distributional
assumptions, is well-defined in general Hilbert spaces, and is quite applicable to real world data.
Currently, the following Kriging variants are implemented:

- Simple Kriging
- Ordinary Kriging
- Universal Kriging
- External Drift Kriging

## Installation

Get the latest stable release with Julia's package manager:

```julia
Pkg.add("KrigingEstimators")
```

## Usage

This package is part of the [GeoStats.jl](https://github.com/juliohm/GeoStats.jl) framework.

For examples of usage, please check the main documentation.

## Asking for help

If you have any questions, please [open an issue](https://github.com/juliohm/KrigingEstimators.jl/issues).

[travis-img]: https://travis-ci.org/juliohm/KrigingEstimators.jl.svg?branch=master
[travis-url]: https://travis-ci.org/juliohm/KrigingEstimators.jl

[julia-pkg-img]: http://pkg.julialang.org/badges/KrigingEstimators_0.6.svg
[julia-pkg-url]: http://pkg.julialang.org/?pkg=KrigingEstimators

[codecov-img]: https://codecov.io/gh/juliohm/KrigingEstimators.jl/branch/master/graph/badge.svg
[codecov-url]: https://codecov.io/gh/juliohm/KrigingEstimators.jl
