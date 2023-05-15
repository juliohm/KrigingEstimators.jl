# KrigingEstimators.jl

[![][build-img]][build-url] [![][codecov-img]][codecov-url]

This package provides high-performance implementations of Kriging estimators introduced by
[Matheron 1971](https://books.google.com/books/about/The_Theory_of_Regionalized_Variables_and.html?id=TGhGAAAAYAAJ),
and further developed by various geostatisticians. As the most general form of estimation with covariance/variogram
models, Kriging does **not** require distributional assumptions, is well-defined in general Hilbert spaces, and is
unbiased. These properties make Kriging a quite useful method for estimation with real world data.
Currently, the following Kriging variants are implemented:

- Simple Kriging
- Ordinary Kriging
- Universal Kriging
- External Drift Kriging

## Installation

Get the latest stable release with Julia's package manager:

```julia
] add KrigingEstimators
```

## Usage

This package is part of the [GeoStats.jl](https://github.com/JuliaEarth/GeoStats.jl) framework.

For examples of usage, please check the main documentation.

## Asking for help

If you have any questions, please [contact our community](https://juliaearth.github.io/GeoStats.jl/stable/about/community.html).

[build-img]: https://img.shields.io/github/actions/workflow/status/JuliaEarth/KrigingEstimators.jl/CI.yml?branch=master&style=flat-square
[build-url]: https://github.com/JuliaEarth/KrigingEstimators.jl/actions

[codecov-img]: https://img.shields.io/codecov/c/github/JuliaEarth/KrigingEstimators.jl?style=flat-square
[codecov-url]: https://codecov.io/gh/JuliaEarth/KrigingEstimators.jl
