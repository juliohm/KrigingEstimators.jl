# ------------------------------------------------------------------
# Copyright (c) 2018, Júlio Hoffimann Mendes <juliohm@stanford.edu>
# Licensed under the ISC License. See LICENCE in the project root.
# ------------------------------------------------------------------

__precompile__()

module KrigingEstimators

importall GeoStatsBase
using GeoStatsDevTools

using Combinatorics: multiexponents
using StatsBase: sample

include("estimators.jl")
include("solvers.jl")

export
  # estimators
  SimpleKriging,
  OrdinaryKriging,
  UniversalKriging,
  ExternalDriftKriging,
  fit!,
  weights,
  estimate,

  # solvers
  Kriging,
  SeqGaussSim

end
