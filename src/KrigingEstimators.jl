# ------------------------------------------------------------------
# Copyright (c) 2018, Júlio Hoffimann Mendes <juliohm@stanford.edu>
# Licensed under the ISC License. See LICENCE in the project root.
# ------------------------------------------------------------------

__precompile__()

module KrigingEstimators

using Distances
using Combinatorics: multiexponents
using SpecialFunctions: besselk
using RecipesBase

# won't be neeeded in Julia v0.7
using Parameters

# extend result_type and pairwise for theoretical variograms
import Distances: result_type, pairwise

# variograms & estimators
include("variograms.jl")
include("estimators.jl")

# plot recipes
include("plotrecipes/variograms.jl")

export
  # variograms
  AbstractVariogram,
  GaussianVariogram,
  ExponentialVariogram,
  MaternVariogram,
  SphericalVariogram,
  CubicVariogram,
  PentasphericalVariogram,
  PowerVariogram,
  SineHoleVariogram,
  CompositeVariogram,
  isstationary,
  pairwise,

  # estimators
  KrigingEstimator,
  SimpleKriging,
  OrdinaryKriging,
  UniversalKriging,
  ExternalDriftKriging,
  fit!,
  weights,
  estimate

end
