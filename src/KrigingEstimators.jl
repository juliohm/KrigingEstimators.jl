# ------------------------------------------------------------------
# Licensed under the ISC License. See LICENCE in the project root.
# ------------------------------------------------------------------

module KrigingEstimators

using Reexport
using LinearAlgebra
using Combinatorics: multiexponents

# export variogram models
@reexport using Variography

include("estimators.jl")

@deprecate estimate predict
@deprecate fit! fit

export
  KrigingEstimator,
  SimpleKriging,
  OrdinaryKriging,
  UniversalKriging,
  ExternalDriftKriging,
  fit, status, predict, weights

end
