# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

module KrigingEstimators

using Meshes
using GeoStatsBase
using Variography

using LinearAlgebra: Factorization, Symmetric
using LinearAlgebra: bunchkaufman, cholesky, issuccess, ⋅
using Combinatorics: multiexponents
using Distributions: Normal
using Unitful

import GeoStatsBase: fit, predict, predictprob, status
import GeoStatsBase: solve

include("estimators.jl")
include("ui.jl")

export
  # estimators
  KrigingEstimator,
  SimpleKriging,
  OrdinaryKriging,
  UniversalKriging,
  ExternalDriftKriging,
  fit,
  predict,
  status,
  weights,
  combine,

  # UI elements
  kriging_ui

end
