# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

module KrigingEstimators

using Meshes
using GeoTables
using GeoStatsBase
using Variography

using LinearAlgebra: Factorization, Symmetric
using LinearAlgebra: bunchkaufman, cholesky, issuccess, ⋅
using Combinatorics: multiexponents
using Distributions: Normal
using Unitful

import GeoStatsBase: fit, predict, predictprob, status

include("estimators.jl")

export
  # estimators
  KrigingEstimator,
  SimpleKriging,
  OrdinaryKriging,
  UniversalKriging,
  ExternalDriftKriging,
  variogram,
  fit,
  predict,
  status,
  weights,
  combine

end
