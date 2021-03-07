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

import GeoStatsBase: fit, predict, status
import GeoStatsBase: solve

include("estimators.jl")

export
  # estimators
  KrigingEstimator,
  SimpleKriging,
  OrdinaryKriging,
  UniversalKriging,
  ExternalDriftKriging,
  fit, predict, status, weights

end
