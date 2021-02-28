# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

module KrigingEstimators

using Meshes
using GeoStatsBase
using Variography

using LinearAlgebra: Factorization, Symmetric
using LinearAlgebra: bunchkaufman, cholesky, issuccess, ⋅
using Distances: Euclidean
using Distributions: Normal
using Combinatorics: multiexponents
using StaticArrays: MVector

import GeoStatsBase: fit, predict, status
import GeoStatsBase: solve

include("estimators.jl")
include("solvers/kriging.jl")

export
  # estimators
  KrigingEstimator,
  SimpleKriging,
  OrdinaryKriging,
  UniversalKriging,
  ExternalDriftKriging,
  fit, predict, status, weights,

  # solvers
  Kriging

end
