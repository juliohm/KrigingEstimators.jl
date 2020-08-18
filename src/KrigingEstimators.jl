# ------------------------------------------------------------------
# Licensed under the ISC License. See LICENSE in the project root.
# ------------------------------------------------------------------

module KrigingEstimators

using GeoStatsBase
using Variography

using LinearAlgebra: Factorization, Symmetric
using LinearAlgebra: lu, cholesky, issuccess, ⋅
using Distances: Euclidean
using Distributions: Normal
using Combinatorics: multiexponents
using StaticArrays: MVector

import GeoStatsBase: fit, predict, status
import GeoStatsBase: preprocess, solve, solvesingle

include("estimators.jl")

include("solvers/kriging.jl")
include("solvers/sgsim.jl")

export
  # estimators
  KrigingEstimator,
  SimpleKriging,
  OrdinaryKriging,
  UniversalKriging,
  ExternalDriftKriging,
  fit, predict, status, weights,

  # solvers
  Kriging, KrigingParam,
  SeqGaussSim, SeqGaussSimParam

end
