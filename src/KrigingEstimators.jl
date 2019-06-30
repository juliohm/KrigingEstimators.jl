# ------------------------------------------------------------------
# Licensed under the ISC License. See LICENCE in the project root.
# ------------------------------------------------------------------

module KrigingEstimators

using GeoStatsBase
using Variography

using LinearAlgebra: Factorization, lu, cholesky, issuccess, ⋅
using Distances: Euclidean
using Distributions: Normal
using Combinatorics: multiexponents
using StaticArrays: MVector

import GeoStatsBase

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
  weights,

  # solvers
  Kriging, KrigingParam,
  SeqGaussSim, SeqGaussSimParam

end
