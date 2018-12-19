# ------------------------------------------------------------------
# Licensed under the ISC License. See LICENCE in the project root.
# ------------------------------------------------------------------

"""
    KrigingEstimator

A Kriging estimator (e.g. Simple Kriging).
"""
abstract type KrigingEstimator end

"""
    fit!(estimator, X, z)

Build LHS of Kriging system from coordinates `X` with
values `z` and save factorization in `estimator`.
"""
function fit!(estimator::KrigingEstimator, X::AbstractMatrix, z::AbstractVector)
  # update data
  estimator.X = X
  estimator.z = z

  # build and factorize LHS
  status = set_lhs!(estimator, X)

  # pre-allocate memory for RHS
  T = eltype(estimator.LHS)
  m = size(estimator.LHS, 1)
  estimator.RHS = Vector{T}(undef, m)

  status
end

"""
    estimate(estimator, xₒ)

Compute mean and variance for the `estimator` at coordinates `xₒ`.
"""
estimate(estimator::KrigingEstimator, xₒ::AbstractVector) =
  combine(estimator, weights(estimator, xₒ), estimator.z)

"""
    weights(estimator, xₒ)

Compute the weights λ (and Lagrange multipliers ν) for the
`estimator` at coordinates `xₒ`.
"""
function weights(estimator::KrigingEstimator, xₒ::AbstractVector)
  nobs = length(estimator.z)

  # build RHS
  set_rhs!(estimator, xₒ)

  # solve Kriging system
  x = estimator.LHS \ estimator.RHS

  # return weights
  Weights(x[1:nobs], x[nobs+1:end])
end

"""
    set_lhs!(estimator, X)

Set LHS of Kriging system using spatial configuration `X`.
"""
function set_lhs!(estimator::KrigingEstimator, X::AbstractMatrix)
  γ = estimator.γ
  nobs = length(estimator.z)
  ncons = nconstraints(estimator)

  # pre-allocate memory for LHS
  x = view(X, :, 1)
  T = Variography.result_type(γ, x, x)
  m = nobs + ncons
  LHS = Matrix{T}(undef, m, m)

  # set variogram/covariance block
  pairwise!(LHS, γ, X)
  if isstationary(γ)
    for j=1:nobs, i=1:nobs
      LHS[i,j] = sill(γ) - LHS[i,j]
    end
  end

  # set blocks of constraints
  set_constraints_lhs!(estimator, LHS)

  # factorize LHS and save
  estimator.LHS = factorize(estimator, LHS)

  # return factorization status
  issuccess(estimator.LHS)
end

"""
    set_rhs!(estimator, xₒ)

Set RHS of Kriging system at coodinates `xₒ`.
"""
function set_rhs!(estimator::KrigingEstimator, xₒ::AbstractVector)
  X = estimator.X
  γ = estimator.γ

  # RHS variogram/covariance
  RHS = estimator.RHS
  for j in 1:size(X, 2)
    xj = view(X, :, j)
    RHS[j] = isstationary(γ) ? sill(γ) - γ(xj, xₒ) : γ(xj, xₒ)
  end

  set_constraints_rhs!(estimator, xₒ)
end

"""
    nconstraints(estimator)

Return number of constraints for `estimator`.
"""
nconstraints(estimator::KrigingEstimator) = error("not implemented")

"""
    set_constraints_lhs!(estimator, LHS)

Set constraints in LHS of Kriging system.
"""
set_constraints_lhs!(estimator::KrigingEstimator, LHS::AbstractMatrix) = error("not implemented")

"""
    set_constraints_rhs!(estimator, xₒ)

Set constraints in RHS of Kriging system.
"""
set_constraints_rhs!(estimator::KrigingEstimator, xₒ::AbstractVector) = error("not implemented")

"""
    factorize(estimator, LHS)

Factorize LHS of Kriging system with appropriate factorization method.
"""
factorize(estimator::KrigingEstimator, LHS::AbstractMatrix) = error("not implemented")

"""
    Weights(λ, ν)

An object storing Kriging weights `λ` and Lagrange multipliers `ν`.
"""
struct Weights{T<:Real}
  λ::Vector{T}
  ν::Vector{T}
end

"""
    combine(estimator, weights, z)

Combine `weights` with values `z` to produce mean and variance.
"""
function combine(estimator::KrigingEstimator, weights::Weights, z::AbstractVector)
  γ = estimator.γ
  b = estimator.RHS
  λ = weights.λ
  ν = weights.ν

  if isstationary(γ)
    z⋅λ, sill(γ) - b⋅[λ;ν]
  else
    z⋅λ, b⋅[λ;ν]
  end
end

#------------------
# IMPLEMENTATIONS
#------------------
include("estimators/simple_kriging.jl")
include("estimators/ordinary_kriging.jl")
include("estimators/universal_kriging.jl")
include("estimators/external_drift_kriging.jl")
