# ------------------------------------------------------------------
# Licensed under the ISC License. See LICENCE in the project root.
# ------------------------------------------------------------------

"""
    KrigingEstimator

A Kriging estimator (e.g. Simple Kriging).
"""
abstract type KrigingEstimator end

"""
    KrigingState(X, z, LHS, RHS)

A Kriging state stores information needed
to perform estimation at any given location.
"""
mutable struct KrigingState{T<:Real,V,
                            A<:AbstractMatrix{T},
                            B<:AbstractVector{V},
                            F<:Factorization,R}
  X::A
  z::B
  LHS::F
  RHS::Vector{R}
end

"""
    KrigingWeights(λ, ν)

An object storing Kriging weights `λ` and Lagrange multipliers `ν`.
"""
struct KrigingWeights{T<:Real,A<:AbstractVector{T}}
  λ::A
  ν::A
end

"""
    FittedKriging(estimator, state)

An object that can be used for making predictions using the
parameters in `estimator` and the current Kriging `state`.
"""
struct FittedKriging{E<:KrigingEstimator,S<:KrigingState}
  estimator::E
  state::S
end

"""
    status(fittedkrig)

Return the status of the `fittedkrig` object, meaning
the factorization of the Kriging system was successful.
"""
status(fittedkrig::FittedKriging) = issuccess(fittedkrig.state.LHS)

#--------------
# FITTING STEP
#--------------

"""
    fit(estimator, X, z)

Build Kriging system from coordinates `X` and
values `z` and return a fitted estimator.
"""
function fit(estimator::KrigingEstimator, X::AbstractMatrix, z::AbstractVector)
  # build Kriging system
  LHS = lhs(estimator, X)
  RHS = Vector{eltype(LHS)}(undef, size(LHS,1))

  # factorize LHS
  FLHS = factorize(estimator, LHS)

  # record Kriging state
  state = KrigingState(X, z, FLHS, RHS)

  # return fitted estimator
  FittedKriging(estimator, state)
end

"""
    lhs(estimator, X)

Return LHS of Kriging system using spatial configuration `X`.
"""
function lhs(estimator::KrigingEstimator, X::AbstractMatrix)
  γ = estimator.γ
  nobs = size(X, 2)
  ncons = nconstraints(estimator)

  # pre-allocate memory for LHS
  x = view(X,:,1)
  T = Variography.result_type(γ, x, x)
  m = nobs + ncons
  LHS = Matrix{T}(undef, m, m)

  # set variogram/covariance block
  pairwise!(LHS, γ, X)
  if isstationary(γ)
    for j=1:nobs, i=1:nobs
      @inbounds LHS[i,j] = sill(γ) - LHS[i,j]
    end
  end

  # set blocks of constraints
  set_constraints_lhs!(estimator, LHS, X)

  LHS
end

"""
    nconstraints(estimator)

Return number of constraints for `estimator`.
"""
nconstraints(estimator::KrigingEstimator) = error("not implemented")

"""
    set_constraints_lhs!(estimator, LHS, X)

Set constraints in LHS of Kriging system.
"""
set_constraints_lhs!(estimator::KrigingEstimator,
                     LHS::AbstractMatrix, X::AbstractMatrix) = error("not implemented")

"""
    factorize(estimator, LHS)

Factorize LHS of Kriging system with appropriate factorization method.
"""
factorize(estimator::KrigingEstimator, LHS::AbstractMatrix) = error("not implemented")

#-----------------
# PREDICTION STEP
#-----------------

"""
    predict(estimator, xₒ)

Compute mean and variance for the `estimator` at coordinates `xₒ`.
"""
predict(estimator::FittedKriging, xₒ::AbstractVector) =
  combine(estimator, weights(estimator, xₒ), estimator.state.z)

"""
    weights(estimator, xₒ)

Compute the weights λ (and Lagrange multipliers ν) for the
`estimator` at coordinates `xₒ`.
"""
function weights(estimator::FittedKriging, xₒ::AbstractVector)
  nobs = size(estimator.state.X, 2)

  set_rhs!(estimator, xₒ)

  # solve Kriging system
  x = estimator.state.LHS \ estimator.state.RHS

  λ = view(x,1:nobs)
  ν = view(x,nobs+1:length(x))

  KrigingWeights(λ, ν)
end

"""
    set_rhs!(estimator, xₒ)

Set RHS of Kriging system at coodinates `xₒ`.
"""
function set_rhs!(estimator::FittedKriging, xₒ::AbstractVector)
  γ = estimator.estimator.γ
  X = estimator.state.X
  RHS = estimator.state.RHS

  # RHS variogram/covariance
  @inbounds for j in 1:size(X, 2)
    xj = view(X,:,j)
    RHS[j] = isstationary(γ) ? sill(γ) - γ(xj, xₒ) : γ(xj, xₒ)
  end

  set_constraints_rhs!(estimator, xₒ)
end

"""
    set_constraints_rhs!(estimator, xₒ)

Set constraints in RHS of Kriging system.
"""
set_constraints_rhs!(estimator::FittedKriging, xₒ::AbstractVector) = error("not implemented")

"""
    combine(estimator, weights, z)

Combine `weights` with values `z` to produce mean and variance.
"""
function combine(estimator::FittedKriging, weights::KrigingWeights, z::AbstractVector)
  γ = estimator.estimator.γ
  b = estimator.state.RHS
  λ = weights.λ
  ν = weights.ν

  # compute b⋅[λ;ν]
  nobs  = length(λ)
  c₁ = view(b,1:nobs)⋅λ
  c₂ = view(b,nobs+1:length(b))⋅ν
  c = c₁ + c₂

  if isstationary(γ)
    z⋅λ, sill(γ) - c
  else
    z⋅λ, c
  end
end

#-----------------
# IMPLEMENTATIONS
#-----------------
include("estimators/simple_kriging.jl")
include("estimators/ordinary_kriging.jl")
include("estimators/universal_kriging.jl")
include("estimators/external_drift_kriging.jl")
