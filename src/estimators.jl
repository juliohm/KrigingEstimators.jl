# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    KrigingEstimator

A Kriging estimator (e.g. Simple Kriging).
"""
abstract type KrigingEstimator end

"""
    KrigingState(data, LHS, RHS)

A Kriging state stores information needed
to perform estimation at any given location.
"""
mutable struct KrigingState{D<:Data,F<:Factorization,R}
  data::D
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
    fit(estimator, data)

Build Kriging system from `data` and return a fitted estimator.
"""
function fit(estimator::KrigingEstimator, data)
  # build Kriging system
  LHS = lhs(estimator, domain(data))
  RHS = Vector{eltype(LHS)}(undef, size(LHS,1))

  # factorize LHS
  FLHS = factorize(estimator, LHS)

  # record Kriging state
  state = KrigingState(data, FLHS, RHS)

  # return fitted estimator
  FittedKriging(estimator, state)
end

"""
    lhs(estimator, domain)

Return LHS of Kriging system for the elements in the `domain`.
"""
function lhs(estimator::KrigingEstimator, domain)
  γ = estimator.γ
  nobs = nelements(domain)
  ncon = nconstraints(estimator)

  # pre-allocate memory for LHS
  u = first(domain)
  R = Variography.result_type(γ, u, u)
  m = nobs + ncon
  LHS = Matrix{R}(undef, m, m)

  # set variogram/covariance block
  pairwise!(LHS, γ, domain)
  if isstationary(γ)
    σ² = sill(γ)
    for j in 1:nobs, i in 1:nobs
      @inbounds LHS[i,j] = σ² - LHS[i,j]
    end
  end

  # set blocks of constraints
  set_constraints_lhs!(estimator, LHS, domain)

  LHS
end

"""
    nconstraints(estimator)

Return number of constraints for `estimator`.
"""
function nconstraints end

"""
    set_constraints_lhs!(estimator, LHS, X)

Set constraints in LHS of Kriging system.
"""
function set_constraints_lhs! end

"""
    factorize(estimator, LHS)

Factorize LHS of Kriging system with appropriate factorization method.
"""
function factorize end

#-----------------
# PREDICTION STEP
#-----------------

"""
    predict(estimator, var, uₒ)

Compute mean and variance of variable `var` using the
`estimator` at point or geometry `uₒ`.
"""
function predict(fitted::FittedKriging, var, uₒ)
  data = fitted.state.data
  combine(fitted, weights(fitted, uₒ), data[var])
end

"""
    weights(estimator, uₒ)

Compute the weights λ (and Lagrange multipliers ν) for the
`estimator` at point or geometry `uₒ`.
"""
function weights(fitted::FittedKriging, uₒ)
  nobs = nelements(fitted.state.data)

  set_rhs!(fitted, uₒ)

  # solve Kriging system
  s = fitted.state.LHS \ fitted.state.RHS

  λ = view(s, 1:nobs)
  ν = view(s, nobs+1:length(s))

  KrigingWeights(λ, ν)
end

"""
    set_rhs!(estimator, uₒ)

Set RHS of Kriging system at point or geometry `uₒ`.
"""
function set_rhs!(fitted::FittedKriging, uₒ)
  γ = fitted.estimator.γ
  dom = domain(fitted.state.data)
  nel = nelements(dom)
  RHS = fitted.state.RHS

  # RHS variogram/covariance
  RHS[1:nel] .= map(u -> γ(u, uₒ), dom)
  if isstationary(γ)
    RHS[1:nel] .= sill(γ) .- RHS[1:nel]
  end

  set_constraints_rhs!(fitted, uₒ)
end

"""
    set_constraints_rhs!(estimator, xₒ)

Set constraints in RHS of Kriging system.
"""
function set_constraints_rhs! end

"""
    combine(estimator, weights, z)

Combine `weights` with values `z` to produce mean and variance
using the appropriate formulas for the `estimator`.
"""
function combine(fitted::FittedKriging,
                 weights::KrigingWeights,
                 z::AbstractVector)
  γ = fitted.estimator.γ
  b = fitted.state.RHS
  λ = weights.λ
  ν = weights.ν

  # compute b⋅[λ;ν]
  nobs  = length(λ)
  c₁ = view(b, 1:nobs) ⋅ λ
  c₂ = view(b, nobs+1:length(b)) ⋅ ν
  c = c₁ + c₂

  if isstationary(γ)
    z⋅λ, sill(γ) - c
  else
    z⋅λ, c
  end
end

# ----------------
# IMPLEMENTATIONS
# ----------------

include("estimators/sk.jl")
include("estimators/ok.jl")
include("estimators/uk.jl")
include("estimators/edk.jl")
