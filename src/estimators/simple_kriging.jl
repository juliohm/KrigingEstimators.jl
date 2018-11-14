# ------------------------------------------------------------------
# Copyright (c) 2017, Júlio Hoffimann Mendes <juliohm@stanford.edu>
# Licensed under the ISC License. See LICENCE in the project root.
# ------------------------------------------------------------------

"""
    SimpleKriging(X, z, γ, μ)

## Parameters

* X ∈ ℜ^(mxn) - matrix of data locations
* z ∈ ℜⁿ      - vector of observations for X
* γ           - variogram model
* μ ∈ ℜ       - mean of z

### Notes

* Simple Kriging requires stationary variograms
"""
mutable struct SimpleKriging{T<:Real,V} <: KrigingEstimator
  # input fields
  γ::Variogram
  μ::V

  # state fields
  X::Matrix{T}
  z::Vector{V}
  LHS::Factorization
  RHS::Vector

  function SimpleKriging{T,V}(γ, μ; X=nothing, z=nothing) where {T<:Real,V}
    @assert isstationary(γ) "Simple Kriging requires stationary variogram"
    SK = new(γ, μ)
    if X ≠ nothing && z ≠ nothing
      fit!(SK, X, z)
    end

    SK
  end
end

SimpleKriging(X, z, γ, μ) = SimpleKriging{eltype(X),eltype(z)}(γ, μ, X=X, z=z)

nconstraints(estimator::SimpleKriging) = 0

set_constraints_lhs!(estimator::SimpleKriging, LHS::AbstractMatrix) = nothing

set_constraints_rhs!(estimator::SimpleKriging, xₒ::AbstractVector) = nothing

factorize(estimator::SimpleKriging, LHS::AbstractMatrix) = cholesky(LHS)

function combine(estimator::SimpleKriging{T,V},
                 weights::Weights, z::AbstractVector) where {T<:Real,V}
  γ = estimator.γ
  μ = estimator.μ
  b = estimator.RHS
  λ = weights.λ
  y = z .- μ

  μ + y⋅λ, sill(γ) - b⋅λ
end
