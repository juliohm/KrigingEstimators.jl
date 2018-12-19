# ------------------------------------------------------------------
# Licensed under the ISC License. See LICENCE in the project root.
# ------------------------------------------------------------------

"""
    OrdinaryKriging(X, z, γ)

## Parameters

* X ∈ ℜ^(mxn) - matrix of data locations
* z ∈ ℜⁿ      - vector of observations for X
* γ           - variogram model
"""
mutable struct OrdinaryKriging{T<:Real,V} <: KrigingEstimator
  # input fields
  γ::Variogram

  # state fields
  X::Matrix{T}
  z::Vector{V}
  LHS::Factorization
  RHS::Vector

  function OrdinaryKriging{T,V}(γ; X=nothing, z=nothing) where {T<:Real,V}
    OK = new(γ)
    if X ≠ nothing && z ≠ nothing
      fit!(OK, X, z)
    end

    OK
  end
end

OrdinaryKriging(X, z, γ) = OrdinaryKriging{eltype(X),eltype(z)}(γ, X=X, z=z)

nconstraints(estimator::OrdinaryKriging) = 1

function set_constraints_lhs!(estimator::OrdinaryKriging, LHS::AbstractMatrix)
  T = eltype(LHS)
  LHS[end,:]   .= one(T)
  LHS[:,end]   .= one(T)
  LHS[end,end]  = zero(T)

  nothing
end

function set_constraints_rhs!(estimator::OrdinaryKriging, xₒ::AbstractVector)
  RHS = estimator.RHS
  RHS[end] = one(eltype(RHS))

  nothing
end

factorize(estimator::OrdinaryKriging, LHS::AbstractMatrix) = lu(LHS, check=false)
