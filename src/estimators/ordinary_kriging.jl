# ------------------------------------------------------------------
# Licensed under the ISC License. See LICENCE in the project root.
# ------------------------------------------------------------------

"""
    OrdinaryKriging(γ)
    OrdinaryKriging(X, z, γ)

Ordinary Kriging with variogram model `γ`.

Optionally, pass the coordinates `X` and values `z`
to the [`fit`](@ref) function.
"""
struct OrdinaryKriging{G<:Variogram} <: KrigingEstimator
  γ::G
end

OrdinaryKriging(γ) = OrdinaryKriging{typeof(γ)}(γ)

OrdinaryKriging(X, z, γ) = fit(OrdinaryKriging(γ), X, z)

nconstraints(estimator::OrdinaryKriging) = 1

function set_constraints_lhs!(estimator::OrdinaryKriging, LHS::AbstractMatrix, X::AbstractMatrix)
  T = eltype(LHS)
  LHS[end,:]   .= one(T)
  LHS[:,end]   .= one(T)
  LHS[end,end]  = zero(T)

  nothing
end

factorize(estimator::OrdinaryKriging, LHS::AbstractMatrix) = lu(LHS, check=false)

function set_constraints_rhs!(estimator::FittedKriging{E,S},
                              xₒ::AbstractVector) where {E<:OrdinaryKriging,S<:KrigingState}
  RHS = estimator.state.RHS

  RHS[end] = one(eltype(RHS))

  nothing
end
