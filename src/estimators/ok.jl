# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    OrdinaryKriging(γ)
    OrdinaryKriging(data, γ)

Ordinary Kriging with variogram model `γ`.

Optionally, pass the geospatial `data` to the [`fit`](@ref) function.
"""
struct OrdinaryKriging{G<:Variogram} <: KrigingEstimator
  γ::G
end

OrdinaryKriging(γ) = OrdinaryKriging{typeof(γ)}(γ)

OrdinaryKriging(data, γ) = GeoStatsBase.fit(OrdinaryKriging(γ), data)

nconstraints(::OrdinaryKriging) = 1

function set_constraints_lhs!(::OrdinaryKriging, LHS::AbstractMatrix, domain)
  T = eltype(LHS)
  LHS[end,:]   .= one(T)
  LHS[:,end]   .= one(T)
  LHS[end,end]  = zero(T)
  nothing
end

factorize(::OrdinaryKriging, LHS::AbstractMatrix) = bunchkaufman(Symmetric(LHS), check=false)

function set_constraints_rhs!(fitted::FittedKriging{<:OrdinaryKriging}, pₒ)
  RHS = fitted.state.RHS
  RHS[end] = one(eltype(RHS))
  nothing
end
