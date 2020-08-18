# ------------------------------------------------------------------
# Licensed under the ISC License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    ExternalDriftKriging(γ, drifts)
    ExternalDriftKriging(X, z, γ, drifts)

External Drift Kriging with variogram model `γ` and
external `drifts` functions.

Optionally, pass the coordinates `X` and values `z`
to the [`fit`](@ref) function.

### Notes

* External drift functions should be smooth
* Kriging system with external drift is often unstable
* Include a constant drift (e.g. `x->1`) for unbiased estimation
* [`OrdinaryKriging`](@ref) is recovered for `drifts = [x->1]`
* For polynomial mean, see [`UniversalKriging`](@ref)
"""
struct ExternalDriftKriging{G<:Variogram} <: KrigingEstimator
  γ::G
  drifts::Vector{Function}
end

ExternalDriftKriging(γ, drifts) = ExternalDriftKriging{typeof(γ)}(γ, drifts)

ExternalDriftKriging(X, z, γ, drifts) = GeoStatsBase.fit(ExternalDriftKriging(γ, drifts), X ,z)

nconstraints(estimator::ExternalDriftKriging) = length(estimator.drifts)

function set_constraints_lhs!(estimator::ExternalDriftKriging, LHS::AbstractMatrix, X::AbstractMatrix)
  drifts = estimator.drifts
  ndrifts = length(drifts)
  nobs = size(X, 2)
  T = eltype(LHS)

  # set drift blocks
  for i=1:nobs, j=1:ndrifts
    LHS[nobs+j,i] = drifts[j](X[:,i])
    LHS[i,nobs+j] = LHS[nobs+j,i]
  end

  # set zero block
  LHS[nobs+1:end,nobs+1:end] .= zero(T)

  nothing
end

factorize(estimator::ExternalDriftKriging, LHS::AbstractMatrix) = lu(Symmetric(LHS), check=false)

function set_constraints_rhs!(estimator::FittedKriging{E,S},
                              xₒ::AbstractVector) where {E<:ExternalDriftKriging,S<:KrigingState}
  drifts = estimator.estimator.drifts
  RHS = estimator.state.RHS
  nobs = size(estimator.state.X, 2)

  for (j, m) in enumerate(drifts)
    RHS[nobs+j] = m(xₒ)
  end

  nothing
end
