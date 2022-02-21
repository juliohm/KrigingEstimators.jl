# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    ExternalDriftKriging(γ, drifts)
    ExternalDriftKriging(data, γ, drifts)

External Drift Kriging with variogram model `γ` and
external `drifts` functions.

Optionally, pass the geospatial `data` to the [`fit`](@ref) function.

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

ExternalDriftKriging(data, γ, drifts) = GeoStatsBase.fit(ExternalDriftKriging(γ, drifts), data)

nconstraints(estimator::ExternalDriftKriging) = length(estimator.drifts)

function set_constraints_lhs!(estimator::ExternalDriftKriging, LHS::AbstractMatrix, domain)
  drifts = estimator.drifts
  ndrifts = length(drifts)
  nobs = nelements(domain)

  # set external drift blocks
  for i in 1:nobs
    x = coordinates(centroid(domain, i))
    for j in 1:ndrifts
      LHS[nobs+j,i] = drifts[j](x)
      LHS[i,nobs+j] = LHS[nobs+j,i]
    end
  end

  # set zero block
  LHS[nobs+1:end,nobs+1:end] .= zero(eltype(LHS))

  nothing
end

factorize(::ExternalDriftKriging, LHS::AbstractMatrix) = bunchkaufman(Symmetric(LHS), check=false)

function set_constraints_rhs!(fitted::FittedKriging{<:ExternalDriftKriging}, uₒ)
  drifts = fitted.estimator.drifts
  RHS = fitted.state.RHS
  nobs = nelements(fitted.state.data)

  # set external drift
  xₒ = coordinates(centroid(uₒ))
  for (j, m) in enumerate(drifts)
    RHS[nobs+j] = m(xₒ)
  end

  nothing
end
