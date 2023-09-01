# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    UniversalKriging(γ, degree, dim)
    UniversalKriging(data, γ, degree)

Universal Kriging with variogram model `γ` and polynomial
`degree` on a geospatial domain of dimension `dim`.

Optionally, pass the geospatial `data` to the [`fit`](@ref) function.

### Notes

* [`OrdinaryKriging`](@ref) is recovered for 0th degree polynomial
* For non-polynomial mean, see [`ExternalDriftKriging`](@ref)
"""
struct UniversalKriging{G<:Variogram} <: KrigingEstimator
  γ::G
  degree::Int
  dim::Int
  exponents::Matrix{Int}

  function UniversalKriging{G}(γ, degree, dim) where {G<:Variogram}
    @assert degree ≥ 0 "degree must be nonnegative"
    @assert dim > 0 "dimension must be positive"
    exponents = UKexps(degree, dim)
    new(γ, degree, dim, exponents)
  end
end

UniversalKriging(γ, degree, dim) = UniversalKriging{typeof(γ)}(γ, degree, dim)

UniversalKriging(data::AbstractGeoTable, γ, degree) =
  GeoStatsBase.fit(UniversalKriging(γ, degree, embeddim(domain(data))), data)

function UKexps(degree::Int, dim::Int)
  # multinomial expansion
  expmats = [hcat(collect(multiexponents(dim, d))...) for d in 0:degree]
  exponents = hcat(expmats...)

  # sort expansion for better conditioned Kriging matrices
  sorted = sortperm(vec(maximum(exponents, dims=1)), rev=true)

  exponents[:, sorted]
end

nconstraints(estimator::UniversalKriging) = size(estimator.exponents, 2)

function set_constraints_lhs!(estimator::UniversalKriging, LHS::AbstractMatrix, domain)
  exponents = estimator.exponents
  nobs = nelements(domain)
  nterms = size(exponents, 2)

  # set polynomial drift blocks
  for i in 1:nobs
    x = coordinates(centroid(domain, i))
    for j in 1:nterms
      LHS[nobs + j, i] = prod(x .^ exponents[:, j])
      LHS[i, nobs + j] = LHS[nobs + j, i]
    end
  end

  # set zero block
  LHS[(nobs + 1):end, (nobs + 1):end] .= zero(eltype(LHS))

  nothing
end

function set_constraints_rhs!(fitted::FittedKriging{<:UniversalKriging}, uₒ)
  exponents = fitted.estimator.exponents
  RHS = fitted.state.RHS
  nobs = nrow(fitted.state.data)
  nterms = size(exponents, 2)

  # set polynomial drift
  xₒ = coordinates(centroid(uₒ))
  for j in 1:nterms
    RHS[nobs + j] = prod(xₒ .^ exponents[:, j])
  end

  nothing
end
