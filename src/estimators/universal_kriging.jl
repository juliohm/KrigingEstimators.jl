# ------------------------------------------------------------------
# Licensed under the ISC License. See LICENCE in the project root.
# ------------------------------------------------------------------

"""
    UniversalKriging(γ, degree, dim)
    UniversalKriging(X, z, γ, degree)

Universal Kriging with variogram model `γ` and polynomial
`degree` on a spatial domain of dimension `dim`.

Optionally, pass the coordinates `X` and values `z`
to the [`fit`](@ref) function.

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

UniversalKriging(X, z, γ, degree) = fit(UniversalKriging(γ, degree, size(X,1)), X, z)

function UKexps(degree::Int, dim::Int)
  # multinomial expansion
  expmats = [hcat(collect(multiexponents(dim, d))...) for d in 0:degree]
  exponents = hcat(expmats...)

  # sort expansion for better conditioned Kriging matrices
  sorted = sortperm(vec(maximum(exponents, dims=1)), rev=true)

  exponents[:,sorted]
end

nconstraints(estimator::UniversalKriging) = size(estimator.exponents, 2)

function set_constraints_lhs!(estimator::UniversalKriging, LHS::AbstractMatrix, X::AbstractMatrix)
  exponents = estimator.exponents
  nobs = size(X, 2)
  nterms = size(exponents, 2)
  T = eltype(LHS)

  # set polynomial drift blocks
  for i=1:nobs, j=1:nterms
    LHS[nobs+j,i] = prod(X[:,i].^exponents[:,j])
    LHS[i,nobs+j] = LHS[nobs+j,i]
  end

  # set zero block
  LHS[nobs+1:end,nobs+1:end] .= zero(T)

  nothing
end

factorize(estimator::UniversalKriging, LHS::AbstractMatrix) = lu(LHS, check=false)

function set_constraints_rhs!(estimator::FittedKriging{E,S},
                              xₒ::AbstractVector) where {E<:UniversalKriging,S<:KrigingState}
  exponents = estimator.estimator.exponents
  RHS = estimator.state.RHS
  nobs = size(estimator.state.X, 2)
  nterms = size(exponents, 2)

  for j in 1:nterms
    RHS[nobs+j] = prod(xₒ.^exponents[:,j])
  end

  nothing
end
