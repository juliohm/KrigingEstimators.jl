# ------------------------------------------------------------------
# Licensed under the ISC License. See LICENCE in the project root.
# ------------------------------------------------------------------

"""
    UniversalKriging(X, z, γ, degree)

## Parameters

* X ∈ ℜ^(mxn) - matrix of data locations
* z ∈ ℜⁿ      - vector of observations for X
* γ           - variogram model
* degree      - polynomial degree for the mean

### Notes

* [`OrdinaryKriging`](@ref) is recovered for 0th degree polynomial
* For non-polynomial mean, see [`ExternalDriftKriging`](@ref)
"""
mutable struct UniversalKriging{T<:Real,V} <: KrigingEstimator
  # input fields
  γ::Variogram
  degree::Int

  # state fields
  X::Matrix{T}
  z::Vector{V}
  LHS::Factorization
  RHS::Vector
  exponents::Matrix{Int}

  function UniversalKriging{T,V}(γ, degree; X=nothing, z=nothing) where {T<:Real,V}
    @assert degree ≥ 0 "degree must be nonnegative"
    UK = new(γ, degree)
    if X ≠ nothing && z ≠ nothing
      fit!(UK, X, z)
    end

    UK
  end
end

UniversalKriging(X, z, γ, degree) = UniversalKriging{eltype(X),eltype(z)}(γ, degree, X=X, z=z)

function nconstraints(estimator::UniversalKriging)
  X = estimator.X
  degree = estimator.degree
  dim, nobs = size(X)

  # multinomial expansion
  expmats = [hcat(collect(multiexponents(dim, d))...) for d in 0:degree]
  exponents = hcat(expmats...)

  # sort expansion for better conditioned Kriging matrices
  sorted = sortperm(vec(maximum(exponents, dims=1)), rev=true)
  exponents = exponents[:,sorted]

  # update object field
  estimator.exponents = exponents

  # return number of terms
  size(exponents, 2)
end

function set_constraints_lhs!(estimator::UniversalKriging, LHS::AbstractMatrix)
  X = estimator.X
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

function set_constraints_rhs!(estimator::UniversalKriging, xₒ::AbstractVector)
  exponents = estimator.exponents
  nterms = size(exponents, 2)
  nobs = length(estimator.z)

  RHS = estimator.RHS
  for j in 1:nterms
    RHS[nobs+j] = prod(xₒ.^exponents[:,j])
  end

  nothing
end

factorize(estimator::UniversalKriging, LHS::AbstractMatrix) = lu(LHS, check=false)
