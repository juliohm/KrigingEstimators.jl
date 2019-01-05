# ------------------------------------------------------------------
# Licensed under the ISC License. See LICENCE in the project root.
# ------------------------------------------------------------------

"""
    SimpleKriging(γ, μ)
    SimpleKriging(X, z, γ, μ)

Simple Kriging with variogram model `γ` and constant mean `μ`.

Optionally, pass the coordinates `X` and values `z`
to the [`fit`](@ref) function.

### Notes

* Simple Kriging requires stationary variograms
"""
struct SimpleKriging{G<:Variogram,V} <: KrigingEstimator
  # input fields
  γ::G
  μ::V

  function SimpleKriging{G,V}(γ, μ) where {G<:Variogram,V}
    @assert isstationary(γ) "Simple Kriging requires stationary variogram"
    new(γ, μ)
  end
end

SimpleKriging(γ, μ) = SimpleKriging{typeof(γ),typeof(μ)}(γ, μ)

SimpleKriging(X, z, γ, μ) = fit(SimpleKriging(γ, μ), X, z)

nconstraints(estimator::SimpleKriging) = 0

set_constraints_lhs!(estimator::SimpleKriging, LHS::AbstractMatrix, X::AbstractMatrix) = nothing

factorize(estimator::SimpleKriging, LHS::AbstractMatrix) = cholesky(LHS, check=false)

set_constraints_rhs!(estimator::FittedKriging{E,S},
                     xₒ::AbstractVector) where {E<:SimpleKriging,S<:KrigingState} = nothing

function combine(estimator::FittedKriging{E,S},
                 weights::KrigingWeights, z::AbstractVector) where {E<:SimpleKriging,S<:KrigingState}
  γ = estimator.estimator.γ
  μ = estimator.estimator.μ
  b = estimator.state.RHS
  λ = weights.λ
  y = z .- μ

  μ + y⋅λ, sill(γ) - b⋅λ
end
