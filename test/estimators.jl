@testset "Estimators" begin
  # floating point tolerance
  tol = 100eps()

  # create some data
  dim = 3
  nobs = 10
  pset = PointSet(10 * rand(dim, nobs))
  data = georef((z=rand(nobs),), pset)

  γ = GaussianVariogram(sill=1.0, range=1.0, nugget=0.0)
  simkrig = SimpleKriging(data, γ, mean(data.z))
  ordkrig = OrdinaryKriging(data, γ)
  unikrig = UniversalKriging(data, γ, 1)
  drikrig = ExternalDriftKriging(data, γ, [x -> 1.0])

  # Kriging is an interpolator
  for j in 1:nobs
    SKestimate, SKvar = predict(simkrig, :z, pset[j])
    OKestimate, OKvar = predict(ordkrig, :z, pset[j])
    UKestimate, UKvar = predict(unikrig, :z, pset[j])
    DKestimate, DKvar = predict(drikrig, :z, pset[j])

    # estimate checks
    @test SKestimate ≈ data.z[j]
    @test OKestimate ≈ data.z[j]
    @test UKestimate ≈ data.z[j]
    @test DKestimate ≈ data.z[j]

    # variance checks
    @test SKvar + tol ≥ 0
    @test OKvar + tol ≥ 0
    @test UKvar + tol ≥ 0
    @test DKvar + tol ≥ 0
    @test SKvar ≤ OKvar + tol
  end

  # save results on a particular location pₒ
  pₒ = rand(Point3)
  SKestimate, SKvar = predict(simkrig, :z, pₒ)
  OKestimate, OKvar = predict(ordkrig, :z, pₒ)
  UKestimate, UKvar = predict(unikrig, :z, pₒ)
  DKestimate, DKvar = predict(drikrig, :z, pₒ)

  # Kriging is translation-invariant
  h = rand(Vec3)
  pset_h = PointSet([pset[i] + h for i in 1:nelements(pset)])
  data_h = georef((z=data.z,), pset_h)
  simkrig_h = SimpleKriging(data_h, γ, mean(data_h.z))
  ordkrig_h = OrdinaryKriging(data_h, γ)
  unikrig_h = UniversalKriging(data_h, γ, 1)
  drikrig_h = ExternalDriftKriging(data_h, γ, [x -> 1.0])
  SKestimate_h, SKvar_h = predict(simkrig_h, :z, pₒ + h)
  OKestimate_h, OKvar_h = predict(ordkrig_h, :z, pₒ + h)
  UKestimate_h, UKvar_h = predict(unikrig_h, :z, pₒ + h)
  DKestimate_h, DKvar_h = predict(drikrig_h, :z, pₒ + h)
  @test SKestimate_h ≈ SKestimate
  @test SKvar_h ≈ SKvar
  @test OKestimate_h ≈ OKestimate
  @test OKvar_h ≈ OKvar
  @test UKestimate_h ≈ UKestimate
  @test UKvar_h ≈ UKvar
  @test DKestimate_h ≈ DKestimate
  @test DKvar_h ≈ DKvar

  # Kriging estimate is invariant under covariance scaling
  # Kriging variance is multiplied by the same factor
  α = 2.0
  γ_α = GaussianVariogram(sill=α, range=1.0, nugget=0.0)
  simkrig_α = SimpleKriging(data, γ_α, mean(data.z))
  ordkrig_α = OrdinaryKriging(data, γ_α)
  unikrig_α = UniversalKriging(data, γ_α, 1)
  drikrig_α = ExternalDriftKriging(data, γ_α, [x -> 1.0])
  SKestimate_α, SKvar_α = predict(simkrig_α, :z, pₒ)
  OKestimate_α, OKvar_α = predict(ordkrig_α, :z, pₒ)
  UKestimate_α, UKvar_α = predict(unikrig_α, :z, pₒ)
  DKestimate_α, DKvar_α = predict(drikrig_α, :z, pₒ)
  @test SKestimate_α ≈ SKestimate
  @test SKvar_α ≈ α * SKvar
  @test OKestimate_α ≈ OKestimate
  @test OKvar_α ≈ α * OKvar
  @test UKestimate_α ≈ UKestimate
  @test UKvar_α ≈ α * UKvar
  @test DKestimate_α ≈ DKestimate
  @test DKvar_α ≈ α * DKvar

  # Kriging variance is a function of data configuration, not data values
  δ = rand(nobs)
  data_δ = georef((z=data.z .+ δ,), pset)
  simkrig_δ = SimpleKriging(data_δ, γ, mean(data_δ.z))
  ordkrig_δ = OrdinaryKriging(data_δ, γ)
  unikrig_δ = UniversalKriging(data_δ, γ, 1)
  drikrig_δ = ExternalDriftKriging(data_δ, γ, [x -> 1.0])
  SKestimate_δ, SKvar_δ = predict(simkrig_δ, :z, pₒ)
  OKestimate_δ, OKvar_δ = predict(ordkrig_δ, :z, pₒ)
  UKestimate_δ, UKvar_δ = predict(unikrig_δ, :z, pₒ)
  DKestimate_δ, DKvar_δ = predict(drikrig_δ, :z, pₒ)
  @test SKvar_δ ≈ SKvar
  @test OKvar_δ ≈ OKvar
  @test UKvar_δ ≈ UKvar
  @test DKvar_δ ≈ DKvar

  # Ordinary Kriging ≡ Universal Kriging with 0th degree drift
  unikrig_0th = UniversalKriging(data, γ, 0)
  OKestimate, OKvar = predict(ordkrig, :z, pₒ)
  UKestimate, UKvar = predict(unikrig_0th, :z, pₒ)
  @test OKestimate ≈ UKestimate
  @test OKvar ≈ UKvar

  # Ordinary Kriging ≡ Kriging with constant external drift
  driftkrig_const = ExternalDriftKriging(data, γ, [x -> 1.0])
  OKestimate, OKvar = predict(ordkrig, :z, pₒ)
  DKestimate, DKvar = predict(driftkrig_const, :z, pₒ)
  @test OKestimate ≈ DKestimate
  @test OKvar ≈ DKvar

  # Non-stationary variograms are allowed
  γ_ns = PowerVariogram()
  ordkrig_ns = OrdinaryKriging(data, γ_ns)
  unikrig_ns = UniversalKriging(data, γ_ns, 1)
  drikrig_ns = ExternalDriftKriging(data, γ_ns, [x -> 1.0])
  for j in 1:nobs
    OKestimate, OKvar = predict(ordkrig_ns, :z, pset[j])
    UKestimate, UKvar = predict(unikrig_ns, :z, pset[j])
    DKestimate, DKvar = predict(drikrig_ns, :z, pset[j])

    # estimate checks
    @test OKestimate ≈ data.z[j]
    @test UKestimate ≈ data.z[j]
    @test DKestimate ≈ data.z[j]

    # variance checks
    @test OKvar + tol ≥ 0
    @test UKvar + tol ≥ 0
    @test DKvar + tol ≥ 0
  end

  # Floating point precision checks
  X_f = rand(Float32, dim, nobs)
  z_f = rand(Float32, nobs)
  X_d = Float64.(X_f)
  z_d = Float64.(z_f)
  pset_f = PointSet(X_f)
  data_f = georef((z=z_f,), pset_f)
  pset_d = PointSet(X_d)
  data_d = georef((z=z_d,), pset_d)
  pₒ_f = rand(Point{dim,Float32})
  pₒ_d = convert(Point{dim,Float64}, pₒ_f)
  γ_f = GaussianVariogram(sill=1.0f0, range=1.0f0, nugget=0.0f0)
  simkrig_f = SimpleKriging(data_f, γ_f, mean(data_f.z))
  ordkrig_f = OrdinaryKriging(data_f, γ_f)
  unikrig_f = UniversalKriging(data_f, γ_f, 1)
  drikrig_f = ExternalDriftKriging(data_f, γ_f, [x -> 1.0f0])
  γ_d = GaussianVariogram(sill=1.0, range=1.0, nugget=0.0)
  simkrig_d = SimpleKriging(data_d, γ_d, mean(data_d.z))
  ordkrig_d = OrdinaryKriging(data_d, γ_d)
  unikrig_d = UniversalKriging(data_d, γ_d, 1)
  drikrig_d = ExternalDriftKriging(data_d, γ_d, [x -> 1.0])
  SKestimate_f, SKvar_f = predict(simkrig_f, :z, pₒ_f)
  OKestimate_f, OKvar_f = predict(ordkrig_f, :z, pₒ_f)
  UKestimate_f, UKvar_f = predict(unikrig_f, :z, pₒ_f)
  DKestimate_f, DKvar_f = predict(drikrig_f, :z, pₒ_f)
  SKestimate_d, SKvar_d = predict(simkrig_d, :z, pₒ_d)
  OKestimate_d, OKvar_d = predict(ordkrig_d, :z, pₒ_d)
  UKestimate_d, UKvar_d = predict(unikrig_d, :z, pₒ_d)
  DKestimate_d, DKvar_d = predict(drikrig_d, :z, pₒ_d)
  @test isapprox(SKestimate_f, SKestimate_d, atol=1e-4)
  @test isapprox(SKvar_f, SKvar_d, atol=1e-4)
  @test isapprox(OKestimate_f, OKestimate_d, atol=1e-4)
  @test isapprox(OKvar_f, OKvar_d, atol=1e-4)
  @test isapprox(UKestimate_f, UKestimate_d, atol=1e-4)
  @test isapprox(UKvar_f, UKvar_d, atol=1e-4)
  @test isapprox(DKestimate_f, DKestimate_d, atol=1e-4)
  @test isapprox(DKvar_f, DKvar_d, atol=1e-4)

  # ------------------
  # CHANGE OF SUPPORT
  # ------------------

  # create some data
  dim = 2
  nobs = 10
  pset = PointSet(10 * rand(dim, nobs))
  data = georef((z=rand(nobs),), pset)

  # basic estimators
  γ = GaussianVariogram(sill=1.0, range=1.0, nugget=0.0)
  simkrig = SimpleKriging(data, γ, mean(data.z))
  ordkrig = OrdinaryKriging(data, γ)
  unikrig = UniversalKriging(data, γ, 1)
  drikrig = ExternalDriftKriging(data, γ, [x -> 1.0])

  # prediction on a quadrangle
  uₒ = Quadrangle((0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0))
  _, SKvar = predict(simkrig, :z, uₒ)
  _, OKvar = predict(ordkrig, :z, uₒ)
  _, UKvar = predict(unikrig, :z, uₒ)
  _, DKvar = predict(drikrig, :z, uₒ)

  # variance checks
  @test SKvar + tol ≥ 0
  @test OKvar + tol ≥ 0
  @test UKvar + tol ≥ 0
  @test DKvar + tol ≥ 0
  @test SKvar ≤ OKvar + tol

  # ------------------
  # COMPOSITIONAL DATA
  # ------------------

  # create some data
  dim = 2
  nobs = 10
  pset = PointSet(10 * rand(dim, nobs))
  table = (z=rand(Composition{3}, nobs),)
  data = georef(table, pset)

  # basic estimators
  γ = GaussianVariogram(sill=1.0, range=1.0, nugget=0.0)
  simkrig = SimpleKriging(data, γ, mean(data.z))
  ordkrig = OrdinaryKriging(data, γ)
  unikrig = UniversalKriging(data, γ, 1)
  drikrig = ExternalDriftKriging(data, γ, [x -> 1.0])

  # prediction on a quadrangle
  uₒ = Quadrangle((0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0))
  SKpred, SKvar = predict(simkrig, :z, uₒ)
  OKpred, OKvar = predict(ordkrig, :z, uₒ)
  UKpred, UKvar = predict(unikrig, :z, uₒ)
  DKpred, DKvar = predict(drikrig, :z, uₒ)

  # type tests
  @test SKpred isa Composition
  @test OKpred isa Composition
  @test UKpred isa Composition
  @test DKpred isa Composition

  # variance checks
  @test SKvar + tol ≥ 0
  @test OKvar + tol ≥ 0
  @test UKvar + tol ≥ 0
  @test DKvar + tol ≥ 0
  @test SKvar ≤ OKvar + tol

  # create some unitful data
  dim = 3
  nobs = 10
  pset = PointSet(10 * rand(dim, nobs))
  data = georef((z=rand(nobs) * u"K",), pset)

  γ = GaussianVariogram(sill=1.0u"K^2")
  sk = SimpleKriging(data, γ, mean(data.z))
  ok = OrdinaryKriging(data, γ)
  uk = UniversalKriging(data, γ, 1)
  dk = ExternalDriftKriging(data, γ, [x -> 1.0])
  for _k in [sk, ok, uk, dk]
    μ, σ² = predict(dk, :z, Point(0.0, 0.0, 0.0))
    @test unit(μ) == u"K"
    @test unit(σ²) == u"K^2"
  end

  # probabilistic predictions
  dim = 3
  nobs = 10
  pset = PointSet(10 * rand(dim, nobs))
  data = georef((z=rand(nobs),), pset)

  γ = GaussianVariogram(sill=1.0, range=1.0, nugget=0.0)
  sk = SimpleKriging(data, γ, mean(data.z))
  ok = OrdinaryKriging(data, γ)
  uk = UniversalKriging(data, γ, 1)
  dk = ExternalDriftKriging(data, γ, [x -> 1.0])

  skpred = predict(sk, :z, Point3(5.0, 5.0, 5.0))
  okpred = predict(ok, :z, Point3(5.0, 5.0, 5.0))
  ukpred = predict(uk, :z, Point3(5.0, 5.0, 5.0))
  dkpred = predict(dk, :z, Point3(5.0, 5.0, 5.0))
  skprob = predictprob(sk, :z, Point3(5.0, 5.0, 5.0))
  okprob = predictprob(ok, :z, Point3(5.0, 5.0, 5.0))
  ukprob = predictprob(uk, :z, Point3(5.0, 5.0, 5.0))
  dkprob = predictprob(dk, :z, Point3(5.0, 5.0, 5.0))
  @test mean(skprob) ≈ first(skpred)
  @test mean(okprob) ≈ first(okpred)
  @test mean(ukprob) ≈ first(ukpred)
  @test mean(dkprob) ≈ first(dkpred)
  @test var(skprob) ≈ last(skpred)
  @test var(okprob) ≈ last(okpred)
  @test var(ukprob) ≈ last(ukpred)
  @test var(dkprob) ≈ last(dkpred)
end
