@testset "Estimators" begin
  # floating point tolerance
  tol = 100eps()

  # create some data
  dim = 3; nobs = 10
  pset = PointSet(10*rand(dim, nobs))
  data = georef((z=rand(nobs),), pset)

  γ = GaussianVariogram(sill=1., range=1., nugget=0.)
  simkrig = SimpleKriging(data, :z, γ, mean(data[:z]))
  ordkrig = OrdinaryKriging(data, :z, γ)
  unikrig = UniversalKriging(data, :z, γ, 1)
  driftkrig = ExternalDriftKriging(data, :z, γ, [x->1.])

  # Kriging is an interpolator
  for j in 1:nobs
    SKestimate, SKvar = predict(simkrig, pset[j])
    OKestimate, OKvar = predict(ordkrig, pset[j])
    UKestimate, UKvar = predict(unikrig, pset[j])
    DKestimate, DKvar = predict(driftkrig, pset[j])

    # estimate checks
    @test SKestimate ≈ data[:z][j]
    @test OKestimate ≈ data[:z][j]
    @test UKestimate ≈ data[:z][j]
    @test DKestimate ≈ data[:z][j]

    # variance checks
    @test SKvar + tol ≥ 0
    @test OKvar + tol ≥ 0
    @test UKvar + tol ≥ 0
    @test DKvar + tol ≥ 0
    @test SKvar ≤ OKvar + tol
  end

  # save results on a particular location pₒ
  pₒ = rand(Point3)
  SKestimate, SKvar = predict(simkrig, pₒ)
  OKestimate, OKvar = predict(ordkrig, pₒ)
  UKestimate, UKvar = predict(unikrig, pₒ)
  DKestimate, DKvar = predict(driftkrig, pₒ)

  # Kriging is translation-invariant
  h = rand(Vec3)
  pset_h = PointSet([pset[i] + h for i in 1:nelements(pset)])
  data_h = georef((z=data[:z],), pset_h)
  simkrig_h = SimpleKriging(data_h, :z, γ, mean(data_h[:z]))
  ordkrig_h = OrdinaryKriging(data_h, :z, γ)
  unikrig_h = UniversalKriging(data_h, :z, γ, 1)
  driftkrig_h = ExternalDriftKriging(data_h, :z, γ, [x->1.])
  SKestimate_h, SKvar_h = predict(simkrig_h, pₒ + h)
  OKestimate_h, OKvar_h = predict(ordkrig_h, pₒ + h)
  UKestimate_h, UKvar_h = predict(unikrig_h, pₒ + h)
  DKestimate_h, DKvar_h = predict(driftkrig_h, pₒ + h)
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
  α = rand()
  γ_α = GaussianVariogram(sill=α, range=1., nugget=0.)
  simkrig_α = SimpleKriging(data, :z, γ_α, mean(data[:z]))
  ordkrig_α = OrdinaryKriging(data, :z, γ_α)
  unikrig_α = UniversalKriging(data, :z, γ_α, 1)
  driftkrig_α = ExternalDriftKriging(data, :z, γ_α, [x->1.])
  SKestimate_α, SKvar_α = predict(simkrig_α, pₒ)
  OKestimate_α, OKvar_α = predict(ordkrig_α, pₒ)
  UKestimate_α, UKvar_α = predict(unikrig_α, pₒ)
  DKestimate_α, DKvar_α = predict(driftkrig_α, pₒ)
  @test SKestimate_α ≈ SKestimate
  @test SKvar_α ≈ α*SKvar
  @test OKestimate_α ≈ OKestimate
  @test OKvar_α ≈ α*OKvar
  @test UKestimate_α ≈ UKestimate
  @test UKvar_α ≈ α*UKvar
  @test DKestimate_α ≈ DKestimate
  @test DKvar_α ≈ α*DKvar

  # Kriging variance is a function of data configuration, not data values
  δ = rand(nobs)
  data_δ = georef((z=data[:z].+δ,), pset)
  simkrig_δ = SimpleKriging(data_δ, :z, γ, mean(data_δ[:z]))
  ordkrig_δ = OrdinaryKriging(data_δ, :z, γ)
  unikrig_δ = UniversalKriging(data_δ, :z, γ, 1)
  driftkrig_δ = ExternalDriftKriging(data_δ, :z, γ, [x->1.])
  SKestimate_δ, SKvar_δ = predict(simkrig_δ, pₒ)
  OKestimate_δ, OKvar_δ = predict(ordkrig_δ, pₒ)
  UKestimate_δ, UKvar_δ = predict(unikrig_δ, pₒ)
  DKestimate_δ, DKvar_δ = predict(driftkrig_δ, pₒ)
  @test SKvar_δ ≈ SKvar
  @test OKvar_δ ≈ OKvar
  @test UKvar_δ ≈ UKvar
  @test DKvar_δ ≈ DKvar

  # Ordinary Kriging ≡ Universal Kriging with 0th degree drift
  unikrig_0th = UniversalKriging(data, :z, γ, 0)
  OKestimate, OKvar = predict(ordkrig, pₒ)
  UKestimate, UKvar = predict(unikrig_0th, pₒ)
  @test OKestimate ≈ UKestimate
  @test OKvar ≈ UKvar

  # Ordinary Kriging ≡ Kriging with constant external drift
  driftkrig_const = ExternalDriftKriging(data, :z, γ, [x->1.])
  OKestimate, OKvar = predict(ordkrig, pₒ)
  DKestimate, DKvar = predict(driftkrig_const, pₒ)
  @test OKestimate ≈ DKestimate
  @test OKvar ≈ DKvar

  # Non-stationary variograms are allowed
  γ_ns = PowerVariogram()
  ordkrig_ns = OrdinaryKriging(data, :z, γ_ns)
  unikrig_ns = UniversalKriging(data, :z, γ_ns, 1)
  driftkrig_ns = ExternalDriftKriging(data, :z, γ_ns, [x->1.])
  for j in 1:nobs
    OKestimate, OKvar = predict(ordkrig_ns, pset[j])
    UKestimate, UKvar = predict(unikrig_ns, pset[j])
    DKestimate, DKvar = predict(driftkrig_ns, pset[j])

    # estimate checks
    @test OKestimate ≈ data[:z][j]
    @test UKestimate ≈ data[:z][j]
    @test DKestimate ≈ data[:z][j]

    # variance checks
    @test OKvar + tol ≥ 0
    @test UKvar + tol ≥ 0
    @test DKvar + tol ≥ 0
  end

  # Floating point precision checks
  X_f    = rand(Float32, dim, nobs)
  z_f    = rand(Float32, nobs)
  X_d    = Float64.(X_f)
  z_d    = Float64.(z_f)
  pset_f = PointSet(X_f)
  data_f = georef((z=z_f,), pset_f)
  pset_d = PointSet(X_d)
  data_d = georef((z=z_d,), pset_d)
  pₒ_f   = rand(Point{dim,Float32})
  pₒ_d   = convert(Point{dim,Float64}, pₒ_f)
  γ_f         = GaussianVariogram(sill=1f0, range=1f0, nugget=0f0)
  simkrig_f   = SimpleKriging(data_f, :z, γ_f, mean(data_f[:z]))
  ordkrig_f   = OrdinaryKriging(data_f, :z, γ_f)
  unikrig_f   = UniversalKriging(data_f, :z, γ_f, 1)
  driftkrig_f = ExternalDriftKriging(data_f, :z, γ_f, [x->1f0])
  γ_d         = GaussianVariogram(sill=1., range=1., nugget=0.)
  simkrig_d   = SimpleKriging(data_d, :z, γ_d, mean(data_d[:z]))
  ordkrig_d   = OrdinaryKriging(data_d, :z, γ_d)
  unikrig_d   = UniversalKriging(data_d, :z, γ_d, 1)
  driftkrig_d = ExternalDriftKriging(data_d, :z, γ_d, [x->1.])
  SKestimate_f, SKvar_f = predict(simkrig_f, pₒ_f)
  OKestimate_f, OKvar_f = predict(ordkrig_f, pₒ_f)
  UKestimate_f, UKvar_f = predict(unikrig_f, pₒ_f)
  DKestimate_f, DKvar_f = predict(driftkrig_f, pₒ_f)
  SKestimate_d, SKvar_d = predict(simkrig_d, pₒ_d)
  OKestimate_d, OKvar_d = predict(ordkrig_d, pₒ_d)
  UKestimate_d, UKvar_d = predict(unikrig_d, pₒ_d)
  DKestimate_d, DKvar_d = predict(driftkrig_d, pₒ_d)
  @test isapprox(SKestimate_f, SKestimate_d, atol=1e-4)
  @test isapprox(SKvar_f, SKvar_d, atol=1e-4)
  @test isapprox(OKestimate_f, OKestimate_d, atol=1e-4)
  @test isapprox(OKvar_f, OKvar_d, atol=1e-4)
  @test isapprox(UKestimate_f, UKestimate_d, atol=1e-4)
  @test isapprox(UKvar_f, UKvar_d, atol=1e-4)
  @test isapprox(DKestimate_f, DKestimate_d, atol=1e-4)
  @test isapprox(DKvar_f, DKvar_d, atol=1e-4)
end
