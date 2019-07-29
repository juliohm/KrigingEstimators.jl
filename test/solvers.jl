@testset "Solvers" begin
  data1D = readgeotable(joinpath(datadir,"data1D.tsv"), delim='\t', coordnames=[:x])
  data2D = readgeotable(joinpath(datadir,"data2D.tsv"), delim='\t', coordnames=[:x,:y])
  grid1D = RegularGrid{Float64}(100)
  grid2D = RegularGrid{Float64}(100,100)

  @testset "Kriging" begin
    problem1D = EstimationProblem(data1D, grid1D, :value)
    problem2D = EstimationProblem(data2D, grid2D, :value)

    global_kriging = Kriging(
      :value => (variogram=GaussianVariogram(range=35.,nugget=0.),)
    )
    nearest_kriging = Kriging(
      :value => (variogram=GaussianVariogram(range=35.,nugget=0.), maxneighbors=3)
    )
    local_kriging = Kriging(
      :value => (variogram=GaussianVariogram(range=35.,nugget=0.),
                 maxneighbors=3, neighborhood=BallNeighborhood(100.))
    )

    solvers = [global_kriging, nearest_kriging, local_kriging]
    snames  = ["GlobalKriging", "NearestKriging", "LocalKriging"]

    solutions1D = [solve(problem1D, solver) for solver in solvers]
    solutions2D = [solve(problem2D, solver) for solver in solvers]

    # basic checks
    for solution in solutions2D
      M, V = solution[:value]
      @test M[26,26] == 1.
      @test M[51,76] == 0.
      @test M[76,51] == 1.
    end

    if visualtests
      gr(size=(800,400))
      for i in 1:2
        solution, sname = solutions1D[i], snames[i]
        @plottest plot(solution) joinpath(datadir,sname*"1D.png") !istravis
      end
      for (solution, sname) in zip(solutions2D, snames)
        @plottest contourf(solution) joinpath(datadir,sname*"2D.png") !istravis
      end
    end
  end

  @testset "SeqGaussSim" begin
    nreals = 3

    ###################
    ##  CONDITIONAL  ##
    ###################
    problem = SimulationProblem(data2D, grid2D, :value, nreals)

    solver = SeqGaussSim(
      :value => (variogram=GaussianVariogram(range=35.),
                 neighborhood=BallNeighborhood(10.))
    )

    Random.seed!(2017)
    solution = solve(problem, solver)

    # basic checks
    reals = solution[:value]
    @test all(reals[i][26,26] == 1. for i in 1:nreals)
    @test all(reals[i][51,76] == 0. for i in 1:nreals)
    @test all(reals[i][76,51] == 1. for i in 1:nreals)

    if visualtests
      gr(size=(800,400))
      @plottest plot(solution) joinpath(datadir,"SGSCond2D.png") !istravis
    end

    ###################
    ## UNCONDITIONAL ##
    ###################
    problem = SimulationProblem(grid2D, :value => Float64, nreals)

    solver = SeqGaussSim(
      :value => (variogram=GaussianVariogram(range=35.),
                 neighborhood=BallNeighborhood(10.))
    )

    Random.seed!(2017)
    solution = solve(problem, solver)

    if visualtests
      gr(size=(800,400))
      @plottest plot(solution) joinpath(datadir,"SGSUncond2D.png") !istravis
    end
  end

  @testset "CookieCutter" begin
    problem = SimulationProblem(grid2D, (:facies => Int, :property => Float64), 3)

    γ₀ = GaussianVariogram(distance=Ellipsoidal([30.,10.],[0.]))
    γ₁ = GaussianVariogram(distance=Ellipsoidal([10.,30.],[0.]))
    solver = CookieCutter(Dummy(:facies => NamedTuple()),
                          [0 => SeqGaussSim(:property => (variogram=γ₀,neighborhood=BallNeighborhood(10.))),
                           1 => SeqGaussSim(:property => (variogram=γ₁,neighborhood=BallNeighborhood(10.)))])

    Random.seed!(1234)
    solution = solve(problem, solver)

    if visualtests
      gr(size=(800,600))
      @plottest plot(solution) joinpath(datadir,"CookieCutter.png") !istravis
    end
  end
end
