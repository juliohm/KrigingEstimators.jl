@testset "Solvers" begin
  data1D = readgeotable(joinpath(datadir,"data1D.tsv"), coordnames=(:x,))
  data2D = readgeotable(joinpath(datadir,"data2D.tsv"), coordnames=(:x,:y))
  grid1D = RegularGrid(100)
  grid2D = RegularGrid(100,100)

  @testset "Kriging1D" begin
    problem = EstimationProblem(data1D, grid1D, :value)

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

    solvers   = [global_kriging, nearest_kriging, local_kriging]
    solnames  = ["GlobalKriging", "NearestKriging", "LocalKriging"]
    solutions = [solve(problem, solver) for solver in solvers]

    if visualtests
      gr(size=(800,400))
      for i in 1:3
        solution, sname = solutions[i], solnames[i]
        @plottest plot(solution) joinpath(datadir,sname*"1D.png") !istravis
      end
    end
  end

  @testset "GlobalKriging2D" begin
    problem = EstimationProblem(data2D, grid2D, :value)

    solver = Kriging(
      :value => (variogram=GaussianVariogram(range=35.,nugget=0.),)
    )

    solution = solve(problem, solver)

    # basic checks
    inds = LinearIndices(size(grid2D))
    M, V = solution[:value]
    @test isapprox(M[inds[26,26]], 1., atol=1e-6)
    @test isapprox(M[inds[51,76]], 0., atol=1e-6)
    @test isapprox(M[inds[76,51]], 1., atol=1e-6)

    if visualtests
      gr(size=(800,400))
      @plottest contourf(solution) joinpath(datadir,"GlobalKriging2D.png") !istravis
    end
  end

  @testset "NearestKriging2D" begin
    problem = EstimationProblem(data2D, grid2D, :value)

    solver = Kriging(
      :value => (variogram=GaussianVariogram(range=35.,nugget=0.), maxneighbors=3)
    )

    solution = solve(problem, solver)

    # basic checks
    inds = LinearIndices(size(grid2D))
    M, V = solution[:value]
    @test isapprox(M[inds[26,26]], 1., atol=1e-6)
    @test isapprox(M[inds[51,76]], 0., atol=1e-6)
    @test isapprox(M[inds[76,51]], 1., atol=1e-6)

    if visualtests
      gr(size=(800,400))
      @plottest contourf(solution) joinpath(datadir,"NearestKriging2D.png") !istravis
    end
  end

  @testset "LocalKriging2D" begin
    problem = EstimationProblem(data2D, grid2D, :value)

    solver = Kriging(
      :value => (variogram=GaussianVariogram(range=35.,nugget=0.),
                 maxneighbors=3, neighborhood=BallNeighborhood(100.))
    )

    solution = solve(problem, solver)

    # basic checks
    inds = LinearIndices(size(grid2D))
    M, V = solution[:value]
    @test isapprox(M[inds[26,26]], 1., atol=1e-6)
    @test isapprox(M[inds[51,76]], 0., atol=1e-6)
    @test isapprox(M[inds[76,51]], 1., atol=1e-6)

    if visualtests
      gr(size=(800,400))
      @plottest contourf(solution) joinpath(datadir,"LocalKriging2D.png") !istravis
    end
  end

  @testset "SeqGaussSim" begin
    nreals = 3

    @testset "Conditional" begin
      problem = SimulationProblem(data2D, grid2D, :value, nreals)

      solver = SeqGaussSim(
        :value => (variogram=GaussianVariogram(range=35.),
                   neighborhood=BallNeighborhood(10.))
      )

      Random.seed!(2017)
      solution = solve(problem, solver)

      # basic checks
      reals = solution[:value]
      inds = LinearIndices(size(grid2D))
      @test all(reals[i][inds[26,26]] == 1. for i in 1:nreals)
      @test all(reals[i][inds[51,76]] == 0. for i in 1:nreals)
      @test all(reals[i][inds[76,51]] == 1. for i in 1:nreals)

      if visualtests
        gr(size=(800,400))
        # @plottest plot(solution) joinpath(datadir,"SGSCond2D.png") !istravis
      end
    end

    @testset "Unconditional" begin
      problem = SimulationProblem(grid2D, :value => Float64, nreals)

      solver = SeqGaussSim(
        :value => (variogram=GaussianVariogram(range=35.),
                   neighborhood=BallNeighborhood(10.))
      )

      Random.seed!(2017)
      solution = solve(problem, solver)

      if visualtests
        gr(size=(800,400))
        # @plottest plot(solution) joinpath(datadir,"SGSUncond2D.png") !istravis
      end
    end
  end
end