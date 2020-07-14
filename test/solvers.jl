@testset "Solvers" begin
  data1D = readgeotable(joinpath(datadir,"data1D.tsv"), coordnames=(:x,))
  data2D = readgeotable(joinpath(datadir,"data2D.tsv"), coordnames=(:x,:y))
  grid1D = RegularGrid(100)
  grid2D = RegularGrid(100,100)

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
    inds = LinearIndices(size(grid2D))
    for solution in solutions2D
      M, V = solution[:value]
      @test isapprox(M[inds[26,26]], 1., atol=1e-6)
      @test isapprox(M[inds[51,76]], 0., atol=1e-6)
      @test isapprox(M[inds[76,51]], 1., atol=1e-6)
    end

    if visualtests
      gr(size=(800,400))
      for i in 1:2
        solution, sname = solutions1D[i], snames[i]
        @plottest plot(solution) joinpath(datadir,sname*"1D.png") !istravis
      end
      for i in 1:2
        solution, sname = solutions2D[i], snames[i]
        @plottest contourf(solution) joinpath(datadir,sname*"2D.png") !istravis
      end
      # TODO: test local_kriging
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
