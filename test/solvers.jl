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
        @plottest plot(solution) joinpath(datadir,sname*"1D.png") !isCI
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
      @plottest contourf(solution) joinpath(datadir,"GlobalKriging2D.png") !isCI
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
      @plottest contourf(solution) joinpath(datadir,"NearestKriging2D.png") !isCI
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
      @plottest contourf(solution) joinpath(datadir,"LocalKriging2D.png") !isCI
    end
  end
end
