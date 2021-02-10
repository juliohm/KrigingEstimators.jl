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
        @test_ref_plot "data/$(sname)1D.png" plot(solution)
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
    S = solution[:value]
    @test isapprox(S[inds[26,26]], 1., atol=1e-6)
    @test isapprox(S[inds[51,76]], 0., atol=1e-6)
    @test isapprox(S[inds[76,51]], 1., atol=1e-6)

    if visualtests
      gr(size=(800,400))
      @test_ref_plot "data/GlobalKriging2D.png" contourf(solution)
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
    S = solution[:value]
    @test isapprox(S[inds[26,26]], 1., atol=1e-6)
    @test isapprox(S[inds[51,76]], 0., atol=1e-6)
    @test isapprox(S[inds[76,51]], 1., atol=1e-6)

    if visualtests
      gr(size=(800,400))
      @test_ref_plot "data/NearestKriging2D.png" contourf(solution)
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
    S = solution[:value]
    @test isapprox(S[inds[26,26]], 1., atol=1e-6)
    @test isapprox(S[inds[51,76]], 0., atol=1e-6)
    @test isapprox(S[inds[76,51]], 1., atol=1e-6)

    if visualtests
      gr(size=(800,400))
      @test_ref_plot "data/LocalKriging2D.png" contourf(solution)
    end
  end
end
