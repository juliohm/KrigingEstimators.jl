@testset "UI elements" begin
  grid = CartesianGrid(10, 10)
  krig = kriging_ui(grid, GaussianVariogram(), nothing, nothing, nothing)
  @test krig isa OrdinaryKriging
  krig = kriging_ui(grid, GaussianVariogram(), 0.0, nothing, nothing)
  @test krig isa SimpleKriging
  krig = kriging_ui(grid, GaussianVariogram(), nothing, 2, nothing)
  @test krig isa UniversalKriging
  krig = kriging_ui(grid, GaussianVariogram(), nothing, nothing, [x -> 1])
  @test krig isa ExternalDriftKriging
end
