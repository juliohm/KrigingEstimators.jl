@testset "Meuse R gstat compatibility" begin
    """meuse_grid.csv constains (x, y) plus estimation/variance of log(zinc)
    for global and local (25 samples)
    """
    datadir = joinpath(@__DIR__,"data")
    data_meuse = readgeotable(joinpath(datadir,"meuse.csv"), coordnames=(:x,:y))
    grid_meuse = readgeotable(joinpath(datadir,"meuse_grid.csv"), coordnames=(:x,:y))

    dim = 2
    tol = 0.001
    knn = 25

    vmodel = SphericalVariogram(range=874.0,nugget=0.04,sill=0.59+0.04)
    problem = EstimationProblem(data_meuse, grid_meuse, :log_zinc)
    
    @testset "LocalKriging" begin

        local_kriging = Kriging(
            :log_zinc => (variogram=vmodel,maxneighbors=knn)
        )

        solution = solve(problem, local_kriging)
        est_local, var_local = solution[:log_zinc]

        pdomain = domain(problem)
        xₒ = Array{coordtype(pdomain),}(undef,dim)
        for i in 1:nelms(grid_meuse)
            coordinates!(xₒ, pdomain, i)
            tv_est_local = grid_meuse[:pred_25][i]
            tv_var_local = grid_meuse[:var_25][i]
            @test est_local[i] ≈ tv_est_local atol=tol
            @test var_local[i] ≈ tv_var_local atol=tol
        end
    end

    @testset "GlobalKriging" begin
        global_kriging = Kriging(
            :log_zinc => (variogram=vmodel,)
        )

        solution = solve(problem, global_kriging)
        est_global, var_global = solution[:log_zinc]

        pdomain = domain(problem)
        xₒ = Array{coordtype(pdomain),}(undef,2)
        for i in 1:nelms(grid_meuse)
            coordinates!(xₒ, pdomain, i)
            tv_est_global = grid_meuse[:pred][i]
            tv_var_global = grid_meuse[:var][i]
            @test est_global[i] ≈ tv_est_global atol=tol
            @test var_global[i] ≈ tv_var_global atol=tol
        end
    end
end