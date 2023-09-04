using KrigingEstimators
using Meshes
using GeoTables
using GeoStatsBase
using Variography
using CoDa
using Unitful
using LinearAlgebra
using Statistics
using Test, Random

# list of tests
testfiles = ["estimators.jl"]

@testset "KrigingEstimators.jl" begin
  for testfile in testfiles
    include(testfile)
  end
end
