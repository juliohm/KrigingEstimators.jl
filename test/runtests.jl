using KrigingEstimators
using Meshes
using GeoStatsBase
using Variography
using CoDa
using Unitful
using LinearAlgebra
using Statistics
using Test, Random

# list of tests
testfiles = [
  "estimators.jl",
  "ui.jl"
]

@testset "KrigingEstimators.jl" begin
  for testfile in testfiles
    include(testfile)
  end
end
