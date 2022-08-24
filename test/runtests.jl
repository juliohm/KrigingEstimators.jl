using KrigingEstimators
using Meshes
using GeoStatsBase
using Variography
using CoDa
using Unitful
using LinearAlgebra
using Statistics
using Plots; gr(size=(600,400))
using MeshPlots # TODO: replace by MeshViz
using ReferenceTests, ImageIO
using Test, Random

# workaround GR warnings
ENV["GKSwstype"] = "100"

# environment settings
isCI = "CI" ∈ keys(ENV)
islinux = Sys.islinux()
visualtests = !isCI || (isCI && islinux)
datadir = joinpath(@__DIR__,"data")

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
