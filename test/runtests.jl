using KrigingEstimators
using GeoStatsBase
using Variography
using LinearAlgebra
using Statistics
using Plots, VisualRegressionTests
using Test, Pkg, Random

# workaround GR warnings
ENV["GKSwstype"] = "100"

# environment settings
isCI = "CI" ∈ keys(ENV)
islinux = Sys.islinux()
visualtests = !isCI || (isCI && islinux)
if !isCI
  Pkg.add("Gtk")
  using Gtk
end
datadir = joinpath(@__DIR__,"data")

# dummy variables for testing
include("dummy.jl")

# list of tests
testfiles = [
  "estimators.jl",
  "solvers.jl"
]

@testset "KrigingEstimators.jl" begin
  for testfile in testfiles
    include(testfile)
  end
end
