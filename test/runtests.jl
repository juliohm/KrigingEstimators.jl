using KrigingEstimators
using GeoStatsBase
using Variography
using LinearAlgebra
using Statistics
using Plots; gr(size=(600,400))
using ReferenceTests, ImageIO
using Test, Random

# workaround GR warnings
ENV["GKSwstype"] = "100"

# environment settings
isCI = "CI" ∈ keys(ENV)
islinux = Sys.islinux()
visualtests = !isCI || (isCI && islinux)
datadir = joinpath(@__DIR__,"data")

# helper functions for visual regression tests
function asimage(plt)
  io = IOBuffer()
  show(io, "image/png", plt)
  seekstart(io)
  ImageIO.load(io)
end
macro test_ref_plot(fname, plt)
  esc(quote
    @test_reference $fname asimage($plt)
  end)
end

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
