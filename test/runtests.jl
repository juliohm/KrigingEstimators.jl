using KrigingEstimators
using LinearAlgebra
using Statistics
using Test, Pkg, Random

# environment settings
islinux = Sys.islinux()
istravis = "TRAVIS" ∈ keys(ENV)
datadir = joinpath(@__DIR__,"data")
visualtests = !istravis || (istravis && islinux)
if !istravis
  Pkg.add("Gtk")
  using Gtk
end

# list of tests
testfiles = [
  "estimators.jl"
]

@testset "KrigingEstimators.jl" begin
  for testfile in testfiles
    include(testfile)
  end
end
