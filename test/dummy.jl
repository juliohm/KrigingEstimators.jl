import GeoStatsBase: solvesingle

# define a dummy solver for testing
@simsolver Dummy begin end
function solvesingle(problem::SimulationProblem,
                     covars::NamedTuple, solver::Dummy, preproc)
  npts = nelms(domain(problem))
  V = variables(problem)[var]
  vcat(fill(zero(V), npts÷2), fill(one(V), npts÷2))
end
