using Distributions
using Distances
using Test
using StatsBase
using ApproxBayes
using Plots
using Random

#useful to check what the backend is
println(backend())
tests = ["sampling", "parameter", "parameter_parallel", "bayesfactor", "util"]

println("Running tests ...")

@time for t in tests
    fn = "test_$t.jl"
    println("* $fn ...")
    include(fn)
end
