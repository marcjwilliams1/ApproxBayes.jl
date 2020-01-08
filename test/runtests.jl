using Distributions
using Distances
using Test
using StatsBase
using ApproxBayes
using Plots
using Random

tests = ["sampling", "parameter", "parameter_parallel", "bayesfactor", "util"]

println("Running tests ...")

@time for t in tests
    fn = "test_$t.jl"
    println("* $fn ...")
    include(fn)
end
