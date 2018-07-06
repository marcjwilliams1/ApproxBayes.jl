using Distributions
using Distances
using Base.Test
using StatsBase
using ApproxBayes
using Plots

#useful to check what the backend is
println(backend())
tests = ["sampling", "parameter", "bayesfactor", "util"]

println("Running tests ...")

@time for t in tests
    fn = "test_$t.jl"
    println("* $fn ...")
    include(fn)
end
