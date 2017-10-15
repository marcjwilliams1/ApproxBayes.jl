using Distributions
using Distances
using Base.Test
using StatsBase
using ApproxBayes

#run using Pkg.test("ApproxBayes")

tests = ["sampling", "parameter", "ode", "bayesfactor"]

println("Running tests ...")

for t in tests
    fn = "test_$t.jl"
    println("* $fn ...")
    include(fn)
end
