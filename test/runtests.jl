using ApproxBayes
using Distributions
using Distances
using Base.Test

#run using Pkg.test("ApproxBayes")

tests = ["sampling", "bayesfactor", "parameter"]

println("Running tests ...")

for t in tests
    fn = "test_$t.jl"
    println("* $fn ...")
    include(fn)
end
