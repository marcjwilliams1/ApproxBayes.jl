using ApproximateBayesianComputation
using Distributions
using Base.Test

tests = ["sampling"]

println("Running tests ...")

for t in tests
    fn = "test_$t.jl"
    println("* $fn ...")
    include(fn)
end
