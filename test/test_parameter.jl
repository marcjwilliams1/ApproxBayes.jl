# this file will test that both algorithms correctly infer parameters from a normal distribution

function getnormal2(params, constants, targetdata)

  simdata = rand(Normal(params...), 1000)
  ApproxBayes.ksdist(simdata, targetdata), 1
end

println("Test parameters of normal distribution are inferred correctly (mean within 5% of true value)")

srand(1)
p1 = 2.0
p2 = 0.4
targetdata = rand(Normal(p1, p2), 1000)

#setup ABC alogrithm specifications for Rejection algorithm
setup = ABCRejection(getnormal2,
  2,
  0.1,
  Prior([Uniform(0.0, 20.0), Uniform(0.0, 2.0)]);
  maxiterations = 10^6,
  )
# run ABC inference
resrejection = runabc(setup, targetdata);

println("\t Check ABC rejection algorithm correctly infers parameters")

# test that mean value of posterior is within 10% of true value
@test isapprox(mean(resrejection.parameters, 1)[1], p1, rtol = 0.05)
@test isapprox(mean(resrejection.parameters, 1)[2], p2, rtol = 0.05)

#now run SMC algorithm
setup = ABCSMC(getnormal2,
  2,
  0.05,
  Prior([Uniform(0.0, 20.0), Uniform(0.0, 2.0)]),
  )
ressmc = runabc(setup, targetdata, verbose = false);

println("\t Check ABC SMC algorithm correctly infers parameters")
# test that mean value of posterior is within 10% of true value
@test isapprox(mean(ressmc.parameters, weights(ressmc.weights), 1)[1], p1, rtol = 0.05)
@test isapprox(mean(ressmc.parameters, weights(ressmc.weights), 1)[2], p2, rtol = 0.05)

#test weights sum to 1
@test isapprox(sum(ressmc.weights), 1.0, rtol = 0.0001)

println("\t Check no errors arising from plotting")
plotparameterposterior(ressmc)

#test SMC is more efficient than rejection algorithm
println("\t Check ABC SMC is more efficient than ABC rejection")
setup = ABCSMC(getnormal2,
  2,
  0.1,
  Prior([Uniform(0.0, 20.0), Uniform(0.0, 2.0)]),
  )
ressmc = runabc(setup, targetdata, verbose = false);
@test ressmc.accratio/resrejection.accratio > 1.0


#Check fallback to taking 100 particles with shortest distance works
setup = ABCRejection(getnormal2,
  2,
  0.1,
  Prior([Uniform(0.0, 20.0), Uniform(0.0, 2.0)]);
  maxiterations = 10^3,
  )
# run ABC inference
resrejection = runabc(setup, targetdata);

@test resrejection.numsims == setup.maxiterations
