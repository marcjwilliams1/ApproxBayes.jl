# this file will test that both algorithms correctly infer parameters from a normal distribution

function getnormal2(params, constants, targetdata)

  simdata = rand(Normal(params...), 1000)
  ApproxBayes.ksdist(simdata, targetdata), 1

end

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

# test that mean value of posterior is within 10% of true value
@test isapprox(mean(ressmc.parameters, 1)[1], p1, rtol = 0.05)
@test isapprox(mean(ressmc.parameters, 1)[2], p2, rtol = 0.05)


#test SMC is more efficient than rejection algorithm
setup = ABCSMC(getnormal2,
  2,
  0.1,
  Prior([Uniform(0.0, 20.0), Uniform(0.0, 2.0)]),
  )
ressmc = runabc(setup, targetdata, verbose = false);

@test ressmc.accratio/resrejection.accratio > 1.0
