using ApproxBayes
using Distributions

function getnormal(params, constants, targetdata)

  simdata = rand(Normal(params...), 100)
  ApproxBayes.ksdist(simdata, targetdata), 1
end

function getbinomial(params, constants, targetdata)

  simdata = rand(Binomial(round(params[1]), params[2]), 100)
  ApproxBayes.ksdist(simdata, targetdata), 1
end

function getpoisson(params, constants, targetdata)

  simdata = rand(Poisson(params...), 100)
  ApproxBayes.ksdist(simdata, targetdata), 1
end

function getuniformdist(params, constants, targetdata)

  simdata = rand(Uniform(params...), 100)
  ApproxBayes.ksdist(simdata, targetdata), 1
end

#generate sime synthetic data
srand(1)
targetdata = rand(Normal(2, 0.4), 100)

#setup ABC alogrithm specifications for Rejection algorithm
setup = ABCRejection(getnormal,
  2,
  0.1,
  Prior([Uniform, Uniform], [[0, 20], [0, 2.0]]);
  maxiterations = 10^6,
  )
# run ABC inference
@time resrejection = runabc(setup, targetdata);
#print summary of inference
show(resrejection)

#do the same with ABC SMC algorithm
setup = ABCSMC(getnormal,
  2,
  0.1,
  Prior([Uniform, Uniform], [[0, 20], [0, 2.0]]),
  )
@time ressmc = runabc(setup, targetdata, verbose = true);
show(ressmc)


smcefficiency = ressmc.accratio/resrejection.accratio
println()
println("SMC algorithm is $(round(smcefficiency, 2)) times more efficient")
