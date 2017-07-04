using ApproxBayes

srand(1)

targetdata = rand(Normal(2, 0.4), 100)

function getnormal(params, constants, targetdata)

  simdata = rand(Normal(params...), 100)

  ApproxBayes.ksdist(simdata, targetdata), 1

end

function getbinomial(params, constants, targetdata)

  simdata = rand(Binomial(round(params[1]), params[2]), 100)

  ApproxBayes.ksdist(simdata, targetdata)

end

function getpoisson(params, constants, targetdata)

  simdata = rand(Poisson(params...), 100)

  ApproxBayes.ksdist(simdata, targetdata)

end

function getuniformdist(params, constants, targetdata)

  simdata = rand(Uniform(params...), 100)

  ApproxBayes.ksdist(simdata, targetdata)

end

function getnormaloneparam(params, constants, targetdata)

  m = params[1]
  sd = constants[1]

  simdata = rand(Normal(m, sd), 100)

  ApproxBayes.ksdist(simdata, targetdata)

end

function getnormalthreeparam(params, constants, targetdata)

  simdata = rand(Normal(params[1], params[2]), 100) .+ params[3]

  ApproxBayes.ksdist(simdata, targetdata)

end

#@time abcres = ABC.runabc(ABC.ABCRejection(getnormal, 2, 0.1,  ABC.PriorUniform([0 4; 0 2.0]); nparticles = 100, maxiterations = 10^7), targetdata);
cst = [[i] for i in 1:2]
#@time abcres= ABC.runabc(ABC.ABCRejectionModel([getnormal, getbinomial], [2, 2], 0.1, [ABC.PriorUniform([0 20; 0 2.0]), ABC.PriorUniform([0 100; 0 1.0])], cst; nparticles = 100, maxiterations = 10^7), targetdata);

@time res = ABC.runabc(ABC.ABCSMC(getnormal, 2, 0.1,  ABC.PriorUniform([0 20; 0 2.0])), targetdata);

#ABCSMC(getnormal, 2, 0.1,  PriorUniform([0 20; 0 2.0])
#@time res = ABC.runabc(ABC.ABCSMC(getnormal, 2, 0.1,  ABC.PriorUniform([0 20; 0 2.0])), targetdata);

#@time abcres = ABC.runabc(ABC.ABCRejection(getnormal, 2, 0.1,  ABC.PriorUniform([0 4; 0 2.0]); nparticles = 100, maxiterations = 10^7), targetdata);
