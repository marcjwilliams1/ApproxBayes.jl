function getnormal(params, constants, targetdata)
  simdata = rand(Normal(params...), 100)
  ksdist(simdata, targetdata), 1
end

function getuniformdist(params, constants, targetdata)
  params = sort(params)
  simdata = rand(Uniform(params...), 100)
  ksdist(simdata, targetdata), 1
end

srand(1234)
targetdata = rand(Normal(3.5, 0.44), 100)

ABCsetup = ABCSMCModel([getnormal, getuniformdist, getnormal], [2, 2, 2], 0.1, [Prior([Uniform(0.0, 20.0), Uniform(0.0, 2.0)]), Prior([Uniform(0.0, 20.0), Uniform(0.0, 2.0)]), Prior([Uniform(0.0, 20.0), Uniform(0.0, 2.0)])]; nparticles = 100, maxiterations = 10^5)

#test model perturbation kernel
Niterations = 10^6
m = zeros(Int64, Niterations )
mstar = 1
modelprob = [1/3, 1/3, 1/3]
for i in 1:Niterations
  mdoublestar = ApproxBayes.perturbmodel(ABCsetup, mstar, modelprob)
  m[i] = mdoublestar
end

@test isapprox(ABCsetup.modelkern, sum(m.==1)/length(m), rtol = 0.01)

@test isapprox((1 - ABCsetup.modelkern)./(ABCsetup.nmodels - 1), sum(m.==2)/length(m), rtol = 0.01)

#test model perturbation kernel when one model has died out
Niterations = 10^6
m = zeros(Int64, Niterations )
mstar = 1
modelprob = [0.5, 0.0, 0.5]
for i in 1:Niterations
  mdoublestar = ApproxBayes.perturbmodel(ABCsetup, mstar, modelprob)
  m[i] = mdoublestar
end

@test isapprox(ABCsetup.modelkern, sum(m.==1)/length(m), rtol = 0.01)
@test isapprox(0.0, sum(m.==2)/length(m))
@test isapprox(1.0 - ABCsetup.modelkern, sum(m.==3)/length(m), rtol = 0.01)


#test get particle weights
ABCrejresults = runabc(ABCRejectionModel(
          map(x -> x.simfunc, ABCsetup.Models),
          map(x -> x.nparams, ABCsetup.Models),
          ABCsetup.Models[1].ϵ1,
          map(x -> x.prior, ABCsetup.Models);
          nparticles = ABCsetup.Models[1].nparticles,
          maxiterations = ABCsetup.Models[1].maxiterations),
          targetdata);

oldparticles, weightsA = ApproxBayes.setupSMCparticles(ABCrejresults, ABCsetup)
weightsA, modelprob = ApproxBayes.getparticleweights(oldparticles, ABCsetup)

#test if modelprob is the same as from ABCrejresults
@test modelprob[:] ≈ ABCrejresults.modelfreq

for i in 1:ABCsetup.nmodels
  @test (map(x -> x.model, oldparticles).==i) == (weightsA[i, :].> 0.0)
end
