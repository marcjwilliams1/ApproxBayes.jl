function getnormal(params, constants, targetdata)
  simdata = rand(Normal(params...), 100)
  ksdist(simdata, targetdata), 1
end

function getuniformdist(params, constants, targetdata)
  params = sort(params)
  simdata = rand(Uniform(params...), 100)
  ksdist(simdata, targetdata), 1
end

Random.seed!(1234)
targetdata = rand(Normal(3.5, 0.44), 100)

ABCsetup = ABCSMCModel([getnormal, getuniformdist, getnormal], [2, 2, 2], 0.1, [Prior([Uniform(0.0, 20.0), Uniform(0.0, 2.0)]), Prior([Uniform(0.0, 20.0), Uniform(0.0, 2.0)]), Prior([Uniform(0.0, 20.0), Uniform(0.0, 2.0)])]; nparticles = 100, maxiterations = 10^5)

#test model perturbation kernel
Niterations = 10^6
global m = zeros(Int64, Niterations)
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
global m = zeros(Int64, Niterations)
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

#test get proposal function
ABCsetup = ABCSMC(getnormal, 2, 0.1, Prior((Uniform(0.0, 20.0), Exponential(1.0))); nparticles = 100, maxiterations = 10^5)
Niterations = 10^6
p1 = zeros(Float64, Niterations)
p2 = zeros(Float64, Niterations)
@time for i in 1:Niterations
  p1[i], p2[i] = ApproxBayes.getproposal(ABCsetup.prior, 2)
end
@test isapprox(mean(p1), 10.0, rtol = 0.001)
@test isapprox(mean(p2), 1.0, rtol = 0.001)

#test get proposal function with models with different Priors
ABCsetup = ABCSMCModel([getnormal, getuniformdist], [2, 2], 1.0,
[Prior([Exponential(1), Uniform(0.0, 1.0)]), Prior([Uniform(0.0, 3.0), Normal(2.0, 0.1)])]; nparticles = 100, maxiterations = 10^5)
Niterations = 10^6
p1 = Float64[]
p2 = Float64[]
p3 = Float64[]
p4 = Float64[]
@time for i in 1:Niterations
  global m = rand(1:2)
  if m == 1
    x = ApproxBayes.getproposal(ABCsetup.Models[m].prior, 2)
    push!(p1, x[1])
    push!(p2, x[2])
  else
    x = ApproxBayes.getproposal(ABCsetup.Models[m].prior, 2)
    push!(p3, x[1])
    push!(p4, x[2])
  end
end
@test isapprox(mean(p1), 1.0, rtol = 0.01)
@test isapprox(mean(p2), 0.5, rtol = 0.01)
@test isapprox(mean(p3), 1.5, rtol = 0.01)
@test isapprox(mean(p4), 2.0, rtol = 0.01)
