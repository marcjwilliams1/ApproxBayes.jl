function getnormal(params, constants, targetdata)

  simdata = rand(Normal(params...), 100)

  ABC.ksdist(simdata, targetdata), 1

end

function getuniformdist(params, constants, targetdata)

  params = sort(params)

  simdata = rand(Uniform(params...), 100)

  ABC.ksdist(simdata, targetdata), 1

end
cst = [[i] for i in 1:2]
targetdata = rand(Normal(3.5, 0.44), 100)

ABCsetup = ABC.ABCSMCModel([getnormal, getuniformdist], [2, 2], 0.1, [ABC.PriorUniform([0.0 20; 0.0 2.0]), ABC.PriorUniform([0.0 20.0; 0.0 20.0])], cst; nparticles = 100, maxiterations = 10^5)

#test model perturbation kernel
Niterations = 10^5
m = zeros(Int64, Niterations )
mstar = 1
modelprob = [0.5, 0.5]
for i in 1:Niterations
  mdoublestar = ABC.perturbmodel(ABCsetup, mstar, modelprob)
  m[i] = mdoublestar
end

@test ABCsetup.modelkern ≈ sum(m.==1)/length(m)

@test_approx_eq_eps ABCsetup.modelkern sum(m.==1)/length(m)
