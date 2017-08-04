
function getproposal(p::Prior, nparams)

  newparams = zeros(Float64, nparams)

  for i in 1:nparams
    newparams[i] = rand(p.distribution[i])
  end

  return newparams
end

particleperturbationkernel(x0, scale) = rand(Uniform(x0 - scale, x0 + scale))

function perturbparticle(particle)

  newparticle = deepcopy(particle) #make copy of particle

  newparams = zeros(Float64, length(newparticle.params))

  for i in 1:length(newparams)
    newparams[i] = particleperturbationkernel(newparticle.params[i], newparticle.scales[i])
  end

  newparticle.params = newparams

  return newparticle
end

function kernel_prob(p1, p2)

    prob = 1

    for i in 1:length(p1.params)

      prob = prob * pdf(Uniform(p2.params[i] - p2.scales[i], p2.params[i] + p2.scales[i]), p1.params[i])

    end

    return prob
end

function perturbmodel(ABCsetup, mstar, modelprob)

    prob = ABCsetup.modelkern

    mprob = ones(Float64, length(modelprob))
    mprob[modelprob.==0.0] = 0.0

    nsurvivingmodels = sum(mprob)

    mprob[mprob.> 0.0] = (1 - prob) / (nsurvivingmodels - 1)
    mprob[mstar] = prob

    wsample(1:ABCsetup.nmodels, mprob)

end

function getmodelfreq(particles, ABCsetup)

  freq = zeros(Int64, ABCsetup.nmodels)
  models = map(x -> x.model, particles)

  for i in 1:ABCsetup.nmodels
    freq[i] = sum(models.==i)
  end

  return freq
end



function getmodelprob(currmodel, prevmodel, modelprob, ABCsetup)

  prob = ABCsetup.modelkern

  if currmodel == prevmodel
    return prob
  elseif sum(modelprob.>0.0) > 1
    return (1 - prob) / (sum(modelprob.>0.0) - 1)
  else
    return prob
  end

end


function smcweights(particles, oldparticles, prior)

  weights = zeros(Float64, length(particles))

  for i in 1:length(particles)
    numerator = priorprob(particles[i].params, prior)
    denom = 0.0
    for j in 1:length(particles)
      denom = denom + kernel_prob(particles[i], oldparticles[j])
    end

    weights[i] = numerator / (oldparticles[i].weight * denom)

  end

  for i in 1:length(particles)
    particles[i].weight = weights[i]
  end

  return particles, weights
end

function denom_modelfunc(p, prevmodelprob, ABCsetup)

  denom_m = 0.0
  for i in 1:ABCsetup.nmodels
    denom_m = denom_m + prevmodelprob[i] * getmodelprob(p.model, i, prevmodelprob, ABCsetup)
  end

  return denom_m

end

function denom_particlefunc(p, particles, prevmodelprob)

  denom_p = 0.0

  for i in particles
    if p.model == i.model

      denom_p = denom_p + ((i.weight * kernel_prob(p, i)) / prevmodelprob[p.model])

    end
  end

  return denom_p
end


function smcweightsmodel(particles, oldparticles, ABCsetup, prevmodelprob)

  numerator = zeros(Float64, length(particles))
  modeldenominator = zeros(Float64, length(particles))
  particledenominator = zeros(Float64, length(particles))

  #calculate numerator
  for i in 1:length(particles)

    numerator[i] = priorprob(particles[i].params, ABCsetup.Models[particles[i].model].prior)
    particledenominator[i] = denom_particlefunc(particles[i], oldparticles, prevmodelprob)
    modeldenominator[i] = denom_modelfunc(particles[i], prevmodelprob, ABCsetup)

  end

  weights = numerator ./ (modeldenominator .* particledenominator)

  weights = weights ./ sum(weights)

  for i in 1:length(particles)
    particles[i].weight = weights[i]
  end

  return particles, weights
end
