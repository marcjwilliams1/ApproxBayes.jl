function priorprob(parameters::Array{Float64, 1}, prior::Prior)
  pprob = 1
  for i in 1:length(parameters)
      pprob = pprob * pdf(prior.distribution[i], parameters[i])
  end
  return pprob
end

function kernelprob(p1, p2, kernel::Kernel)
    prob = 1
    for i in 1:length(p1.params)
      prob = prob * kernel.pdf_function(p1, p2, kernel.kernel_parameters, i)
    end
    return prob
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

function smcweights(particles, oldparticles, prior, kernel::Kernel)

  weights = zeros(Float64, length(particles))

  for i in 1:length(particles)
    numerator = priorprob(particles[i].params, prior)
    denominator = 0.0
    for j in 1:length(particles)
      denominator = denominator + kernelprob(particles[i], oldparticles[j], kernel)
    end
    weights[i] = numerator / (oldparticles[i].weight * denominator)
  end

  weights = weights ./ sum(weights)
  for i in 1:length(particles)
    particles[i].weight = weights[i]
  end

  return particles, weights
end

#the below functions are related to the SMC with model selection algorithm

function modelperturbation(p, prevmodelprob, ABCsetup)
  denominator_m = 0.0
  for i in 1:ABCsetup.nmodels
    denominator_m = denominator_m + prevmodelprob[i] * getmodelprob(p.model, i, prevmodelprob, ABCsetup)
  end
  return denominator_m
end

function parameterperturbation(p, particles, prevmodelprob, kernel)
  denominator_p = 0.0
  for i in particles
    if p.model == i.model
      denominator_p = denominator_p + ((i.weight * kernelprob(p, i, kernel)) / prevmodelprob[p.model])
    end
  end
  return denominator_p
end

function smcweightsmodel(particles, oldparticles, ABCsetup, prevmodelprob)

  #calculate weights of each particle
  numerator = zeros(Float64, length(particles))
  modeldenominator = zeros(Float64, length(particles))
  particledenominator = zeros(Float64, length(particles))

  for i in 1:length(particles)
    # get prior probabilty of parameter set
    numerator[i] = priorprob(particles[i].params, ABCsetup.Models[particles[i].model].prior)
    # get probability of parameter perturbation
    particledenominator[i] = parameterperturbation(particles[i], oldparticles, prevmodelprob, ABCsetup.Models[particles[i].model].kernel)
    # get probability of model perturbation
    modeldenominator[i] = modelperturbation(particles[i], prevmodelprob, ABCsetup)
  end

  weights = numerator ./ (modeldenominator .* particledenominator)
  weights = weights ./ sum(weights)

  for i in 1:length(particles)
    particles[i].weight = weights[i]
  end

  return particles, weights
end

function getparticleweights(particles, ABCsetup)
  w = zeros(Float64, ABCsetup.nmodels, ABCsetup.nparticles)
  for i in 1:ABCsetup.nparticles
    w[particles[i].model, i] = particles[i].weight
  end
  weights = w ./ sum(w, dims = 2)
  return weights, sum(w, dims = 2)[:]
end
