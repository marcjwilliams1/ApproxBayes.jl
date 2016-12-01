
function getproposal(p::Prior, nparams)

  newparams = zeros(Float64, nparams)

  for i in 1:nparams
    newparams[i] = rand(Uniform(p.p[i,:]...))
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

      prob = prob * pdf(Uniform(p2.params[i] - p2.scales[i], p2.scales[i] + p2.scales[i]), p1.params[i])

    end

    return prob
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
