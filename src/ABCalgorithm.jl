
function runabc(ABCsetup::ABCRejection, targetdata)

  #initalize array of particles
  particles = Array(ParticleRejection, ABCsetup.nparticles)

  i = 1 #set particle indicator to 1
  its = 0 #keep track of number of iterations
  distvec = zeros(Float64, ABCsetup.nparticles) #store distances in an array

  while (i < (ABCsetup.nparticles + 1)) & (its < ABCsetup.maxiterations)

    its += 1

    #get new proposal parameters
    newparams = getproposal(ABCsetup.prior, ABCsetup.nparams)

    #simulate with new parameters
    dist = ABCsetup.simfunc(newparams, ABCsetup.constants, targetdata)

    #if simulated data is less than target tolerance accept particle
    if dist < ABCsetup.ϵ
      particles[i] = ParticleRejection(newparams)
      distvec[i] = dist
      i +=1
    end


  end

  i > ABCsetup.nparticles || error("Only accepted $(i-1) particles with ϵ < $(ABCsetup.ϵ). \n\tDecrease ϵ or increase maxiterations ")

  return ABCrejectionresults(particles, its, ABCsetup, distvec)

end


function runabc(ABCsetup::ABCSMC, targetdata)

  println("start")
  #run first population with parameters sampled from prior
  ABCrejresults = runabc(ABCRejection(ABCsetup.simfunc, ABCsetup.nparams,
                  ABCsetup.ϵ1, ABCsetup.prior; nparticles = ABCsetup.nparticles,
                  maxiterations = ABCsetup.maxiterations, constants = ABCsetup.constants), targetdata);

  oldparticles, weights = setupSMCparticles(ABCrejresults, ABCsetup)
  ϵ = quantile(ABCrejresults.dist, ABCsetup.α) # set new ϵ to αth quantile
  ϵvec = [ϵ] #store epsilon values
  numsims = [ABCrejresults.numsims] #keep track of number of simualtions
  particles = Array(ParticleSMC, ABCsetup.nparticles) #define particles array


  while (ϵ > ABCsetup.ϵT) & (sum(numsims) < ABCsetup.maxiterations)
    #println(ϵ)
    i = 1 #set particle indicator to 1
    particles = Array(ParticleSMC, ABCsetup.nparticles)
    distvec = zeros(Float64, ABCsetup.nparticles)
    its = 1
    while i < ABCsetup.nparticles + 1

      j = wsample(1:ABCsetup.nparticles, weights)
      particle = oldparticles[j]
      newparticle = perturbparticle(particle)

      priorp = priorprob(newparticle.params, ABCsetup.prior)

      if priorp == 0.0
        continue
      end

      #simulate with new parameters
      dist = ABCsetup.simfunc(newparticle.params, ABCsetup.constants, targetdata)

      #if simulated data is less than target tolerance accept particle
      if dist < ϵ
        particles[i] = newparticle
        distvec[i] = dist
        i += 1
      end

      its += 1
    end

    particles, weights = smcweights(particles, oldparticles, ABCsetup.prior)
    particles = getscales(particles)
    oldparticles = deepcopy(particles)
    ϵ = quantile(distvec, ABCsetup.α)
    push!(ϵvec, ϵ)
    push!(numsims, its)

  end

  return ABCSMCresults(particles, numsims, ABCsetup, ϵvec)

end

function setupSMCparticles(ABCrejresults, ABCsetup)

  weights = ones(ABCsetup.nparticles)./ABCsetup.nparticles
  scales = collect((maximum(ABCrejresults.parameters, 1) -
                  minimum(ABCrejresults.parameters, 1) ./2)')

  particles = Array(ParticleSMC, ABCsetup.nparticles)

  for i in 1:length(particles)

    particles[i] = ParticleSMC(ABCrejresults.particles[i].params, weights[1], scales)

  end

  return particles, weights
end

function getscales(particles)

  parameters = hcat(map(x -> x.params, particles)...)'
  scales = collect((maximum(parameters, 1) -
                  minimum(parameters, 1) ./2)')

  for i in 1:length(particles)
    particles[i].scales = scales
  end

  return particles
end
