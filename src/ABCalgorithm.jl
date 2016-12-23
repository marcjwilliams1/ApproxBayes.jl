
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

  out = ABCrejectionresults(particles, its, ABCsetup, distvec)
  return out

end


function runabc(ABCsetup::ABCRejectionModel, targetdata)

    ABCsetup.nmodels > 1 || error("Only 1 model specified, use ABCRejection method to estimate parameters for a single model")

  #initalize array of particles
  particles = Array(ParticleRejectionModel, ABCsetup.Models[1].nparticles)

  i = 1 #set particle indicator to 1
  its = 0 #keep track of number of iterations
  distvec = zeros(Float64, ABCsetup.Models[1].nparticles) #store distances in an array

  while (i < (ABCsetup.Models[1].nparticles + 1)) & (its < ABCsetup.Models[1].maxiterations)

    its += 1

    #sample uniformly from models
    model = rand(1:ABCsetup.nmodels)

    #get new proposal parameters
    newparams = getproposal(ABCsetup.Models[model].prior, ABCsetup.Models[model].nparams)

    #simulate with new parameters
    dist = ABCsetup.Models[model].simfunc(newparams, ABCsetup.Models[model].constants, targetdata)

    #if simulated data is less than target tolerance accept particle
    if dist < ABCsetup.Models[1].ϵ
      particles[i] = ParticleRejectionModel(newparams, model)
      distvec[i] = dist
      i +=1
    end


  end

  i > ABCsetup.Models[1].nparticles || error("Only accepted $(i-1) particles with ϵ < $(ABCsetup.Models[1].ϵ). \n\tDecrease ϵ or increase maxiterations ")

  out = ABCrejectionmodelresults(particles, its, ABCsetup, distvec)
  return out

end


function runabc(ABCsetup::ABCSMC, targetdata)

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

    i = 1 #set particle indicator to 1
    particles = Array(ParticleSMC, ABCsetup.nparticles)
    distvec = zeros(Float64, ABCsetup.nparticles)
    its = 1
    while i < ABCsetup.nparticles + 1

      j = wsample(1:ABCsetup.nparticles, weights)
      particle = oldparticles[j]
      newparticle = perturbparticle(particle)

      priorp = priorprob(newparticle.params, ABCsetup.prior)

      if priorp == 0.0 #return to beginning of loop if prior probability is 0
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

    println("Finished population with tolerance $(round(ϵ, 2))\n")

  end

  out = ABCSMCresults(particles, numsims, ABCsetup, ϵvec)

  return out

end





function runabc(ABCsetup::ABCSMCModel, targetdata)

  #run first population with parameters sampled from prior
  ABCrejresults = runabc(ABCRejectionModel(
            map(x -> x.simfunc, ABCsetup.Models),
            map(x -> x.nparams, ABCsetup.Models),
            ABCsetup.Models[1].ϵ1,
            map(x -> x.prior, ABCsetup.Models),
            map(x -> x.constants, ABCsetup.Models);
            nparticles = ABCsetup.Models[1].nparticles,
            maxiterations = ABCsetup.Models[1].maxiterations),
            targetdata);

  oldparticles, weights = setupSMCparticles(ABCrejresults, ABCsetup)
  ϵ = quantile(ABCrejresults.dist, ABCsetup.α) # set new ϵ to αth quantile
  ϵvec = [ϵ] #store epsilon values
  numsims = [ABCrejresults.numsims] #keep track of number of simualtions
  particles = Array(ParticleSMC, ABCsetup.nparticles) #define particles array

  while (ϵ > ABCsetup.ϵT) & (sum(numsims) < ABCsetup.maxiterations)

    i = 1 #set particle indicator to 1
    particles = Array(ParticleSMC, ABCsetup.nparticles)
    distvec = zeros(Float64, ABCsetup.nparticles)
    its = 1
    while i < ABCsetup.nparticles + 1

      j = wsample(1:ABCsetup.nparticles, weights)
      particle = oldparticles[j]
      newparticle = perturbparticle(particle)

      priorp = priorprob(newparticle.params, ABCsetup.prior)

      if priorp == 0.0 #return to beginning of loop if prior probability is 0
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

    println("Finished population with tolerance $(round(ϵ, 2))\n")

  end

  out = ABCSMCresults(particles, numsims, ABCsetup, ϵvec)

  return out

end
