#File which defined all the algorithms. Each algorithm takes in an ABCtype


"""
    runabc(ABCsetup::ABCtype, targetdata; progress = false, verbose = false, parallel = true)

Run ABC with ABCsetup defining the algorithm and inputs to algorithm, targetdata is the data we wish to fit the model to and will be used as an input for the simulation function defined in ABCsetup. If progress is set to `true` a progress meter will be shown. Inference will be run in parallel via multithreading if `parallel = true`. The environmental variable JULIA_NUM_THREADS needs to be set prior to launching a julia session.
"""
function runabc(ABCsetup::ABCRejection, targetdata; progress = false, verbose = false, parallel = false)

  #initalize array of particles
  particlesall = Array{ParticleRejection}(undef, ABCsetup.maxiterations)

  if progress
    p = Progress(ABCsetup.nparticles, 1, "Running ABC rejection algorithm...", 30)
  end

  if parallel
    Printf.@printf("Preparing to run in parallel on %i processors\n", nthreads())

    particles = Array{ParticleRejection}(undef, ABCsetup.maxiterations)
    distvec = zeros(Float64, ABCsetup.maxiterations) #store distances in an array
    i = Atomic{Int64}(0)
    cntr = Atomic{Int64}(0)
    @threads for its = 1:ABCsetup.maxiterations

      if i[] > ABCsetup.nparticles
        break
      end

      #get new proposal parameters
      newparams = getproposal(ABCsetup.prior, ABCsetup.nparams)
      #simulate with new parameters
      dist, out = ABCsetup.simfunc(newparams, ABCsetup.constants, targetdata)
      #keep track of all particles incase we don't reach nparticles with dist < ϵ
      particlesall[its] = ParticleRejection(newparams, dist, out)

      #if simulated data is less than target tolerance accept particle
      if dist < ABCsetup.ϵ
        particles[its] = ParticleRejection(newparams, dist, out)
        distvec[its] = dist
        atomic_add!(i, 1)
      end
      atomic_add!(cntr,1)

    end
    # Remove particles that are still #undef and corresponding distances
    idx = [isassigned(particles,ii) for ii in eachindex(particles)]
    particles = particles[idx]
    distvec = distvec[idx]
    i = length(particles)    # Number of accepted particles
    its = cntr[]    # Total number of simulations

  else
    Printf.@printf("Preparing to run in serial on %i processor\n", 1)

    particles = Array{ParticleRejection}(undef, ABCsetup.nparticles)
    distvec = zeros(Float64, ABCsetup.nparticles) #store distances in an array
    i = 1 #set particle indicator to 1
    its = 0 #keep track of number of iterations
    while (i < (ABCsetup.nparticles + 1)) & (its < ABCsetup.maxiterations)

      its += 1
      #get new proposal parameters
      newparams = getproposal(ABCsetup.prior, ABCsetup.nparams)
      #simulate with new parameters
      dist, out = ABCsetup.simfunc(newparams, ABCsetup.constants, targetdata)
      #keep track of all particles incase we don't reach nparticles with dist < ϵ
      particlesall[its] = ParticleRejection(newparams, dist, out)

      #if simulated data is less than target tolerance accept particle
      if dist < ABCsetup.ϵ
        particles[i] = ParticleRejection(newparams, dist, out)
        distvec[i] = dist
        i +=1
        if progress
          next!(p)
        end
      end
    end
    i -= 1    # Correct to total number of particels
  end

  if i < ABCsetup.nparticles
    @warn "Only accepted $(i) particles with ϵ < $(ABCsetup.ϵ). \n\tYou may want to increase ϵ or increase maxiterations. \n\t Resorting to taking the $(ABCsetup.nparticles) particles with smallest distance"
    d = map(p -> p.distance, particlesall)
    particles = particlesall[sortperm(d)[1:ABCsetup.nparticles]]
    distvec = map(p -> p.distance, particles)
  elseif i > ABCsetup.nparticles
    particles = particles[1:ABCsetup.nparticles]
    distvec = distvec[1:ABCsetup.nparticles]
  end

  out = ABCrejectionresults(particles, its, ABCsetup, distvec)
  return out
end


function runabc(ABCsetup::ABCRejectionModel, targetdata; progress = false, verbose = false)

  ABCsetup.nmodels > 1 || error("Only 1 model specified, use ABCRejection method to estimate parameters for a single model")

  #initalize array of particles
  particles = Array{ParticleRejectionModel}(undef,
                ABCsetup.Models[1].nparticles)

  i = 1 #set particle indicator to 1
  its = 0 #keep track of number of iterations
  distvec = zeros(Float64, ABCsetup.Models[1].nparticles) #store distances in an array

  if progress
    p = Progress(ABCsetup.Models[1].nparticles, 1, "Running ABC rejection algorithm...", 30)
  end

  while (i < (ABCsetup.Models[1].nparticles + 1)) & (its < ABCsetup.Models[1].maxiterations)

    its += 1
    #sample uniformly from models
    model = rand(1:ABCsetup.nmodels)
    #get new proposal parameters
    newparams = getproposal(ABCsetup.Models[model].prior, ABCsetup.Models[model].nparams)
    #simulate with new parameters
    dist, out = ABCsetup.Models[model].simfunc(newparams, ABCsetup.Models[model].constants, targetdata)

    #if simulated data is less than target tolerance accept particle
    if dist < ABCsetup.Models[1].ϵ
      particles[i] = ParticleRejectionModel(newparams, model, dist, out)
      distvec[i] = dist
      i +=1
      if progress
        next!(p)
      end
    end
  end

  i > ABCsetup.Models[1].nparticles || error("Only accepted $(i) particles with ϵ < $(ABCsetup.Models[1].ϵ). \n\tDecrease ϵ or increase maxiterations ")

  out = ABCrejectionmodelresults(particles, its, ABCsetup, distvec)
  return out
end

function runabc(ABCsetup::ABCSMC, targetdata; verbose = false, progress = false, parallel = false)

  #run first population with parameters sampled from prior
  if verbose
    println("##################################################")
    println("Use ABC rejection to get first population")
  end
  ABCrejresults = runabc(ABCRejection(ABCsetup.simfunc, ABCsetup.nparams,
                  ABCsetup.ϵ1, ABCsetup.prior; nparticles = ABCsetup.nparticles,
                  maxiterations = ABCsetup.maxiterations,
                  constants = ABCsetup.constants), targetdata,
                  progress = progress, parallel = parallel);

  oldparticles, weights = setupSMCparticles(ABCrejresults, ABCsetup)
  ABCsetup.kernel.kernel_parameters = (maximum(ABCrejresults.parameters, dims = 1) - minimum(ABCrejresults.parameters, dims = 1) ./2)[:]
  ϵ = quantile(ABCrejresults.dist, ABCsetup.α) # set new ϵ to αth quantile
  ϵvec = [ϵ] #store epsilon values
  numsims = [ABCrejresults.numsims] #keep track of number of simualtions
  particles = Array{ParticleSMC}(undef, ABCsetup.nparticles) #define particles array

  if verbose
    println("Running ABC SMC... \n")
  end

  if parallel
    Printf.@printf("Preparing to run in parallel on %i processors\n", nthreads())
  else
    Printf.@printf("Preparing to run in serial on %i processor\n", 1)
  end

  popnum = 1
  finalpop = false

  if sum(ABCrejresults.dist .< ABCsetup.ϵT) == ABCsetup.nparticles
      @warn "Target ϵ reached with ABCRejection algorithm, no need to use ABC SMC algorithm, returning ABCRejection output..."
      return ABCrejresults
  end

  while (ϵ > ABCsetup.ϵT) & (sum(numsims) < ABCsetup.maxiterations)

    if progress
      p = Progress(ABCsetup.nparticles, 1, "ABC SMC population $(popnum), new ϵ: $(round(ϵ, digits = 2))...", 30)
    end

    if parallel

      # Arrays initialised with length maxiterations to enable use of unique index ii
      particles = Array{ParticleSMC}(undef, ABCsetup.maxiterations)
      distvec = zeros(Float64, ABCsetup.maxiterations)
      i = Atomic{Int64}(0)
      its = Atomic{Int64}(0)

      @threads for ii = 1:ABCsetup.maxiterations

        if i[] > ABCsetup.nparticles
          break
        end

        j = wsample(1:ABCsetup.nparticles, weights)
        particle = oldparticles[j]
        newparticle = perturbparticle(particle, ABCsetup.kernel)
        priorp = priorprob(newparticle.params, ABCsetup.prior)
        if priorp == 0.0 #return to beginning of loop if prior probability is 0
          continue
        end

        #simulate with new parameters
        dist, out = ABCsetup.simfunc(newparticle.params, ABCsetup.constants, targetdata)

        #if simulated data is less than target tolerance accept particle
        if dist < ϵ
          particles[ii] = newparticle
          distvec[ii] = dist
          particles[ii].other = out
          particles[ii].distance = dist
          atomic_add!(i, 1)
        end

        atomic_add!(its,1)
      end

      # Remove particles that are still #undef and corresponding distances
      idx = [isassigned(particles,ii) for ii in eachindex(particles)]
      particles = particles[idx][1:ABCsetup.nparticles]
      distvec = distvec[idx][1:ABCsetup.nparticles]
      its = its[]

    else
      particles = Array{ParticleSMC}(undef, ABCsetup.nparticles)
      distvec = zeros(Float64, ABCsetup.nparticles)
      i = 1 #set particle indicator to 1
      its = 1

      while i < ABCsetup.nparticles + 1

        j = wsample(1:ABCsetup.nparticles, weights)
        particle = oldparticles[j]
        newparticle = perturbparticle(particle, ABCsetup.kernel)
        priorp = priorprob(newparticle.params, ABCsetup.prior)
        if priorp == 0.0 #return to beginning of loop if prior probability is 0
          continue
        end

        #simulate with new parameters
        dist, out = ABCsetup.simfunc(newparticle.params, ABCsetup.constants, targetdata)

        #if simulated data is less than target tolerance accept particle
        if dist < ϵ
          particles[i] = newparticle
          particles[i].other = out
          particles[i].distance = dist
          distvec[i] = dist
          i += 1
          if progress
            next!(p)
          end
        end
        its += 1
      end
    end

    particles, weights = smcweights(particles, oldparticles, ABCsetup.prior, ABCsetup.kernel)
    ABCsetup.kernel.kernel_parameters = ABCsetup.kernel.calculate_kernel_parameters(particles)
    oldparticles = particles

    if finalpop == true
      break
    end

    ϵ = quantile(distvec, ABCsetup.α)

    if ϵ < ABCsetup.ϵT
      ϵ = ABCsetup.ϵT
      push!(ϵvec, ϵ)
      push!(numsims, its)
      popnum = popnum + 1
      finalpop = true
      continue
    end

    push!(ϵvec, ϵ)
    push!(numsims, its)

    if ((( abs(ϵvec[end - 1] - ϵ )) / ϵvec[end - 1]) < ABCsetup.convergence) == true
      if verbose
        println("\nNew ϵ is within $(round(ABCsetup.convergence * 100, digits=2))% of previous population, stop ABC SMC\n")
      end
      break
    end

    popnum = popnum + 1

    if verbose
      println("##################################################")
      show(ABCSMCresults(particles, numsims, ABCsetup, ϵvec))
      println("##################################################\n")
    end
  end
  if sum(numsims) >= ABCsetup.maxiterations
    if verbose
      println("\nPassed maxiterations=$(ABCsetup.maxiterations), stop ABC SMC\n")
    end
  end

  out = ABCSMCresults(particles, numsims, ABCsetup, ϵvec)
  return out
end

"""
    runabc(ABCsetup::ABCtype, targetdata; progress = false, verbose = false)

When the SMC algorithms are used, a print out at the end of each population will be made if verbose = true.
"""
function runabc(ABCsetup::ABCSMCModel, targetdata; verbose = false, progress = false)

  ABCsetup.nmodels > 1 || error("Only 1 model specified, use ABCSMC method to estimate parameters for a single model")

  #run first population with parameters sampled from prior
  if verbose
    println("##################################################")
    println("Use ABC rejection to get first population")
  end
  ABCrejresults = runabc(ABCRejectionModel(
            map(x -> x.simfunc, ABCsetup.Models),
            map(x -> x.nparams, ABCsetup.Models),
            ABCsetup.Models[1].ϵ1,
            map(x -> x.prior, ABCsetup.Models),
            constants = map(x -> x.constants, ABCsetup.Models),
            nparticles = ABCsetup.Models[1].nparticles,
            maxiterations = ABCsetup.Models[1].maxiterations),
            targetdata);

  oldparticles, weights = setupSMCparticles(ABCrejresults, ABCsetup)
  ϵ = quantile(ABCrejresults.dist, ABCsetup.α) # set new ϵ to αth quantile
  ϵvec = [ϵ] #store epsilon values
  numsims = [ABCrejresults.numsims] #keep track of number of simulations
  particles = Array{ParticleSMCModel}(undef, ABCsetup.nparticles) #define particles array
  weights, modelprob = getparticleweights(oldparticles, ABCsetup)
  ABCsetup = modelselection_kernel(ABCsetup, oldparticles)

  modelprob = ABCrejresults.modelfreq

  if verbose
    println("Run ABC SMC \n")
  end

  popnum = 1
  finalpop = false

  if verbose
    println(ABCSMCmodelresults(oldparticles, numsims, ABCsetup, ϵvec))
  end

  if sum(ABCrejresults.dist .< ABCsetup.ϵT) == ABCsetup.nparticles
      @warn "Target ϵ reached with ABCRejection algorithm, no need to use ABC SMC algorithm, returning ABCRejection output..."
      return ABCrejresults
  end

  while (ϵ >= ABCsetup.ϵT) & (sum(numsims) <= ABCsetup.maxiterations)

    i = 1 #set particle indicator to 1
    particles = Array{ParticleSMCModel}(undef, ABCsetup.nparticles)
    distvec = zeros(Float64, ABCsetup.nparticles)
    its = 1

    if progress
      p = Progress(ABCsetup.nparticles, 1, "ABC SMC population $(popnum), new ϵ: $(round(ϵ, digits = 2))...", 30)
    end
    while i < ABCsetup.nparticles + 1

      #draw model from previous model probabilities
      mstar = wsample(1:ABCsetup.nmodels, modelprob)
      #perturb model
      mdoublestar = perturbmodel(ABCsetup, mstar, modelprob)
      # sample particle with correct model
      j = wsample(1:ABCsetup.nparticles, weights[mdoublestar, :])
      particletemp = oldparticles[j]
      #perturb particle
      newparticle = perturbparticle(particletemp, ABCsetup.Models[mdoublestar].kernel)
      #calculate priorprob
      priorp = priorprob(newparticle.params, ABCsetup.Models[mdoublestar].prior)

      if priorp == 0.0 #return to beginning of loop if prior probability is 0
        continue
      end

      #simulate with new parameters
      dist, out = ABCsetup.Models[mdoublestar].simfunc(newparticle.params, ABCsetup.Models[mdoublestar].constants, targetdata)

      #if simulated data is less than target tolerance accept particle
      if dist < ϵ
        particles[i] = newparticle
        particles[i].other = out
        particles[i].distance = dist
        distvec[i] = dist
        i += 1
        if progress
          next!(p)
        end
      end

      its += 1
    end

    particles, weights = smcweightsmodel(particles, oldparticles, ABCsetup, modelprob)
    weights, modelprob = getparticleweights(particles, ABCsetup)
    ABCsetup = modelselection_kernel(ABCsetup, particles)
    oldparticles = particles

    if finalpop == true
      break
    end

    ϵ = quantile(distvec, ABCsetup.α)

    if ϵ < ABCsetup.ϵT
      ϵ = ABCsetup.ϵT
      push!(ϵvec, ϵ)
      push!(numsims, its)
      popnum = popnum + 1
      finalpop = true
      continue
    end

    push!(ϵvec, ϵ)
    push!(numsims, its)

    if verbose
      println("##################################################")
      println(ABCSMCmodelresults(particles, numsims, ABCsetup, ϵvec))
      println("##################################################\n")
    end

    if ((( abs(ϵvec[end - 1] - ϵ )) / ϵvec[end - 1]) < ABCsetup.convergence) == true
      println("New ϵ is within $(round(ABCsetup.convergence * 100, digits = 2))% of previous population, stop ABC SMC")
      break
    end

    popnum = popnum + 1
  end

  out = ABCSMCmodelresults(particles, numsims, ABCsetup, ϵvec)
  return out
end
