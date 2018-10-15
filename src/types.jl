abstract type ABCtype end
abstract type Particle end

"""
    Prior(distributions)

    Create Prior type for ABC algorithm specifying priors for each parameters. This is an array of Distribution types from Distribution.jl, each element corresponding to a parameter.
"""
mutable struct Prior
  distribution
  Prior(distributionarray) = new(tuple(distributionarray...))
end

mutable struct ParticleRejection <: Particle
  params::Array{Float64, 1}
  distance::Float64
  other::Any
end

mutable struct ParticleRejectionModel <: Particle
  params::Array{Float64, 1}
  model::Int64
  distance::Float64
  other::Any
end

mutable struct ParticleSMC <: Particle
  params::Array{Float64, 1}
  weight::Float64
  distance::Float64
  other::Any
end

mutable struct ParticleSMCModel <: Particle
  params::Array{Float64, 1}
  weight::Float64
  model::Int64
  distance::Float64
  other::Any
end

"""
    ABCRejection(sim_func::Function, nparams::Int64, ϵ::Float64, prior::Prior; <keyword arguments>)

Create an ABCRejection type which will simulate data with sim_func. nparams is the number of parameters inputted into sim_func, ϵ is the target tolerance and prior sets the priors for the parameters. sim_func needs to take in 3 values, the parameters (in an array), constants (array) and target data in that order and needs to return 2 values, the first being the distance between the target data and simulated data and the second can be anything but is useful if for example you want to record some additional information about the simulations.
...
## Arguments
- `maxiterations = 10^5`: Maximum number of samples before the ABC algorithm terminates.
- `constants = []`: Any constants needed to simulate from sim_func
- `nparticles = 100`: Number of particles (ie samples) of ABC algorithm
...
"""
mutable struct ABCRejection <: ABCtype

  simfunc::Function
  nparams::Int64
  ϵ::Float64
  nparticles::Int64
  constants::Array{Any,1}
  maxiterations::Int64
  prior::Prior

  ABCRejection(sim_func::Function, nparams::Int64, ϵ::Float64, prior::Prior;
    maxiterations = 10000,
    constants = [],
    nparticles = 100,
    ) =
  new(sim_func, nparams, ϵ, nparticles, constants, maxiterations, prior)

end

"""
    ABCRejection(sim_func::Function, nparams::Int64, ϵT::Float64, prior::Prior; <keyword arguments>)

Create an ABCSMC type which will simulate data with sim_func. nparams is the number of parameters inputted into sim_func, ϵT is the target tolerance and prior sets the priors for the parameters. sim_func needs to take in 3 values, the parameters (in an array), constants (array) and target data in that order and needs to return 2 values, the first being the distance between the target data and simulated data and the second can be anything but is useful if for example you want to record some additional information about the simulations.
...
## Arguments
- `maxiterations = 10^5`: Maximum number of samples before the ABC algorithm terminates.
- `constants = []`: Any constants needed to simulate from sim_func
- `nparticles = 100`: Number of particles (ie samples) of ABC algorithm
- `α = 0.3`: The αth quantile of population i is chosen as the ϵ for population i + 1
- `ϵ1 = 10^5`: Starting ϵ for first ABC SMC populations
- `convergence = 0.05`: ABC SMC stops when ϵ in population i + 1 is within 0.05 of populations i
- `kernel = uniformkernel`: Parameter perturbation kernel, default is a uniform distribution. `gaussiankernel` is also an option that is already available in ApproxBayes.jl. Alternatively you can code up your own kernel function. See kernels.jl for examples.
...
"""
mutable struct ABCSMC <: ABCtype

  simfunc::Function
  nparams::Int64
  ϵ1::Float64
  ϵT::Float64
  nparticles::Int64
  constants::Array{Any,1}
  maxiterations::Int64
  prior::Prior
  α::Float64
  convergence::Float64
  kernel::Kernel

  ABCSMC(sim_func::Function, nparams::Int64, ϵT::Float64, prior::Prior;
    maxiterations = 10^5,
    constants = [],
    nparticles = 100,
    α = 0.3,
    ϵ1 = 10000.0,
    convergence = 0.05,
    kernel = ApproxBayes.uniformkernel
    ) =
  new(sim_func, nparams, ϵ1, ϵT, nparticles, constants, maxiterations, prior, α, convergence, kernel)

end

"""
    ABCRejectionModel(sim_func::Array{Function, 1}, nparams::Array{Int64, 1}, ϵ::Float64, prior::Array{Prior, 1}; <keyword arguments>)

Create an ABCRejectionModel type which will create a type to run ABC with model selection. Each model is specified with a function, first input is an array of functions. nparams and priors are arrays for the number of parameters and priors for each model. each sim_func needs to take in 3 values, the parameters (in an array), constants (array) and target data in that order and needs to return 2 values, the first being the distance between the target data and simulated data and the second can be anything but is useful if for example you want to record some additional information about the simulations.
...
## Arguments
- `maxiterations = 10^5`: Maximum number of samples before the ABC algorithm terminates.
- `constants = [[]]`: Any constants needed to simulate from sim_func, needs to be an array of arrays, each one corresponding to a model function.
- `nparticles = 100`: Number of particles (ie samples) of ABC algorithm
...
"""
mutable struct ABCRejectionModel <: ABCtype

  Models::Array{ABCRejection, 1}
  nmodels::Int64

  ABCRejectionModel(sim_func::Array{Function, 1}, nparams::Array{Int64, 1}, ϵ::Float64, prior::Array{Prior, 1};
    constants = repeat([[]], outer = length(sim_func)),
    maxiterations = 10000,
    nparticles = 100,
    ) =
  new([ABCRejection(sim_func[i], nparams[i], ϵ, prior[i],  maxiterations = maxiterations, constants = constants[i], nparticles = nparticles) for i in 1:length(sim_func)], length(sim_func))

end

"""
    ABCSMCModel(sim_func::Array{Function, 1}, nparams::Array{Int64, 1}, ϵT::Float64, prior::Array{Prior, 1}; <keyword arguments>)

Create an ABCSMCModel type which will create a type to run the ABC SMC with model selection algorithm. Each model is specified with a function, first input is an array of functions. nparams and priors are arrays for the number of parameters and priors for each model, ϵT is the target tolerance. Each sim_func needs to take in 3 values, the parameters (in an array), constants (array) and target data in that order and needs to return 2 values, the first being the distance between the target data and simulated data and the second can be anything but is useful if for example you want to record some additional information about the simulations.
...
## Arguments
- `maxiterations = 10^5`: Maximum number of samples before the ABC algorithm terminates.
- `constants = []`: Any constants needed to simulate from sim_func
- `nparticles = 100`: Number of particles (ie samples) of ABC algorithm
- `α = 0.3`: The αth quantile of population i is chosen as the ϵ for population i + 1
- `ϵ1 = 10^5`: Starting ϵ for first ABC SMC populations
- `convergence = 0.05`: ABC SMC stops when ϵ in population i + 1 is within 0.05 of populations i
- `modelkern = 0.7`: Probability model stays the same in model perturbation kernel, ie 70% of the time the model perturbation kernel will leave the model the same.
...
"""
mutable struct ABCSMCModel <: ABCtype

  Models::Array{ABCSMC, 1}
  nmodels::Int64
  modelkern::Float64
  nparticles::Int64
  α::Float64
  ϵT::Float64
  maxiterations::Int64
  convergence::Float64
  other::Any

  function ABCSMCModel(sim_func::Array{Function, 1}, nparams::Array{Int64, 1}, ϵT::Float64, prior::Array{Prior, 1};
    constants = repeat([[]], outer = length(sim_func)),
    maxiterations = 10^5,
    nparticles = 100,
    α = 0.3,
    ϵ1 = 10000.0,
    modelkern = 0.7,
    convergence = 0.05,
    other = [],
    kernels = [ApproxBayes.uniformkernel for i in 1:length(sim_func)]
    )
    smcarray = [ABCSMC(sim_func[i], nparams[i], ϵT, prior[i],  maxiterations = maxiterations, constants = constants[i], nparticles = nparticles, α = α, ϵ1 = ϵ1, convergence = convergence, kernel = deepcopy(kernels[i])) for i in 1:length(sim_func)]
    nmodels = length(sim_func)
    new(smcarray, nmodels, modelkern, nparticles, α, ϵT, maxiterations, convergence, other)
  end

end

mutable struct ABCrejectionresults

  parameters::Array{Float64,2}
  accratio::Float64
  numsims::Int64
  dist::Array{Float64,1}
  particles::Array{ParticleRejection, 1}
  setup::ABCRejection

   function ABCrejectionresults(particles, its, ABCsetup, dist)
      parameters = hcat(map(x -> x.params, particles)...)'
      accratio = ABCsetup.nparticles/its
      new(parameters, accratio, its, dist, particles, ABCsetup)
   end
end

mutable struct ABCrejectionmodelresults

   parameters
   accratio::Float64
   numsims::Int64
   dist::Array{Float64,1}
   particles::Array{ParticleRejectionModel, 1}
   modelfreq::Array{Float64, 1}
   setup::ABCRejectionModel

   function ABCrejectionmodelresults(particles, its, ABCsetup, dist)
     parameters = []
     modelfreq = []
     for i in 1:ABCsetup.nmodels
        push!(modelfreq, sum(map(x -> x.model, particles) .== i))
        models = map(x -> x.model, particles)
        parameters = push!(parameters,
        hcat(map(x -> x.params, particles[map(x -> x.model, particles) .== i])...)')
     end
     accratio = ABCsetup.Models[1].nparticles/its
     modelfreq = modelfreq ./ sum(modelfreq)

     new(parameters, accratio, its, dist, particles, modelfreq, ABCsetup)
   end
end

mutable struct ABCSMCresults

  parameters::Array{Float64,2}
  weights::Array{Float64, 1}
  accratio::Float64
  numsims::Array{Int64,1}
  ϵ::Array{Float64,1}
  particles::Array{ParticleSMC, 1}
  setup::ABCSMC

   function ABCSMCresults(particles, numsims, ABCsetup, epsvec)
      parameters = hcat(map(x -> x.params, particles)...)'
      weights = map(x -> x.weight, particles)
      accratio = ABCsetup.nparticles/sum(numsims)
      new(parameters, weights, accratio, numsims, epsvec, particles, ABCsetup)
   end
end


mutable struct ABCSMCmodelresults

   parameters
   weights
   accratio::Float64
   numsims::Array{Int64, 1}
   ϵ::Array{Float64,1}
   particles::Array{ParticleSMCModel, 1}
   modelfreq::Array{Float64, 1}
   modelprob::Array{Float64, 1}
   setup::ABCSMCModel

   function ABCSMCmodelresults(particles, its, ABCsetup, dist)

     parameters = []
     weights = []
     modelfreq = []
     modelweights = map(x -> x.weight, particles)
     models = map(x -> x.model, particles)
     modelprob = []
     for i in 1:ABCsetup.nmodels
        push!(modelfreq, sum(map(x -> x.model, particles) .== i))
        push!(modelprob, sum(modelweights[models .== i]))
        parameters = push!(parameters, hcat(map(x -> x.params, particles[map(x -> x.model, particles) .== i])...)')
        weights = push!(weights, map(x -> x.weight, particles[map(x -> x.model, particles) .== i]))
     end
     accratio = ABCsetup.nparticles ./ sum(its)
     modelfreq = modelfreq ./ sum(modelfreq)

     new(parameters, weights, accratio, its, dist, particles, modelfreq, modelprob, ABCsetup)
   end
end
