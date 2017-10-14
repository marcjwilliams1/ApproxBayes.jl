@compat abstract type ABCtype end
@compat abstract type Particle end

"""
    Prior(Array)

    Create Prior type for ABC algorithm specifying priors for each parameters. This is an array of Distribution types from Distribution.jl, each elements corresponding to a parameter.
"""
type Prior
  distribution::Array{Distributions.Distribution{Distributions.Univariate,Distributions.Continuous},1}
end

type ParticleRejection <: Particle
  params::Array{Float64, 1}
  distance::Float64
  other::Any
end

type ParticleRejectionModel <: Particle
  params::Array{Float64, 1}
  model::Int64
  distance::Float64
  other::Any
end

type ParticleSMC <: Particle
  params::Array{Float64, 1}
  weight::Float64
  scales::Array{Float64, 1}
  distance::Float64
  other::Any
end

type ParticleSMCModel <: Particle
  params::Array{Float64, 1}
  weight::Float64
  scales::Array{Float64, 1}
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
type ABCRejection <: ABCtype

  simfunc::Function
  nparams::Int64
  ϵ::Float64
  nparticles::Int64
  constants::Array{Any,1}
  maxiterations::Int64
  prior::Prior

  ABCRejection(sim_func::Function, nparams::Int64, ϵ::Float64, prior::Prior;
    maxiterations = 10000,
    constants = [1.0],
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
- `scalefactor = 2`: : Parameter for perturbation kernel for parameter values. Larger values means space will be explored more slowly but fewer particles will be perturbed outside prior range.
...
"""
type ABCSMC <: ABCtype

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
  scalefactor::Int64

  ABCSMC(sim_func::Function, nparams::Int64, ϵT::Float64, prior::Prior;
    maxiterations = 10^5,
    constants = [1.0],
    nparticles = 100,
    α = 0.3,
    ϵ1 = 10000.0,
    convergence = 0.05,
    scalefactor = 2,
    ) =
  new(sim_func, nparams, ϵ1, ϵT, nparticles, constants, maxiterations, prior, α, convergence, scalefactor)

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
type ABCRejectionModel <: ABCtype

  Models::Array{ABCRejection, 1}
  nmodels::Int64

  ABCRejectionModel(sim_func::Array{Function, 1}, nparams::Array{Int64, 1}, ϵ::Float64, prior::Array{Prior, 1};
    constants = repeat([[1]], inner = length(sim_func)),
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
- `scalefactor = 2`: Parameter for perturbation kernel for parameter values. Larger values means space will be explored more slowly but fewer particles will be perturbed outside prior range.
- `modelkern = 0.7`: Probability model stays the same in model perturbation kernel, ie 70% of the time the model perturbation kernel will leave the model the same.
...
"""
type ABCSMCModel <: ABCtype

  Models::Array{ABCSMC, 1}
  nmodels::Int64
  modelkern::Float64
  nparticles::Int64
  α::Float64
  ϵT::Float64
  maxiterations::Int64
  convergence::Float64
  scalefactor::Int64

  ABCSMCModel(sim_func::Array{Function, 1}, nparams::Array{Int64, 1}, ϵT::Float64, prior::Array{Prior, 1};
    constants = repeat([[1]], inner = length(sim_func)),
    maxiterations = 10^5,
    nparticles = 100,
    α = 0.3,
    ϵ1 = 10000.0,
    modelkern = 0.7,
    convergence = 0.05,
    scalefactor = 2,
    ) =
  new([ABCSMC(sim_func[i], nparams[i], ϵT, prior[i],  maxiterations = maxiterations, constants = constants[i], nparticles = nparticles, α = α, ϵ1 = ϵ1, convergence = convergence) for i in 1:length(sim_func)], length(sim_func), modelkern, nparticles, α, ϵT, maxiterations, convergence, scalefactor)

end

type ABCrejectionresults

  parameters::Array{Float64,2}
  accratio::Float64
  numsims::Int64
  dist::Array{Float64,1}
  particles::Array{ParticleRejection, 1}

   function ABCrejectionresults(particles, its, ABCsetup, dist)
      parameters = hcat(map(x -> x.params, particles)...)'
      accratio = ABCsetup.nparticles/its
      new(parameters, accratio, its, dist, particles)
   end
end

type ABCrejectionmodelresults

   parameters
   accratio::Float64
   numsims::Int64
   dist::Array{Float64,1}
   particles::Array{ParticleRejectionModel, 1}
   modelfreq::Array{Float64, 1}

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

     new(parameters, accratio, its, dist, particles, modelfreq)
   end
end

type ABCSMCresults

  parameters::Array{Float64,2}
  accratio::Float64
  numsims::Array{Int64,1}
  ϵ::Array{Float64,1}
  particles::Array{ParticleSMC, 1}

   function ABCSMCresults(particles, numsims, ABCsetup, epsvec)
      parameters = hcat(map(x -> x.params, particles)...)'
      accratio = ABCsetup.nparticles/sum(numsims)
      new(parameters, accratio, numsims, epsvec, particles)
   end
end


type ABCSMCmodelresults

   parameters
   accratio::Float64
   numsims::Array{Int64, 1}
   ϵ::Array{Float64,1}
   particles::Array{ParticleSMCModel, 1}
   modelfreq::Array{Float64, 1}
   modelprob::Array{Float64, 1}

   function ABCSMCmodelresults(particles, its, ABCsetup, dist)

     parameters = []
     modelfreq = []
     modelweights = map(x -> x.weight, particles)
     models = map(x -> x.model, particles)
     modelprob = []
     for i in 1:ABCsetup.nmodels
        push!(modelfreq, sum(map(x -> x.model, particles) .== i))
        push!(modelprob, sum(modelweights[models .== i]))
        parameters = push!(parameters, hcat(map(x -> x.params, particles[map(x -> x.model, particles) .== i])...)')
     end
     accratio = ABCsetup.nparticles ./ sum(its)
     modelfreq = modelfreq ./ sum(modelfreq)

     new(parameters, accratio, its, dist, particles, modelfreq, modelprob)
   end
end

type SimData
  params::Array{Float64, 1}
  dist::Float64
end
