abstract type ABCtype end
abstract type Particle end

type Prior

  distribution::Array{Sampleable, 1}
  lims::Array{Array{Real, 1}, 1}

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

  ABCSMC(sim_func::Function, nparams::Int64, ϵT::Float64, prior::Prior;
    maxiterations = 10^5,
    constants = [1.0],
    nparticles = 100,
    α = 0.3,
    ϵ1 = 10000.0,
    convergence = 0.05
    ) =
  new(sim_func, nparams, ϵ1, ϵT, nparticles, constants, maxiterations, prior, α, convergence)

end

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

type ABCSMCModel <: ABCtype

  Models::Array{ABCSMC, 1}
  nmodels::Int64
  modelkern::Float64
  nparticles::Int64
  α::Float64
  ϵT::Float64
  maxiterations::Int64
  convergence::Float64

  ABCSMCModel(sim_func::Array{Function, 1}, nparams::Array{Int64, 1}, ϵT::Float64, prior::Array{Prior, 1};
    constants = repeat([[1]], inner = length(sim_func)),
    maxiterations = 10^5,
    nparticles = 100,
    α = 0.3,
    ϵ1 = 10000.0,
    modelkern = 0.7,
    convergence = 0.05
    ) =
  new([ABCSMC(sim_func[i], nparams[i], ϵT, prior[i],  maxiterations = maxiterations, constants = constants[i], nparticles = nparticles, α = α, ϵ1 = ϵ1, convergence = convergence) for i in 1:length(sim_func)], length(sim_func), modelkern, nparticles, α, ϵT, maxiterations, convergence)

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
        parameters = push!(parameters, hcat(map(x -> x.params, particles[map(x -> x.model, particles) .== i])...)')

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
   dist::Array{Float64,1}
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
