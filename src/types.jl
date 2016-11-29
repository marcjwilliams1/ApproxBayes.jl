abstract ABCtype
abstract Prior
abstract Particle

type ABCRejection <: ABCtype

  sim_func::Function
  nparams::Int64
  ϵ::Float64
  nparticles::Int64
  constants::Array{Float64,1}
  maxiterations::Int64
  prior::Prior

  ABCRejection(sim_func::Function, nparams::Int64, ϵ::Float64, prior::Prior,
    maxiterations = 10000,
    constants = [1.0],
    nparticles = 100,
    ) =
  new(sim_func, nparams, ϵ, nparticles, constants, maxiterations, prior)

end


type PriorUniform <: Prior

  p::Array{Float64, 2}

end


type ParticleRejection <: Particle

  params::Array{Float64, 1}

end

type SimData

  params::Array{Float64, 1}
  dist::Float64

end
