abstract ABCtype
abstract Prior
abstract Particle

type ABCRejection <: ABCtype

  sim_func::Function
  nparams::Int64
  ϵ::Float64
  nparticles::Int64
  constants::Float64

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
