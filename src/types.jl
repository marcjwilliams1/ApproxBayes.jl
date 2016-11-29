abstract ABCtype
abstract Prior
abstract Particle

type ABCrejection <: ABCtype

  sim_func::Function
  nparams::Int64
  ϵ::Float64
  nparticles::Int64

end


type PriorUniform <: Prior

  p::Array{Float64, 2}

end


type Particlerejection <: Particle

  params::Array{Float64, 1}

end

type SimData

  params::Array{Float64, 1}
  dist::Float64

end
