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

  ABCRejection(sim_func::Function, nparams::Int64, ϵ::Float64, prior::Prior;
    maxiterations = 10000,
    constants = [1.0],
    nparticles = 100,
    ) =
  new(sim_func, nparams, ϵ, nparticles, constants, maxiterations, prior)

end

type ABCrejectionresults

  parameters::Array{Float64,2}
  accratio::Float64
  numsims::Float64

   function ABCrejectionresults(particles, its, ABCsetup)

      parameters = hcat(map(x -> x.params, particles)...)'
      accratio = ABCsetup.nparticles/its

      for i in 1:size(parameters, 2)
         println("Parameter $i")
         println("\t Mean = $(mean(parameters[:,i]))")
         println("\t stdev = $(std(parameters[:,i]))")
      end

      println("")
      println("Number of simulation = $its")
      println("Acceptance Ratio = $(accratio)")

      new(parameters, accratio, its)
   end
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
