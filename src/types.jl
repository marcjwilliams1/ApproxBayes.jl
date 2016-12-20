abstract ABCtype
abstract Prior
abstract Particle

type ParticleRejection <: Particle

  params::Array{Float64, 1}

end

type ParticleSMC <: Particle

  params::Array{Float64, 1}
  weight::Float64
  scales::Array{Float64, 1}

end

type ABCRejection <: ABCtype

  simfunc::Function
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

type ABCRejectionModel <: ABCtype

  Models::Array{ABCRejection, 1}

  ABCRejectionModel(sim_func::Array{Function, 1}, nparams::Array{Int64, 1}, ϵ::Float64, prior::Array{Prior, 1}, constants;
    maxiterations = 10000,
    nparticles = 100,
    ) =
  new(sim_func, nparams, ϵ, nparticles, constants, maxiterations, prior)
  new([ABCRejectionModel(sim_func[i], nparams[i], ϵ, nparticles, constants[i], maxiterations, prior[i]) for i in 1:length(sim_func)])

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


type PriorUniform <: Prior

  p::Array{Float64, 2}

end

type SimData

  params::Array{Float64, 1}
  dist::Float64

end


type ABCSMC <: ABCtype

  simfunc::Function
  nparams::Int64
  ϵ1::Float64
  ϵT::Float64
  nparticles::Int64
  constants::Array{Float64,1}
  maxiterations::Int64
  prior::Prior
  α::Float64

  ABCSMC(sim_func::Function, nparams::Int64, ϵT::Float64, prior::Prior;
    maxiterations = 10^5,
    constants = [1.0],
    nparticles = 100,
    α = 0.3,
    ϵ1 = 10000.0
    ) =
  new(sim_func, nparams, ϵ1, ϵT, nparticles, constants, maxiterations, prior, α)

end

type ABCSMCModel <: ABCtype

  Models::Array{ABCSMC, 1}

  ABCSMCModel(sim_func::Function, nparams::Array{Int64, 1}, ϵT::Float64, prior::Array{Prior, 1}, constants;
    maxiterations = 10^5,
    nparticles = 100,
    α = 0.3,
    ϵ1 = 10000.0
    ) =
  new(sim_func, nparams, ϵ1, ϵT, nparticles, constants, maxiterations, prior, α)
  new([ABCSMC(sim_func[i], nparams[i], ϵ1, ϵT, nparticles, constants[i], maxiterations, prior[i], α) for i in 1:length(sim_func)])

end

function show(ABCresults::ABCrejectionresults)

  upperci = zeros(Float64, size(ABCresults.parameters, 2))
  lowerci = zeros(Float64, size(ABCresults.parameters, 2))
  parametermeans = zeros(Float64, size(ABCresults.parameters, 2))
  parametermedians = zeros(Float64, size(ABCresults.parameters, 2))

  for i in 1:size(ABCresults.parameters, 2)
    parametermeans[i] = mean(ABCresults.parameters[:, i])
    parametermedians[i] = median(ABCresults.parameters[:, i])
    (lowerci[i], upperci[i]) = quantile(ABCresults.parameters[:, i], [0.025,0.975])
  end

  @printf("Number of simulations: %.2e\n", ABCresults.numsims)
  @printf("Acceptance ratio: %.2e\n\n", ABCresults.accratio)

  print("Median (95% intervals):\n")
  for i in 1:length(parametermeans)
      @printf("Parameter %d: %.2f (%.2f,%.2f)\n", i, parametermedians[i], lowerci[i], upperci[i])
  end

end

function show(ABCresults::ABCSMCresults)

  upperci = zeros(Float64, size(ABCresults.parameters, 2))
  lowerci = zeros(Float64, size(ABCresults.parameters, 2))
  parametermeans = zeros(Float64, size(ABCresults.parameters, 2))
  parametermedians = zeros(Float64, size(ABCresults.parameters, 2))

  for i in 1:size(ABCresults.parameters, 2)
    parametermeans[i] = mean(ABCresults.parameters[:, i])
    parametermedians[i] = median(ABCresults.parameters[:, i])
    (lowerci[i], upperci[i]) = quantile(ABCresults.parameters[:, i], [0.025,0.975])
  end

  @printf("Total Number of simulations: %.2e\n", sum(ABCresults.numsims))
  println("Cumulative number of simulations = $(cumsum(ABCresults.numsims))")
  @printf("Acceptance ratio: %.2e\n", ABCresults.accratio)
  println("Tolerance schedule = $(round(ABCresults.ϵ, 2))\n")

  print("Median (95% intervals):\n")
  for i in 1:length(parametermeans)
      @printf("Parameter %d: %.2f (%.2f,%.2f)\n", i, parametermedians[i], lowerci[i], upperci[i])
  end

end
