function copyparticle(particle::ParticleSMC)
  return ParticleSMC(copy(particle.params),
  copy(particle.weight), copy(particle.distance), [1])
end

function copyparticle(particle::ParticleSMCModel)
  return ParticleSMCModel(copy(particle.params),
  copy(particle.weight), copy(particle.model), copy(particle.distance), [1])
end

function copyparticle(particle::ParticleRejection)
  return ParticleRejection(copy(particle.params), copy(particle.distance),
  [1])
end

function copyparticle(particle::ParticleRejectionModel)
  return ParticleRejectionModel(copy(particle.params), copy(particle.model),
  copy(particle.distance), [1])
end

function ksdist(x::AbstractVector{T}, y::AbstractVector{S}) where {T <: Real, S <: Real}

  #adapted from HypothesisTest.jl
  n_x, n_y = length(x), length(y)
  sort_idx = sortperm([x; y])
  pdf_diffs = [ones(n_x)/n_x; -ones(n_y)/n_y][sort_idx]
  cdf_diffs = cumsum(pdf_diffs)
  δp = maximum(cdf_diffs)
  δn = -minimum(cdf_diffs)
  δ = max(δp, δn)

  return δ
end

function setupSMCparticles(ABCrejresults::ABCrejectionresults, ABCsetup)
  #convert to SMC type after using ABC rejection
  weights = ones(ABCsetup.nparticles)./ABCsetup.nparticles

  particles = Array{ParticleSMC}(undef, ABCsetup.nparticles)

  for i in 1:length(particles)
    particles[i] = ParticleSMC(ABCrejresults.particles[i].params, weights[1], ABCrejresults.particles[i].distance,
    ABCrejresults.particles[i].other)
  end

  return particles, weights
end

function setupSMCparticles(ABCrejresults::ABCrejectionmodelresults, ABCsetup)
  #convert to SMC type after using ABC rejection
  weights = ones(ABCsetup.Models[1].nparticles)./ABCsetup.Models[1].nparticles
  particles = Array{ParticleSMCModel}(undef, ABCsetup.Models[1].nparticles)

  for i in 1:length(particles)
    particles[i] = ParticleSMCModel(ABCrejresults.particles[i].params,
    weights[1], ABCrejresults.particles[i].model,
    ABCrejresults.particles[i].distance,
    ABCrejresults.particles[i].other)
  end

  return particles, weights
end

function show(io::IO, ABCresults::ABCrejectionresults)

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

function show(io::IO, ABCresults::ABCSMCresults)

  upperci = zeros(Float64, size(ABCresults.parameters, 2))
  lowerci = zeros(Float64, size(ABCresults.parameters, 2))
  parametermeans = zeros(Float64, size(ABCresults.parameters, 2))
  parametermedians = zeros(Float64, size(ABCresults.parameters, 2))

  for i in 1:size(ABCresults.parameters, 2)
    parametermeans[i] = mean(ABCresults.parameters[:, i],
    weights(ABCresults.weights))
    parametermedians[i] = median(ABCresults.parameters[:, i],
    weights(ABCresults.weights))
    (lowerci[i], upperci[i]) = quantile(ABCresults.parameters[:, i],
    weights(ABCresults.weights),
    [0.025,0.975])
  end

  @printf("Total number of simulations: %.2e\n", sum(ABCresults.numsims))
  println("Cumulative number of simulations = $(cumsum(ABCresults.numsims))")
  @printf("Acceptance ratio: %.2e\n", ABCresults.accratio)
  println("Tolerance schedule = $(round.(ABCresults.ϵ, digits = 2))\n")

  print("Median (95% intervals):\n")
  for i in 1:length(parametermeans)
      @printf("Parameter %d: %.2f (%.2f,%.2f)\n", i, parametermedians[i], lowerci[i], upperci[i])
  end
end


function show(io::IO, ABCresults::ABCrejectionmodelresults)

  @printf("Number of simulations: %.2e\n", ABCresults.numsims)
  @printf("Acceptance ratio: %.2e\n\n", ABCresults.accratio)
  print("Model frequencies:\n")
  for j in 1:length(ABCresults.modelfreq)
    @printf("\tModel %d: %.2f\n", j, ABCresults.modelfreq[j])
  end

  print("\nParameters:\n\n")

  for j in 1:length(ABCresults.parameters)
    print("Model $j\n")

    upperci = zeros(Float64, size(ABCresults.parameters[j], 2))
    lowerci = zeros(Float64, size(ABCresults.parameters[j], 2))
    parametermeans = zeros(Float64, size(ABCresults.parameters[j], 2))
    parametermedians = zeros(Float64, size(ABCresults.parameters[j], 2))

    for i in 1:size(ABCresults.parameters[j], 2)
      parametermeans[i] = mean(ABCresults.parameters[j][:, i])
      parametermedians[i] = median(ABCresults.parameters[j][:, i])
      (lowerci[i], upperci[i]) = quantile(ABCresults.parameters[j][:, i], [0.025,0.975])
    end

    print("\tMedian (95% intervals):\n")
    for i in 1:length(parametermeans)
        @printf("\tParameter %d: %.2f (%.2f,%.2f)\n", i, parametermedians[i], lowerci[i], upperci[i])
    end
  end
end

function show(io::IO, ABCresults::ABCSMCmodelresults)

  @printf("Total number of simulations: %.2e\n", sum(ABCresults.numsims))
  println("Cumulative number of simulations = $(cumsum(ABCresults.numsims))")
  @printf("Acceptance ratio: %.2e\n\n", ABCresults.accratio)
  println("Tolerance schedule = $(round.(ABCresults.ϵ, digits = 2))\n")

  print("Model probabilities:\n")
  for j in 1:length(ABCresults.modelprob)
    @printf("\tModel %d: %.3f\n", j, ABCresults.modelprob[j])
  end

  print("\nParameters:\n\n")

  for j in 1:length(ABCresults.parameters)
    print("Model $j\n")

    upperci = zeros(Float64, size(ABCresults.parameters[j], 2))
    lowerci = zeros(Float64, size(ABCresults.parameters[j], 2))
    parametermeans = zeros(Float64, size(ABCresults.parameters[j], 2))
    parametermedians = zeros(Float64, size(ABCresults.parameters[j], 2))

    for i in 1:size(ABCresults.parameters[j], 2)
      parametermeans[i] = mean(ABCresults.parameters[j][:, i],
      weights(ABCresults.weights[j]))
      parametermedians[i] = median(ABCresults.parameters[j][:, i],
      weights(ABCresults.weights[j]))
      (lowerci[i], upperci[i]) = quantile(ABCresults.parameters[j][:, i],
      weights(ABCresults.weights[j]), [0.025,0.975])
    end

    print("\tMedian (95% intervals):\n")
    for i in 1:length(parametermeans)
        @printf("\tParameter %d: %.2f (%.2f,%.2f)\n", i, parametermedians[i], lowerci[i], upperci[i])
    end
  end
end

"""
    writeoutput(results; <keyword arguments>)

Write the results of an ABC inference to a text file. For model selection algorithms a text file with the parameters of each model will be written and a text file with model probabilities.
...
## Arguments
- `dir = ""`: Directory where the text file will be written to.
- `file= ""`: Filename to write to, default depends on the type of inference.
...
"""
function writeoutput(results::ABCSMCresults; dir = "", file = "SMC-output.txt")
  distance = map(x -> x.distance, results.particles)
  wts = map(x -> x.weight, results.particles)
  nparams =  size(results.parameters)[2]

  head = map(x -> "parameter$x\t", 1:nparams)
  append!(head, ["distance\t", "weights\n"])

  out = hcat(results.parameters, distance, wts)

  f = open(joinpath(dir, file), "w")
  write(f, "## ABC SMC algorithm\n")
  write(f, "## Number of simulations: $(sum(results.numsims))\n")
  write(f, "## Acceptance ratio: $(round(results.accratio, digits = 4))\n")
  write(f, "## Tolerance schedule : $(results.ϵ)\n")
  write(f, head...)
  writedlm(f, out)
  close(f)

end

function writeoutput(results::ABCSMCmodelresults; dir = "", file = "SMCModel-output")

  for i in 1:length(results.modelprob)
      prtcles = results.particles[map(x -> x.model, results.particles).==i]
      distance = map(x -> x.distance, prtcles)
      wts = map(x -> x.weight, prtcles)
      nparams =  size(results.parameters[i])[2]

      head = map(x -> "parameter$x\t", 1:nparams)
      append!(head, ["distance\t", "weights\n"])

      out = hcat(results.parameters[i], distance, wts)

      #write parameters to file
      f = open(joinpath(dir, "$(file)model$i.txt"), "w")
      write(f, "## ABC SMC algorithm\n")
      write(f, "## Model: $i\n")
      write(f, "## Number of simulations: $(sum(results.numsims))\n")
      write(f, "## Acceptance ratio: $(round(results.accratio, digits = 4))\n")
      write(f, "## Tolerance schedule : $(results.ϵ)\n")
      write(f, head...)
      writedlm(f, out)
      close(f)
  end

  #write model probabilities to file
  head = ["Model\t", "Probability\n"]
  f = open(joinpath(dir, "$(file)modelprobabilities.txt"), "w")
  write(f, head...)
  writedlm(f, hcat(collect(1:length(results.modelprob)), results.modelprob))
  close(f)

end

function writeoutput(results::ABCrejectionresults; dir = "", file = "Rejection-output.txt")
  distance = map(x -> x.distance, results.particles)
  nparams =  size(results.parameters)[2]

  head = map(x -> "parameter$x\t", 1:nparams)
  append!(head, ["distance\n"])

  out = hcat(results.parameters, distance)

  f = open(joinpath(dir, file), "w")
  write(f, "## ABC Rejection algorithm\n")
  write(f, "## Number of simulations: $(results.numsims)\n")
  write(f, "## Acceptance ratio: $(round(results.accratio, digits = 4))\n")
  write(f, head...)
  writedlm(f, out)
  close(f)

end

function writeoutput(results::ABCrejectionmodelresults; dir = "", file = "RejectionModel-output")

    for i in 1:length(results.modelfreq)
      prtcles = results.particles[map(x -> x.model, results.particles).==i]
      distance = map(x -> x.distance, prtcles)
      nparams =  size(results.parameters[i])[2]

      head = map(x -> "parameter$x\t", 1:nparams)
      append!(head, ["distance\n"])
      out = hcat(results.parameters[i], distance)

      f = open(joinpath(dir, "$(file)model$i.txt"), "w")
      write(f, "## ABC Rejection algorithm\n")
      write(f, "## Number of simulations: $(results.numsims)\n")
      write(f, "## Acceptance ratio: $(round(results.accratio, digits = 4))\n")
      write(f, head...)
      writedlm(f, out)
      close(f)
  end

  #write model probabilities to file
  head = ["Model\t", "Probability\n"]
  f = open(joinpath(dir, "$(file)modelprobabilities.txt"), "w")
  write(f, head...)
  writedlm(f, hcat(collect(1:length(results.modelfreq)), results.modelfreq))
  close(f)

end
