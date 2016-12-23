# some useful function for calculating distances and summary statistics


function ksdist{T<:Real, S<:Real}(x::AbstractVector{T}, y::AbstractVector{S})

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


function setupSMCparticles(ABCrejresults, ABCsetup)

  weights = ones(ABCsetup.nparticles)./ABCsetup.nparticles
  scales = collect((maximum(ABCrejresults.parameters, 1) -
                  minimum(ABCrejresults.parameters, 1) ./2)')

  particles = Array(ParticleSMC, ABCsetup.nparticles)

  for i in 1:length(particles)

    particles[i] = ParticleSMC(ABCrejresults.particles[i].params, weights[1], scales)

  end

  return particles, weights
end

function getscales(particles)

  parameters = hcat(map(x -> x.params, particles)...)'
  scales = collect((maximum(parameters, 1) -
                  minimum(parameters, 1) ./2)')

  for i in 1:length(particles)
    particles[i].scales = scales
  end

  return particles
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


function show(ABCresults::ABCrejectionmodelresults)

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
