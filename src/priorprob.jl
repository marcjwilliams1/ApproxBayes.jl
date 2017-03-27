function priorprob(parameters::Array{Float64, 1}, prior::PriorUniform)

  pprob = 1

  for i in 1:length(parameters)
      pprob = pprob * pdf(Uniform(prior.p[i, 1], prior.p[i, 2]), parameters[i])
  end

  return pprob

end
