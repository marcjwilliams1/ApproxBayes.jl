function priorprob(parameters::Array{Float64, 1}, prior::Prior)

  pprob = 1

  for i in 1:length(parameters)
      pprob = pprob * pdf(prior.distribution[i], parameters[i])
  end

  return pprob

end
