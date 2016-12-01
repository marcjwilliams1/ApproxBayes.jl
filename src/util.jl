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
