function getscales(particles, ABCsetup::ABCSMC)
  #calculate the range of parameter values (ie the scale) to use for the
  #perturbation kernel
  parameters = hcat(map(x -> x.params, particles)...)'
  scales = ((maximum(parameters, dims = 1) -
                  minimum(parameters, dims = 1)) ./ABCsetup.scalefactor)[:]

  for i in 1:length(particles)
    particles[i].scales = scales
  end

  return particles
end

particleperturbationkernel(x0, scale) = rand(Uniform(x0 - scale, x0 + scale))

function perturbparticle(particle, kernel::Kernel)
  newparticle = copyparticle(particle)
  newparams = zeros(Float64, length(newparticle.params))
  for i in 1:length(newparams)
    newparams[i] = kernel.perturb_function(newparticle.params[i], kernel.kernel_parameters)
  end
  newparticle.params = newparams
  return newparticle
end

function kernelprob(p1, p2)
    prob = 1
    for i in 1:length(p1.params)
      prob = prob * pdf(Uniform(p2.params[i] - p2.scales[i], p2.params[i] + p2.scales[i]), p1.params[i])
    end
    return prob
end

particleperturbationkernel(x0, scale) = rand(Uniform(x0 - scale, x0 + scale))
pdf_function(p1, p2, scales, i) = pdf(Uniform(p2.params[i] - scales[i], p2.params[i] + scales[i]), p1.params[i])
function calculate_kernel_paramaters(particles)
  #calculate the range of parameter values (ie the scale) to use for the
  #perturbation kernel
  parameters = hcat(map(x -> x.params, particles)...)'
  scales = ((maximum(parameters, dims = 1) -
                  minimum(parameters, dims = 1)) ./2)[:]

  return scales
end

uniformkernel = Kernel(
    particleperturbationkernel,
    pdf_function,
    calculate_kernel_paramaters,
    [],
    []
)