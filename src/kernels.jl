perturbation_function(prevparameter, scale) = rand(Uniform(prevparameter - scale, prevparameter + scale))

pdf_function(newparticle, prevparticle, scales, i) = pdf(Uniform(prevparticle.params[i] - scales[i], prevparticle.params[i] + scales[i]), newparticle.params[i])
function calculate_kernel_parameters(particles)
  #calculate the range of parameter values (ie the scale) to use for the
  #perturbation kernel
  parameters = hcat(map(x -> x.params, particles)...)'
  scales = ((maximum(parameters, dims = 1) -
                  minimum(parameters, dims = 1)) ./ 2.0)[:]
  return scales
end

uniformkernel = Kernel(
    perturbation_function,
    pdf_function,
    calculate_kernel_parameters
)

perturbation_function_gauss(currentparameter, stdev) = rand(Normal(currentparameter, stdev))
pdf_function_gauss(newparticle, prevparticle, stdev, i) = pdf(Normal(prevparticle.params[i], stdev[i]), newparticle.params[i])
function calculate_kernel_parameters_gauss(particles)
  #calculate the standard deviation of previous population parameters,
  # returns an array for the standard deviation of each parameter
  std(hcat(map(x -> x.params, particles)...), dims = 2)
end

gaussiankernel = Kernel(
    perturbation_function_gauss,
    pdf_function_gauss,
    calculate_kernel_parameters_gauss
)

function modelselection_kernel(ABCsetup, particles)
    #calculate the parameters for the perturbation kernel for each model

    #get the indeces of parameters which are part of
    modelindex = trues(ABCsetup.nparticles, ABCsetup.nmodels)
    for i in 1:ABCsetup.nmodels
        modelindex[:, i] = map(x -> x.model, particles) .== i
    end

    modelfreq = sum(modelindex, dims = 1)

    for i in 1:ABCsetup.nmodels
      if modelfreq[i] == 0
        ABCsetup.Models[i].kernel.kernel_parameters = [0.0]
      elseif modelfreq[i] == 1
        ABCsetup.Models[i].kernel.kernel_parameters = ABCsetup.Models[i].kernel.kernel_parameters
      else
        ABCsetup.Models[i].kernel.kernel_parameters =
        ABCsetup.Models[i].kernel.calculate_kernel_parameters(particles[modelindex[:, i]])
      end
    end
    return ABCsetup
end
