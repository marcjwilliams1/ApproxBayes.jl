"""
    Kernel(perturbation_function::Function,
    pdf_function::Function,
    calculate_kernel_parameters::Function)

Create a parameter perturbation kernel. Required inputs are 3 functions.
First is the `perturbation_function` which should take 2 parameters,
the parameter to be perturbed and any kernel specific parameter
(for example the standard deviation of a normal distribution if this is the kernel of choice).
Second function is the `pdf_function`, that requires 4 inputs:
    1) the newparticle
    2) the old particle
    3) kernel specific parameters and
    4) an index i.
The third function is `calculate_kernel_parameters` which given an array of
particles should calculate the kernel specific parameters for the next population.
Should you wish to keep the same parameters throughout you can just write a function that returns a number(s).
"""
mutable struct Kernel
    perturbation_function::Function
    pdf_function::Function
    calculate_kernel_parameters::Function
    kernel_parameters

    Kernel(perturbation_function::Function,
    pdf_function::Function,
    calculate_kernel_parameters::Function) =
    new(perturbation_function,
    pdf_function,
    calculate_kernel_parameters,
    [])
end

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

    #get the indeces of particles that are a member of each model
    modelindex = trues(ABCsetup.nparticles, ABCsetup.nmodels)
    for i in 1:ABCsetup.nmodels
        modelindex[:, i] = map(x -> x.model, particles) .== i
    end

    #calculate the frequency of each model
    modelfreq = sum(modelindex, dims = 1)

    #calculate the kernel parameters for each model
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
