# Approximate Bayesian Computation

[![Build Status](https://travis-ci.org/marcjwilliams1/ApproxBayes.jl.svg?branch=master)](https://travis-ci.org/marcjwilliams1/ApproxBayes.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/marcjwilliams1/ApproxBayes.jl?branch=master&svg=true)](https://ci.appveyor.com/project/marcjwilliams1/approxbayes-jl/branch/master)
[![Coverage Status](https://coveralls.io/repos/github/marcjwilliams1/ApproxBayes.jl/badge.svg?branch=master)](https://coveralls.io/github/marcjwilliams1/ApproxBayes.jl?branch=master)
[![codecov.io](http://codecov.io/github/marcjwilliams1/ApproxBayes.jl/coverage.svg?branch=master)](http://codecov.io/github/marcjwilliams1/ApproxBayes.jl?branch=master)


Package to implement Approximate Bayesian computation algorithms in the [Julia](https://julialang.org/) programming language. Package implements basic ABC rejection sampler and sequential monte carlo algorithm (ABC SMC) as in Toni. et al 2009 as well as model selection versions of both (Toni. et al 2010).

## Getting Started
To download the package, once you're in a Julia session type the following command:
```julia
Pkg.add("ApproxBayes")
```

## Examples
Below is a simple example using the package to infer the mean of a normal distribution. The first step is to create an ABC type which stores the information required to run an analysis. The first input is the simulation function which returns a distance between the simulated and target data sets, the second input is the number of parameters and the the third is the desired tolerance. The final required input is the prior distributions for the parameters, this specified as by creating an a ```Prior``` type which is an array of distribution types from [Distributions.jl](https://github.com/JuliaStats/Distributions.jl/) of the same length as the number of parameters. There are some more optional parameters that are specific the the different algorithms.

First we'll load ```ApproxBayes``` and ```Distributions``` packages.

```julia
using ApproxBayes
using Distributions
```

Now we'll set up the simulation function, we'll use the Kolmogorov Distance as our distance measure. The simulation needs to return 2 values the first being the distance, the second value is useful if additional information from the simulation needs to be stored, here this is not the case so we'll simply return 1, for example sometimes we might want to keep the raw data generated from each simulation.
```julia
function normaldist(params, constants, targetdata)
  simdata = rand(Normal(params...), 1000)
  ApproxBayes.ksdist(simdata, targetdata), 1
end
```

Now we can generate some target data, we'll take 100 samples from a normal distirbution with mean = 2.0 and variance = 0.4.
```julia
using Random
Random.seed!(1)
p1 = 2.0
p2 = 0.4
targetdata = rand(Normal(p1, p2), 1000)
```

Now we can setup an ABCrejection type and run the inference.
```julia
setup = ABCRejection(normaldist, #simulation function
  2, # number of parameters
  0.1, #target ϵ
  Prior([Uniform(0.0, 20.0), Uniform(0.0, 2.0)]); # Prior for each of the parameters
  maxiterations = 10^6, #Maximum number of iterations before the algorithm terminates
  )

# run ABC inference
rejection = runabc(setup, targetdata)
```

We can do the same with ABC SMC algorithm.
```julia
setup = ABCSMC(normaldist, #simulation function
  2, # number of parameters
  0.1, #target ϵ
  Prior([Uniform(0.0, 20.0), Uniform(0.0, 2.0)]), #Prior for each of the parameters
  )

smc = runabc(setup, targetdata, verbose = true, progress = true)
```

### Parallelism
Parallelism is provided via multithreading. To use multithreading you'll need to set the JULIA_NUM_THREADS environmental variable before running julia (one way of doing this exporting the variable in the terminal eg `export JULIA_NUM_THREADS=1`). Then when running an ABCRejection or ABCSMC inference in parallel set the `parallel` keyword to true. For example the normal distribution example above would be run in parallel as follows:

```julia
setup = ABCSMC(normaldist, #simulation function
  2, # number of parameters
  0.1, #target ϵ
  Prior([Uniform(0.0, 20.0), Uniform(0.0, 2.0)]), #Prior for each of the parameters
  )

smc = runabc(setup, targetdata, verbose = true, progress = true, parallel = true)
```

### Optional arguments
There are more optional arguments for each of the algorithms, to see these simply use ```?ABCSMC``` in a Julia session. If verbose and progress are set to true then a progress meter will be displayed and at the end of each population a summary will be printed.

There are more examples provided in the examples directory and used as tests in the test directory. ApproxBayes.jl is also available as an option to perform Bayesian inference with differential equations in [DiffEqBayes.jl](https://github.com/JuliaDiffEq/DiffEqBayes.jl).

### Perturbation kernels
One requirement for the ABC SMC is to have a perturbation kernel. This kernel takes a sampled particle and perturbs the parameter vector in some way to explore the parameter space. Two default kernels are supplied by ApproxBayes.jl, a uniform kernel and a gaussian kernel. Both are adaptive in that the parameters specific to the kernel change as the distance decreases. For example, in the gaussian kernel the variance is calculated from the variance of the previous population. If you want to write your own kernel, take a look at `src/kernels.jl` for examples.

### Convenience functions
Also provided are some convenience functions for plotting and saving the output.

- `writeoutput(abcresults)`: This will write the output to a text file should you wish to some additional analysis or plotting using some other tools or languages.
- `plot`: Plotting recipes for use with [Plots.jl](https://github.com/JuliaPlots/Plots.jl) are provided. Just use `plot` on any ABC return type. This will plot histograms of the posterior distributions. For the model selection algorithm `plot(result::ABCSMCmodelresults)` will plot the model posterior probabilities, a second argument indexing a particular model will plot the parameter posterior distributions for that model, ie `plot(result::ABCSMCmodelresults, 1)` will plot the posterior distribution of parameters for model 1. You'll need to add the [Plots.jl](https://github.com/JuliaPlots/Plots.jl) packages yourself as it is not bundled in with `ApproxBayes.jl`.

## Acknowledgments
Some of the code was inspired by [ABC-SysBio](http://www.theosysbio.bio.ic.ac.uk/resources/abc-sysbio/).
