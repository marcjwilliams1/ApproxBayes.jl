# Approximate Bayesian Computation

Package to implement Approximate Bayesian computation algorithms in the [Julia](https://julialang.org/) programming language. Package implements basic ABC rejection sampler and sequential monte carlo algorithm (ABC SMC) as in Toni. et al 2009 as well as model selection versions of both (Toni. et al 2010).

## Getting Started
Package has been tested extensively with [Julia](https://julialang.org/) v0.5.1 but should work with later versions. If there any problems please report an issue.

To download the package, once you're in a Julia session type the following command:
```
Pkg.clone("https://github.com/marcjwilliams1/ApproxBayes.jl")
```

## Examples
Below is a simple example using the package to infer the mean of a normal distribution. The first step is to create an ABC type which stores the information required to run an analysis. The first input is the simulation function which returns a distance between the simulated and target data sets, the second input is the number of parameters and the the third is the desired tolerance. The final required input is the prior distributions for the parameters, this specified as by creating an a ```Prior``` type which is an array of distribution types from [Distributions.jl]() of the same length as the number of parameters. There are some more optional parameters that are specific the the different algorithms.

First we'll load ```ApproxBayes``` and ```Distributions``` packages.

```
using ApproxBayes
using Distributions
```

Now we'll set up the simulation function, we'll use the Kolmogorov Distance as our distance measure. The simulation needs to return 2 values the first being the distance, the second value is useful if additional information from the simulation needs to be stored, here this is not the case so we'll simply return 1.
```
function normaldist(params, constants, targetdata)

  simdata = rand(Normal(params...), 1000)
  ApproxBayes.ksdist(simdata, targetdata), 1
end
```

Now we can generate some target data, we'll take 100 samples from a normal distirbution with mean = 2.0 and variance = 0.4.
```
srand(1)
p1 = 2.0
p2 = 0.4
targetdata = rand(Normal(p1, p2), 1000)
```

Now we can setup an ABCrejection type and run the inference.
```
setup = ABCRejection(normaldist,
  2,
  0.1,
  Prior([Uniform(0.0, 20.0), Uniform(0.0, 2.0)]);
  maxiterations = 10^6,
  )

# run ABC inference
rejection = runabc(setup, targetdata);

#print summary
show(rejection)
```

We can do the same with ABC SMC algorithm.
```
setup = ABCSMC(normaldist,
  2,
  0.1,
  Prior([Uniform(0.0, 20.0), Uniform(0.0, 2.0)]),
  )

smc = runabc(setup, targetdata, verbose = true, progress = true);
#print summary
show(smc)
```

If verbose and progress are set to true then a progress meter will be displayed and at the end of each population a summary will be printed.

## Acknowledgments
Acknowledge ABC sysbio
