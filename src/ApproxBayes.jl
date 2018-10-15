module ApproxBayes

using Distributions
using ProgressMeter
using StatsBase
using RecipesBase
using Printf
using Distances
using DelimitedFiles
using Random
using Statistics
using Base.Threads

import Base.show

export
  # types
  ABCtype,
  Prior,
  Particle,
  abctype,
  ParticleRejection,
  ABCRejection,
  ABCSMC,
  ABCRejectionModel,
  ABCSMCModel,
  Kernel,
  gaussiankernel,
  uniformkernel,

  #functions
  ksdist,
  runabc,
  writeoutput

### source files
include("kernels.jl")
include("types.jl")
include("util.jl")
include("ABCalgorithm.jl")
include("sampling.jl")
include("calculateweights.jl")
include("plots.jl")

end
