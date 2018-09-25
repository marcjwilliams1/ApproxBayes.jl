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
  SimData,
  ABCRejection,
  ABCSMC,
  ABCRejectionModel,
  ABCSMCModel,

  #functions
  ksdist,
  runabc,
  writeoutput

### source files
include("types.jl")
include("util.jl")
include("ABCalgorithm.jl")
include("sampling.jl")
include("calculateweights.jl")
#include("util.jl")
include("plots.jl")

end
