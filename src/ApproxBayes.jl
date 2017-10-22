module ApproxBayes

using Distributions
using ProgressMeter
using Gadfly
using StatsBase
using Compat
using DataFrames
using Colors

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
  plotresults,
  plotparameterposterior,
  plotmodelposterior

### source files
include("types.jl")
include("ABCalgorithm.jl")
include("sampling.jl")
include("util.jl")
include("plots.jl")

end
