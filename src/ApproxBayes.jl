module ApproxBayes

using Distributions
using ProgressMeter
using Gadfly
using Compat

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
  plotresults

### source files
include("types.jl")
include("ABCalgorithm.jl")
include("sampling.jl")
include("util.jl")
include("priorprob.jl")
include("plots.jl")

end
