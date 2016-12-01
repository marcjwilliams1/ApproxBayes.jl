module ABC


using Distributions

import Base.show

export
  # types
  ABCtype,
  Prior,
  Particle,
  abctype,
  PriorUniform,
  ParticleRejection,
  SimData,
  ABCRejection,
  ABCSMC,

  #functions
  ksdist,
  runabc





### source files

include("types.jl")
include("ABCalgorithm.jl")
include("sampling.jl")
include("util.jl")
include("priorprob.jl")


end
