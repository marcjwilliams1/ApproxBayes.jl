module ABC


using Distributions

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
