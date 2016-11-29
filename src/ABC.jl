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

  #functions
  ksdist,
  runabc





### source files

include("types.jl")
include("ABCalg.jl")
include("sampling.jl")
include("util.jl")


end
