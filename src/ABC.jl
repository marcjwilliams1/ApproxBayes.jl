module ABC


using Distributions

export
  # types
  abctype



### source files

include("types.jl")
include("ABCalg.jl")
include("sampling.jl")


end
