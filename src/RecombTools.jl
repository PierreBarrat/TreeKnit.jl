module RecombTools

using TreeTools
using Distributions
using StatsBase

include("objects.jl")
include("tools.jl")
include("MCC.jl")
include("resolving.jl")
include("reading.jl")
include("Iterating/iterating.jl")
include("SplitGraph/SplitGraph.jl")

using RecombTools.Iterating
using RecombTools.SplitGraph

end