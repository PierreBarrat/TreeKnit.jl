module RecombTools

using TreeTools
using Distributions
using StatsBase

include("objects.jl")
include("tools.jl")
include("MCC.jl")
include("resolving.jl")
include("reading.jl")
include("SplitGraph/SplitGraph.jl")
include("Iterating/iterating.jl")

using RecombTools.SplitGraph
using RecombTools.Iterating

end