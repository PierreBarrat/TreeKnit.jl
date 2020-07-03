module RecombTools

using TreeTools
using Distributions
using StatsBase


import Base.getindex, Base.==, Base.setindex!

include("objects.jl")
include("tools.jl")
include("MCC.jl")
include("resolving.jl")
include("reading.jl")
include("SplitGraph/SplitGraph.jl")
include("Iterating/Iterating.jl")
include("Interfacing/Interfacing.jl")
include("ArtificialData/ArtificialData.jl")
include("ArtificialData/simulate.jl")

using RecombTools.SplitGraph
using RecombTools.Iterating

end