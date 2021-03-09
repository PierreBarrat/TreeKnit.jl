module RecombTools

using TreeTools, ARGTools
using Distributions
using StatsBase
using DataFrames
using CSV
using Debugger


import Base.getindex, Base.==, Base.setindex!

include("tools.jl")
include("MCC.jl")
# include("resolving.jl")
include("reading.jl")
include("SplitGraph/SplitGraph.jl")
include("artificialdata.jl")

using RecombTools.SplitGraph

end