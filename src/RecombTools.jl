module RecombTools

using TreeTools, ARGTools
using Distributions
using StatsBase
using DataFrames
using CSV
using Debugger
using Setfield
using JSON3


import Base.getindex, Base.==, Base.setindex!

include("tools.jl")
include("MCC.jl")
include("resolving.jl")
include("reading.jl")
include("SplitGraph/SplitGraph.jl")
using RecombTools.SplitGraph
include("mut_crossmap.jl")
include("main.jl")
include("artificialdata.jl")





end