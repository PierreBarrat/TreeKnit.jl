module RecombTools

using TreeTools, ARGTools
using TreeAlgs, TreeAlgs.CompatibilityTree
using Distributions
using StatsBase
using DataFrames
using CSV
using Debugger
using Setfield
using JSON3
using Parameters


import Base.getindex, Base.==, Base.setindex!

include("tools.jl")
include("MCC.jl")
export naive_mccs
export name_mcc_clades!, adjust_branchlength!, reduce_to_mcc, reduce_to_mcc!

include("resolving.jl")
export resolve!, resolve_from_mccs!, resolve_polytomy, resolve_polytomies

include("reading.jl")
include("SplitGraph/SplitGraph.jl")
using RecombTools.SplitGraph

include("mut_crossmap.jl")
include("objects.jl")
export OptArgs

include("main.jl")
export runopt

include("artificialdata.jl")





end