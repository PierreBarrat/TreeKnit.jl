module RecombTools

using CSV
using DataFrames
using Debugger
using Distributions
using JSON3
using Parameters
using StatsBase
using Setfield

using ARGTools
using TreeAlgs, TreeAlgs.CompatibilityTree
using TreeTools 


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