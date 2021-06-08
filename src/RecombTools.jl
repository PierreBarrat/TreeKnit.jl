module RecombTools


# To remove eventually
using DataFrames
using Distributions
#External modules
using Parameters
using StatsBase
using Setfield
# Personal modules
using ARGTools
using TreeTools


include("mcc_base.jl")
export naive_mccs, reduce_to_mcc, reduce_to_mcc!

include("mcc_splits.jl")

include("mcc_tools.jl")

include("resolving.jl")
export resolve!, resolve_from_mccs!

include("reading.jl")

include("SplitGraph/SplitGraph.jl")
using RecombTools.SplitGraph

include("objects.jl")
export OptArgs

include("main.jl")
export computeMCCs, computeMCCs!

include("artificialdata.jl")





end
