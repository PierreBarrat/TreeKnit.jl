module RecombTools


# To remove eventually
using DataFrames
# External modules
using LightGraphs
using Parameters
using Setfield
# Personal modules
using TreeTools




include("mcc_base.jl")
export naive_mccs, reduce_to_mcc, reduce_to_mcc!

include("mcc_splits.jl")

include("mcc_tools.jl")

include("mcc_IO.jl")
export read_mccs, write_mccs

include("resolving.jl")
export resolve!

include("SplitGraph/SplitGraph.jl")
using RecombTools.SplitGraph

include("objects.jl")
export OptArgs

include("main.jl")
export computeMCCs, computeMCCs!

include("Flu.jl")
export Flu

# TreeTools re-exports for docs
import TreeTools: node2tree, parse_newick
export node2tree
export parse_newick


end # module
