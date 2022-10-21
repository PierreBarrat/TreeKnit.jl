module TreeKnit

# External modules
using Comonicon
using Logging
using LoggingExtras
using Parameters
using Random
using Setfield
# Personal modules
using TreeTools




include("mcc_base.jl")
export naive_mccs

include("mcc_splits.jl")

include("mcc_tools.jl")
export map_mccs, map_mccs!

include("mcc_IO.jl")
export read_mccs, write_mccs

include("resolving.jl")
export resolve!

include("SplitGraph/SplitGraph.jl")
using TreeKnit.SplitGraph

include("objects.jl")
export OptArgs

include("main.jl")
export computeMCCs, inferARG

include("SimpleReassortmentGraph/SimpleReassortmentGraph.jl")
import TreeKnit.SimpleReassortmentGraph: SRG
export SRG

include("cli.jl")

# TreeTools re-exports for docs
import TreeTools: node2tree, parse_newick
export node2tree
export parse_newick


end # module
