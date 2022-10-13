module TreeKnit

# External modules
using Comonicon
using Logging
using LoggingExtras
using Parameters
using Random
using Combinatorics
using Setfield
# Personal modules
using TreeTools

include("objects.jl")
export OptArgs

include("mcc_base.jl")
export naive_mccs

include("mcc_splits.jl")

include("mcc_tools.jl")

include("MultiTreeKnit/multitree_objects.jl")
include("MultiTreeKnit/MultiTreeKnit.jl")

include("mcc_measures.jl")

include("mcc_branch_labels.jl")


include("mcc_IO.jl")
export read_mccs, write_mccs

include("resolving.jl")
export resolve!

include("SplitGraph/SplitGraph.jl")
using TreeKnit.SplitGraph

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
