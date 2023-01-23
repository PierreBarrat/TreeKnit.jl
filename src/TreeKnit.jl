module TreeKnit

# External modules
using Combinatorics
using Comonicon
using Dagger
using Dates
using JSON3
using Logging
using LoggingExtras
using Parameters
using Random
using Setfield

# Personal modules
using TreeTools

import Base: get, getindex, print, copy
include("objects.jl")
export OptArgs, MCC_set

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

include("main.jl")
export run_treeknit!, inferARG

include("SimpleReassortmentGraph/SimpleReassortmentGraph.jl")
import TreeKnit.SimpleReassortmentGraph: SRG
export SRG

const date_format = "HH:MM" # for logging
include("cli.jl")

# TreeTools re-exports for docs
import TreeTools: node2tree, parse_newick
export parse_newick_string


end # module
