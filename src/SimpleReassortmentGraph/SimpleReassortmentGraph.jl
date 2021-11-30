module SimpleReassortmentGraph

using Parameters
using Random
using TreeKnit
using TreeTools

import Base: ==
import Base: isequal, push!, write

const SRG = SimpleReassortmentGraph
export SRG

include("objects.jl")
include("misc.jl")
include("prunegraft.jl")
include("construct.jl")
include("trees.jl")
include("IO.jl")

end # module
