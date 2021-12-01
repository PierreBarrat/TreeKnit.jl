using TreeKnit
using TreeKnit.SplitGraph
using Test
using TreeTools

println("##### splitgraph #####")

include("$(dirname(pathof(TreeKnit)))/..//test/splitgraph/basic/test.jl")
include("$(dirname(pathof(TreeKnit)))/..//test/splitgraph/3solutions/test.jl")
include("$(dirname(pathof(TreeKnit)))/..//test/splitgraph/w_resolution/test.jl")
