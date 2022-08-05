using TreeKnit
using TreeKnit.SplitGraph
using Test
using TreeTools

println("##### MultiTreeKnit benchmark #####")

include("$(dirname(pathof(TreeKnit)))/..//test/MultiTreeKnit/test_constraints.jl")
include("$(dirname(pathof(TreeKnit)))/..//test/MultiTreeKnit/test_benchmark.jl")
include("$(dirname(pathof(TreeKnit)))/..//test/MultiTreeKnit/test_measures.jl")