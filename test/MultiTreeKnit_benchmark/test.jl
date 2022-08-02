using TreeKnit
using TreeKnit.SplitGraph
using Test
using TreeTools

println("##### MultiTreeKnit benchmark #####")

include("$(dirname(pathof(TreeKnit)))/..//test/MultiTreeKnit_benchmark/test_constraints.jl")
include("$(dirname(pathof(TreeKnit)))/..//test/MultiTreeKnit_benchmark/test_benchmark.jl")
include("$(dirname(pathof(TreeKnit)))/..//test/MultiTreeKnit_benchmark/test_measures.jl")