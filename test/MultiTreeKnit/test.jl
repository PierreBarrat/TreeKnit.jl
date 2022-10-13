using TreeKnit
using TreeKnit.SplitGraph
using Test
using TreeTools

println("##### MultiTreeKnit #####")

include("$(dirname(pathof(TreeKnit)))/..//test/MultiTreeKnit/test_measures.jl")
include("$(dirname(pathof(TreeKnit)))/..//test/MultiTreeKnit/test_constraints.jl")
include("$(dirname(pathof(TreeKnit)))/..//test/MultiTreeKnit/test_infer_tree_pairs.jl")
include("$(dirname(pathof(TreeKnit)))/..//test/MultiTreeKnit/test_parallelTK.jl")

using Pkg
if "ARGTools" in keys(Pkg.project().dependencies)
    include("$(dirname(pathof(TreeKnit)))/..//test/MultiTreeKnit/test_random_multitrees.jl")
end