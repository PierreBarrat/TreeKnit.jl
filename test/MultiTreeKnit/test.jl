using TreeKnit
using TreeKnit.SplitGraph
using TreeKnit.MTK
using Test
using TreeTools

println("##### MultiTreeKnit #####")
include("$(dirname(pathof(TreeKnit)))/..//test/MultiTreeKnit/test_infer_tree_pairs.jl")
include("$(dirname(pathof(TreeKnit)))/..//test/MultiTreeKnit/test_parallel_MTK.jl")

using Pkg
if "ARGTools" in keys(Pkg.project().dependencies) && "MTKTools" in keys(Pkg.project().dependencies) && "StatsBase" in keys(Pkg.project().dependencies)
    include("$(dirname(pathof(TreeKnit)))/..//test/MultiTreeKnit/test_random_args.jl")
end