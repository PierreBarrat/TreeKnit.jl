using Pkg
using TreeKnit
using TreeKnit.SplitGraph
using Test
using TreeTools

println("##### run_treeknit!() #####")
include("$(dirname(pathof(TreeKnit)))/..//test/main/test_run_treeknit.jl")
include("$(dirname(pathof(TreeKnit)))/..//test/main/test_run_parallel_treeknit.jl")
include("$(dirname(pathof(TreeKnit)))/..//test/main/test_output.jl")

using Pkg
if "ARGTools" in keys(Pkg.project().dependencies) && "MTKTools" in keys(Pkg.project().dependencies) && "StatsBase" in keys(Pkg.project().dependencies)
    include("$(dirname(pathof(TreeKnit)))/..//test/main/test_random_args.jl")
end
