include("$(dirname(pathof(TreeKnit)))/..//test/SRG/fix_shared_singletons.jl")
include("$(dirname(pathof(TreeKnit)))/..//test/SRG/construct_arg.jl")
include("$(dirname(pathof(TreeKnit)))/..//test/SRG/test3.jl")

using Pkg
if "ARGTools" in keys(Pkg.project().dependencies)
	include("$(dirname(pathof(TreeKnit)))/..//test/SRG/random_args.jl")
end
