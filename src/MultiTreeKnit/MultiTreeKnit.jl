module MultiTreeKnit

using TreeKnit
using TreeTools
using Dagger
using Combinatorics

include("multitree_objects.jl")
export MCC_set
include("multitree_constraints.jl")
include("multitree_functions.jl")

let verbose::Bool = false, vverbose::Bool = false
	global v() = verbose
	global vv() = vverbose
	global set_verbose(v) = (verbose = v)
	global set_vverbose(v) = (vverbose = v)
end

const MTK = MultiTreeKnit
export MTK

end ##module 
