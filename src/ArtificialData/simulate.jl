module SimulateARG

using RecombTools.ARGTools
using Random
using Distributions

let n::Int64, nc::Int64
	global reset_popsize(n0) = (n=n0; nc=nc)
	global inc_n() = (n+=1)
	global inc_nc() = (nc+=1)
	global get_n() = n
	global get_nc() = nc
end
const global K = 2 # Two trees


struct SimParam
	N::Int64 # Global pop. size
	r::Float64 # Per individual recombination rate
	n0::Int64 # Initial sample size
	Tmax::Int64 # Maximal number of generations
end
mutable struct SimState
	arg::ARG
	found_color_root::Array{Bool,1} # If the root of a given color is found, stop considering this color.
	eligible_for_reassortment::Array{String,1} # Nodes of the arg eligible for reassortment
end

"""
	initiate(param::SimParam)

Create `param.n0` `ARGNode` structures with uninitialized parents. 
"""
function initiate(param::SimParam)
	reset_popsize(param.n0)
	arg = ARG(degree=K)
	for i in 1:param.n0
		an = ARGNode(degree=K, 
			label="$i_0")
		arg.nodes[an.label] = an
	end
	return SimState(arg, 
			found_color_root=zeros(Bool,K),
			eligible_for_reassortment=collect(keys(an.nodes)))
end

"""
	choose_event(r, n::Int64)

Choose type of the next event. Return the time to the discrete time to next event as well as its type `:coa` or `:split`. 
"""
function choose_event(r, n::Int, nc::Int)
	iTr = r*nc
	iTc = n*(n-1)/2.
	t = Distributions.rand(Distribution.Exponential(1./(iTr + iTc)) )
	if rand() < iTc/iTr
		etype = :coa
	else
		etype = :split
	end
	return (t,etype)
end

end