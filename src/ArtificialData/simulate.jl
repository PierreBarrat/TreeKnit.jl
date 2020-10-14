module SimulateARG

using RecombTools, RecombTools.ARGTools
using TreeTools
using Random
using Distributions

export simulate, SimState, SimParam

let verbose::Bool = false, vverbose::Bool = false
	global v() = verbose
	global vv() = vverbose
	global set_verbose(v) = (verbose = v)
	global set_vverbose(v) = (vverbose = v)
end

let disc_t::Int64, exact_t::Float64
	global reset_discrete_t() = (disc_t = 1)
	global inc_discrete_t() = (disc_t += 1)
	global get_discrete_t() = (disc_t)
	global reset_exact_t() = (exact_t=0.)
	global inc_exact_t(t) = (exact_t += t)
	global get_exact_t() = (exact_t)
end
const global K = 2 # Two trees

# What matters is n/(N*r) where n is the size of the current ancestry
# If 2/(Nr) is small, it's likely that a split happens before a coalescence before the last two ancestors coalesce. 
struct SimParam
	N::Int64 # Global pop. size
	r::Float64 # Per individual recombination rate
	ρ::Float64 # Population reassortment rate (i.e. N*r)
	n0::Int64 # Initial sample size
	Tmax::Int64 # Maximal number of generations
end
mutable struct SimState
	arg::ARG
	node_times::Dict{String, Float64} # Exact time of each node
	found_color_root::Array{Bool,1} # Roots found for different colors. 
	pop_per_color::Array{Int64,1} # Number of current ancestors that have a given color
	eligible_for_reassortment::Array{String,1} # Nodes of the arg eligible for reassortment
	eligible_for_coalescence::Array{String,1} # 
end
mutable struct SimSnapshot # Used to record history of the simulation
	N::Real # Population size
	nc::Int64 # Number of nodes eligible for coalescence
	nr::Int64 # Number of nodes eligible for reassortment
	pop_per_color::Array{Int64,1} # Number of nodes per color
	event::Union{Symbol, Nothing} # :coa or :split 
	discrete_t::Int64 # Number of the event
	exact_t::Float64 # Time of the event
end

"""
	simulate(param::SimParam)
	simulate(N,r,n0)
"""
function simulate(param::SimParam; 
	verbose=false, 
	vverbose = false,
	prune_singletons=true,
	popvar=t->1, 
	output_history = false,
	simtype = :kingman)
	set_verbose(verbose)
	set_vverbose(vverbose)
	# 
	simstate = initiate(param)
	reset_discrete_t()
	reset_exact_t()
	simhistory = SimSnapshot[SimSnapshot(param.N, length(simstate.eligible_for_coalescence), 
		length(simstate.eligible_for_reassortment), 
		simstate.pop_per_color,
		nothing,
		get_discrete_t(),
		get_exact_t())]
	sim = true
	while sim
		N = param.N * popvar(get_exact_t())
		v() && println("Population size: $N")
		τ, etype = choose_event(param.r, N, 
			length(simstate.eligible_for_coalescence), 
			length(simstate.eligible_for_reassortment),
			simtype)
		v() && println("Event : $etype - Time $τ")
		if etype == :coa
			do_coalescence!(simstate, τ)
		elseif etype == :split
			do_split!(simstate, τ)
		end
		inc_discrete_t()
		inc_exact_t(τ)
		vv() && println("Eligible for coalescence: $(length(simstate.eligible_for_coalescence)) nodes - $(simstate.eligible_for_coalescence)")
		vv() && println("Eligible for reassortment: $(length(simstate.eligible_for_reassortment)) nodes - $(simstate.eligible_for_reassortment)")
		# Storing Snapshot
		push!(simhistory, SimSnapshot(N, length(simstate.eligible_for_coalescence), 
			length(simstate.eligible_for_reassortment), 
			simstate.pop_per_color,
			etype,
			get_discrete_t(),
			get_exact_t()))
		# Stop ? 
		if halt_condition(simstate, param.Tmax)
			sim = false
		elseif get_discrete_t() > param.Tmax # No more recomb to finish simulation
			v() && println("Stopping: Maximum number of iterations reached ($(param.Tmax)).")
			ρ = 0
		end
	end
	set_roots_ancestry!(simstate.arg)
	prune_singletons && ARGTools.prune_singletons!(simstate.arg)
	ARGTools.check_arg(simstate.arg)
	if output_history
		return simstate.arg, simhistory
	else
		return simstate.arg
	end
end
simulate(N,r,n0; 
	Tmax=1e6,
	verbose=false, 
	vverbose=false, 
	popvar=t->1, 
	prune_singletons=true,
	output_history=false,
	simtype=:kingman) = simulate(SimParam(N,r,N*r,n0,Tmax), 
										verbose=verbose, 
										vverbose=vverbose, 
										popvar=popvar, 
										prune_singletons=prune_singletons,
										output_history = output_history,
										simtype=simtype)


"""
	initiate(param::SimParam)

Create `param.n0` `ARGNode` structures with uninitialized parents. 
"""
function initiate(param::SimParam)
	arg = ARG(degree=K)
	for i in 1:param.n0
		an = ARGNode(degree=K, 
			anc = Array{Union{ARGNode,Nothing}}(nothing, 0),
			label="$(i)_0",
			isroot = zeros(Bool, K),
			isleaf = true)
		arg.nodes[an.label] = an
		arg.leaves[an.label] = an
	end

	return SimState(arg, 
			Dict(x=>0. for x in keys(arg.nodes)),
			zeros(Bool,K),
			param.n0*ones(Int64,K),
			collect(keys(arg.nodes)),
			collect(keys(arg.nodes)))
end

"""
	halt_condition(simstate::SimState, Tmax::Int64)
"""
function halt_condition(simstate::SimState, Tmax::Int64)
	if (&)(simstate.found_color_root...) 
		v() && println("Stopping: Roots for all colors have been found.")
		return true
	end
	return false
end

"""
	choose_event(r::Real, N::Real, n::Int, nr::Int)

Choose type of the next event. Return the time to the time to next event as well as its type `:coa` or `:split`.  
r, n and nr are resp. the population reassortment rate, the number of nodes available for coalescence and the number of nodes available for reassortment. 
"""
function choose_event(r::Real, N::Real, n::Int, nr::Int, simtype::Symbol)
	iTr = r*nr
	if simtype == :kingman 
		iTc = n*(n-1) /2. /N
	elseif simtype == :yule
		iTc = (n-1)/2. /N
	end

	t = Distributions.rand(Distributions.Exponential(1. /(iTr + iTc)) )
	if rand() <= iTc/(iTr + iTc)
		etype = :coa
	else
		etype = :split
	end
	return (t, etype)
end

"""
	do_coalescence!(simstate::SimState, t)

Find two nodes to coalesce in `simstate.arg`. 
"""
function do_coalescence!(simstate::SimState, t)
	# Chosing two nodes
	if length(simstate.eligible_for_coalescence) == 1
		@error "Can't coalesce singleton population"
	end
	i = rand(1:length(simstate.eligible_for_coalescence))
	j = rand(1:length(simstate.eligible_for_coalescence))
	while j==i 
		j = rand(1:length(simstate.eligible_for_coalescence))
	end
	n1 = simstate.arg.nodes[simstate.eligible_for_coalescence[i]]
	n2 = simstate.arg.nodes[simstate.eligible_for_coalescence[j]]
	t1 = t + get_exact_t() - simstate.node_times[n1.label]
	t2 = t + get_exact_t() - simstate.node_times[n2.label]
	# Coalesce
	new_node = do_coalescence!(simstate.arg, n1, n2, t1, t2, simstate)
	for (i,c) in enumerate(n1.color .& n2.color)
		simstate.pop_per_color[i] -= c # If n1 and n2 have color c, remove one 
		if simstate.pop_per_color[i] == 1 && !simstate.found_color_root[i]
			# new_node is the root for color i
			simstate.found_color_root[i] = true
			new_node.isroot[i] = true
			v() && println("Found root for color $i.")
			simstate.arg.root[i] = new_node
		end
	end
	# `simstate.eligible_for_coalescence` 
	deleteat!(simstate.eligible_for_coalescence, (min(i,j), max(i,j)))
	push!(simstate.eligible_for_coalescence, new_node.label)
	# simstate.eligible_for_reassortment
	# Note: if `new_node` is a root, I remove it from eligible_for_reassortment. This greatly simplifies the code... 
	idx = findall(x->in(x, (n1.label, n2.label)), simstate.eligible_for_reassortment)
	deleteat!(simstate.eligible_for_reassortment, idx)
	(new_node.degree > 1 && !(|)(new_node.isroot...)) && push!(simstate.eligible_for_reassortment, new_node.label)
	# Node times
	simstate.node_times[new_node.label] = simstate.node_times[n1.label] + n1.data[1].tau
	delete!(simstate.node_times, n1.label)
	delete!(simstate.node_times, n2.label)
	# `pop_per_color`
	# If new_node is the root for some colors and does not have any other color, remove it from eligible_for_coalescence
	if sum(new_node.isroot) == new_node.degree
		idx = findall(x->x==new_node.label, simstate.eligible_for_coalescence)
		deleteat!(simstate.eligible_for_coalescence, idx)
	end
	# Add new node to ARG
	simstate.arg.nodes[new_node.label] = new_node
end
"""	
	do_coalescence!(arg::ARG, n1::ARGNode, n2::ARGNode, t1, t2)

Coalesce `n1` and `n2` into a single `ARGNode`, and adds it to `arg`. 
"""
function do_coalescence!(arg::ARG, n1::ARGNode, n2::ARGNode, t1, t2, simstate)
	vv() && println("Attempting to coalesce $(n1.label) and $(n2.label)")
	# Parent node
	new_label = "internal_$(get_discrete_t())"
	new_color = n1.color .| n2.color
	new_color .*= (!).(simstate.found_color_root) # Colors of roots don't propagate up
	new_node = ARGNode(children = [n1,n2],
		anc = Array{Union{ARGNode,Nothing}}(nothing, 0),
		color = new_color,
		degree = sum(new_color),
		label=new_label,
		isroot=zeros(Bool, length(new_color)),
		isleaf=false)
	vv() && println("New node $(new_label) with color $new_color")
	# Children nodes
	for (n,t) in zip((n1,n2),(t1,t2))
		push!(n.anc, new_node)	# Pushing in case n is the root node of some color. In this case it already has an ancestor (nothing)
		push!(n.data, TreeTools.EvoData(tau=t))
		push!(n.anccolor, copy(n.color) .* new_node.color)
	end
	#
	return new_node
end

"""
"""
function do_split!(simstate::SimState, t)
	# Choose a node
	n = simstate.arg.nodes[sample(simstate.eligible_for_reassortment)]
	# Choose a color split
	nc1 = rand(1:n.degree-1)
	tmp = shuffle(findall(n.color))
	c1 = tmp[1:nc1]
	c2 = tmp[(nc1+1):end]
	# Split it backwards
	a1, a2 = do_split!(n, c1, c2, t + get_exact_t() - simstate.node_times[n.label])
	# Temp check -- 
	# for c in findall(simstate.found_color_root)
	# 	if a1.color[c] || a2.color[c] || n.color[c]
	# 		println("Split of node $(n.label) with color $(n.color).")
	# 		println(n.isroot)
	# 		println("$(a1.label) with color $(a1.color)")
	# 		println("$(a2.label) with color $(a2.color)")
	# 	end
	# end
	# node times
	for a in (a1,a2)
		simstate.node_times[a.label] = t + get_exact_t()
	end
	delete!(simstate.node_times, n.label)
	# eligible_for_reassortment
	deleteat!(simstate.eligible_for_reassortment, findfirst(x->x==n.label, simstate.eligible_for_reassortment))
	for a in (a1,a2)
		a.degree > 1 && push!(simstate.eligible_for_reassortment, a.label)
	end
	# eligible for coalescence
	deleteat!(simstate.eligible_for_coalescence, findfirst(x->x==n.label, simstate.eligible_for_coalescence))
	push!(simstate.eligible_for_coalescence, a1.label, a2.label)
	# Add ancestors to arg
	simstate.arg.nodes[a1.label] = a1
	simstate.arg.nodes[a2.label] = a2
	return a1, a2
end
function do_split!(n::ARGNode, c1::Array{Int64,1}, c2::Array{Int64,1}, t)
	vv() && println("Attempting to split node $(n.label) with color $(n.color).")
	vv() && println("Colors $c1 going left and $c2 going right.")
	# Parent nodes
	new_label1 = "internal1_$(get_discrete_t())"
	new_label2 = "internal2_$(get_discrete_t())"
	new_clr1 = ARGTools._color(c1, length(n.color))
	new_clr2 = ARGTools._color(c2, length(n.color))
	a1 = ARGNode(children = [n],
		anc = Array{Union{ARGNode,Nothing}}(nothing, 0),
		color = new_clr1,
		degree = sum(new_clr1),
		label = new_label1,
		isroot=zeros(Bool, length(new_clr1)),
		isleaf = false)
	a2 = ARGNode(children = [n],
		anc = Array{Union{ARGNode,Nothing}}(nothing, 0),
		color = new_clr2,
		degree = sum(new_clr2),
		label = new_label2,
		isroot=zeros(Bool, length(new_clr2)),
		isleaf = false)
	# Child
	push!(n.anc, a1)
	push!(n.anccolor, copy(new_clr1))
	if length(n.anccolor[end]) != 2 
		@warn "$(n.color)"
	end
	push!(n.data, EvoData(tau=t))
	push!(n.anc, a2)
	push!(n.anccolor, copy(new_clr2))
	if length(n.anccolor[end]) != 2 
		@warn "$(n.color)"
	end
	push!(n.data, EvoData(tau=t))
	return a1, a2
end

"""
	set_roots_ancestry!(arg::ARG)

For nodes that are roots for a color `c`, removes all ancestry corresponding to `c` and adds `nothing` as an ancestor for `c`.
"""
function set_roots_ancestry!(arg::ARG)
	for (c,ar) in enumerate(arg.root)
		# Cut branches for indices `todel` and color `c`
		todel = Int64[]
		for (i,clr) in enumerate(ar.anccolor)
			clr[c] && push!(todel, i)
		end
		!isempty(todel) && ARGTools.cut_branch!(ar, todel, c)
		# Adding `nothing` as ancestor
		push!(ar.anc, nothing)
		push!(ar.anccolor, ARGTools._color(c, arg.degree))
		push!(ar.data, EvoData(tau=missing))
	end
	ARGTools.prune_lone_nodes!(arg)
end


end