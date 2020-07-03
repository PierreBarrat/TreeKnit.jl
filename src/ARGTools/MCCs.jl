"""
	MCCs_from_arg(arg)

For each arg leave, find the set of leaves that are connected to it by full branches. 
"""
function MCCs_from_arg(arg)
	mccs = Dict{AbstractString,Array{AbstractString,1}}()
	visited = Dict(x => false for x in values(arg.nodes))
	for (s1,n1) in arg.nodes
		if n1.isleaf && !visited[n1]
			mccs[s1] = String[s1]
			visited[n1] = true
			for (s2,n2) in arg.nodes
				if n2.isleaf && !visited[n2] && is_linked_pair(n1,n2)
					push!(mccs[s1], s2)
					visited[n2] = true
				end
			end
		end
	end
	return mccs
end


"""
	coherent_subtrees(r::ARGNode)
"""
function coherent_subtrees(r::ARGNode)
	for c in r.color
		if !c
			@error "Node should exist for all segments."
		end
	end
	# 
	if r.isleaf 
		out = [r.label]
	else
		out = String[]
		for c in r.children
			if length(c.anc) == 1 && sum(c.color) == length(c.color)
				append!(out, coherent_subtrees(c))
			end
		end
	end
	return out
end

function coherent_subtrees(arg)
	out = Array{String,1}[]
	for n in values(arg.nodes)
		if sum(n.color) == length(n.color) && length(n.anc) > 1
			sb = coherent_subtrees(n)
			!isempty(sb) && push!(out, coherent_subtrees(n))
		end
	end
	return out
end

function common_splits(arg::ARG)
	@warn "This function doesn't make sense -- Common arg branch can lead to different splits"
	t = trees_from_ARG(arg, prune_singletons=false)[1]
	out = []
	for (s,n) in arg.nodes
		if !prod(n.isroot) && n.degree == arg.degree && length(n.anc) == 1
			sp = RecombTools.Split(Set(TreeTools.node_leavesclade_labels(t.lnodes[s])))
			if !in(sp, out)
				push!(out, sp)
			end
		end
	end
	return out
end

"""
Is the branch above `n` common to all trees?
"""
function is_common_branch(n::ARGNode)
	return prod(.!n.isroot) && n.degree == length(n.color) && length(n.anc) == 1
end

"""
Goal: Can I link the two nodes with fully common branches?  
Method:  
- go up from `n1` until a non-common branch is found. Store nodes found. 
- go up from `n2` until a non-common branch is found (false) or a previously stored node is found (true). 
"""
function is_linked_pair(n1::ARGNode, n2::ARGNode)
	if n1 == n2
		return true
	end
	d = Dict{String,Bool}()
	a = n1
	while is_common_branch(a)
		a = a.anc[1]
		d[a.label] = true
	end
	#
	# println(d)
	a = n2
	while is_common_branch(a)
		a = a.anc[1]
		if get(d, a.label, false)
			return true
		end
	end
	return false
end