###############################################################################################################
########################################## Basic resolve functions ############################################
###############################################################################################################

"""
	resolve!(t::Tree, S::SplitList; conflict=:fail, usemask=false, tau=0.)

Add splits in `S` to `t` by introducing internal nodes.
New nodes are assigned a time `tau` (`0` by default).
If `conflict != :ignore`, will fail if a split `s` in `S` is not compatible with `t`.
Otherwise, silently skip the conflicting splits.
Add shared identity to nodes or `shared_map` if given as input.
"""
function resolve!(
	t::Tree{T}, S::SplitList;
	conflict=:fail, usemask=false, tau=0., safe=false, shared_map = nothing
) where T
	# Label for created nodes
	label_i = parse(Int64, TreeTools.create_label(t, "RESOLVED")[10:end])
	#
	tsplits = SplitList(t)
	for (i,s) in enumerate(S)
		if !safe && !in(s, tsplits; usemask)
			if iscompatible(s, tsplits; usemask)
				roots = TreeTools.blca([t.lleaves[x] for x in leaves(S,i)]...)
				R = lca(roots)
				# Creating a new node with `roots` as children and `r` as ancestor.
				nr = TreeNode(T(), label="RESOLVED_$(label_i)")
				label_i += 1
				for r in roots
					TreeTools.prunenode!(r)
					TreeTools.graftnode!(nr,r)
				end
				TreeTools.graftnode!(R, nr, tau=tau)
				if typeof(R.data) != TreeTools.EmptyData && haskey(R.data.dat, "shared_branch") && all([haskey(r.data.dat, "shared_branch") for r in roots])
					if R.data.dat["shared_branch"] && any([r.data.dat["shared_branch"] for r in roots])
						nr.data.dat["shared_branch"] = true
					else
						nr.data.dat["shared_branch"] = false
					end
				end
				if !isnothing(shared_map)
					if shared_map[R.label] && any([shared_map[r.label] for r in roots])
						shared_map[nr.label] = true
					else
						shared_map[nr.label] = false
					end
				end
				push!(tsplits.splits, s)
			elseif conflict != :ignore
				error("Tried to resolve tree with an incompatible split.")
			end
		end
	end
	node2tree!(t, t.root)
end

"""
	resolve(trees::Dict{T, <:Tree}, splits::Dict{T, <:SplitList}; kwargs...) where T

Resolve `trees[s]` with splits in `splits[s]` by calling `resolve!`.
Return dictionary of resolved trees. `trees` and `splits` must share keys.
This is meant to be used for dictionaries of trees/splits indexed by flu segments.
"""
function resolve(trees::Dict{T, <:Tree}, splits::Dict{T, <:SplitList}; kwargs...) where T
	resolved_trees = Dict(k => copy(t) for (k,t) in trees)
	for (s,S) in splits
		resolve!(resolved_trees[s], S; kwargs...)
	end
	return resolved_trees
end


"""
	resolve!(S::SplitList, t::Tree, S2::SplitList)

Add splits of `Sk` in `S` list if they resolve `t` in list.
"""
function resolve!(Snew, S::Vector{SplitList{String}}, t::Vector{Tree{T}}, Sk::SplitList) where T
	c = 0
	for sk in Sk
		count = 0
		resolve_list = []
		for i in 1:length(t)
			ri = t[i].lleaves[Sk.leaves[sk.dat[1]]]
			for l in 2:length(sk.dat)
				ri = lca(ri, t[i].lleaves[Sk.leaves[sk.dat[l]]])
			end
			si = S[i].splitmap[ri.label]
			if si==sk
				count +=1
			elseif si != sk && !in(sk, Snew[i]) && arecompatible(si, sk)
				# Consider the set of splits just below ri that are subsplits of sk
				# If I join those, I should get exactly sk
				# Otherwise, can't use sk to resolve ri
				stmp = Split(0)
				for n in ri.child
					if n.isleaf
						l = findfirst(==(n.label), Sk.leaves)
						if in(l, sk.dat)
							TreeTools.joinsplits!(stmp, Split([l]))
						end
					else
						if TreeTools.is_sub_split(S[i].splitmap[n.label], sk)
							TreeTools.joinsplits!(stmp, S[i].splitmap[n.label])
						end
					end
				end

				if stmp == sk
					count +=1
					push!(resolve_list, i)
				else
					break
				end
			end
			## only resolve if compatible with all other trees
			if count==length(t)
				for i in resolve_list
					push!(S[i].splits, sk)
					push!(Snew[i].splits, sk)
					c += 1
				end
			end
		end
	end

	return c
end

function resolve!(S::Vector{SplitList{String}}, t::Vector{Tree{T}}) where T
	nit = 0
	nitmax = 20
	flag = true
	Snew = [SplitList(s.leaves) for s in S]
	while flag && nit < nitmax
		flag = false
		for k in 1:length(S)
			c = resolve!(Snew[1:end .!= k], S[1:end .!= k], t[1:end .!= k], S[k])
			c != 0 && (flag = true)
		end
		nit += 1
	end

	if nit == nitmax
		@warn "Maximum number of iterations reached"
	end

	return Snew
end

###############################################################################################################
################################### Resolve trees using eachother splits ######################################
###############################################################################################################

"""
	resolve!(t1::Tree, t2::Tree, tn::Vararg{Tree}; tau=0., shared_maps=nothing)

Resolve `t1` using splits of `t2` and inversely. Every split of `t2` a tree that is
compatible with `t1` is introduced in `t1` with branch length `tau` (and inversely). 
Return new splits in each tree.
If more than two trees are given, only introduce a split from `ti` in tree `tj` if 
the split is compatible with all trees.
Optionally add a `shared_maps` dictionary (`Dict{String, Bool}`) with the tree 
node label names as keys and if the corresponding branch is shared or not as value. 
The resolve function will update this dictionary with any new nodes/branches.
"""
function resolve!(t1::Tree, t2::Tree, tn::Vararg{Tree}; tau=0., shared_maps=nothing)
	S = [SplitList(t) for t in (t1,t2, tn...)]
	Snew = resolve!(S, [t1, t2, tn...])
	if !isnothing(shared_maps) && length(S) == 2
		for (t, s, d) in zip((t1,t2), S, shared_maps)
			resolve!(t, s; conflict=:fail, usemask=false, tau, shared_map=d)
		end
	else
		for (t, s) in zip((t1,t2,tn...), S)
			resolve!(t, s; conflict=:fail, usemask=false, tau)
		end
	end

	return Snew
end



###############################################################################################################
############################## Resolve trees using inferred compatible clades #################################
###############################################################################################################

"""
	resolve!(t1::Tree{T}, t2::Tree{T}, MCCs; tau = 0.)

Resolve `t1` using `t2` and inversely using the list of MCCs.
New branches have a length `tau`.
Return the list of resolved splits in each tree.
"""

function resolve!(t1::Tree{T}, t2::Tree{T}, MCCs; tau = 0., strict=false) where T 
	resolvable_splits = TreeKnit.new_splits(MCCs, t1, t2; strict)
	resolve!(t1, resolvable_splits[1]; conflict=:fail, tau)
	resolve!(t2, resolvable_splits[2]; conflict=:fail, tau)

	return resolvable_splits
end


