###############################################################################################################
########################################## Basic resolve functions ############################################
###############################################################################################################

"""
	resolve!(t::Tree, S::SplitList; conflict=:fail, usemask=false, tau=0.)

Add splits in `S` to `t` by introducing internal nodes. New nodes are assigned a time `tau` (`0` by default).
If `conflict != :ignore`, will fail if a split `s` in `S` is not compatible with `t`. Otherwise, silently skip the conflicting splits.
"""
function resolve!(
	t::Tree{T}, S::SplitList;
	conflict=:fail, usemask=false, tau=0., safe=false
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
					prunenode!(r)
					graftnode!(nr,r)
				end
				graftnode!(R, nr, tau=tau)
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

Resolve `trees[s]` with splits in `splits[s]` by calling `resolve!`. Return dictionary of resolved trees. `trees` and `splits` must share keys. This is meant to be used for dictionaries of trees/splits indexed by flu segments.
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
					c += 1
				else
					break
				end
			end
			if count==length(t)
				for i in resolve_list
					push!(S[i].splits, sk)
					push!(Snew[i].splits, sk)
				end
			end
		end
	end

	return c
end

function resolve!(S::Vector{SplitList{String}}, t::Vector{Tree{T}}) where T
	print("new resolve function")
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
	resolve!(t1::Tree, t2::Tree; tau=0.)

Resolve `t1` using splits of `t2` and inversely. Every split of `t2` a tree that is compatible with `t1` is introduced in `t1` with branch length `tau` (and inversely). Return new splits in each tree.
"""
function resolve!(t1::Tree, t2::Tree, tn::Vararg{Tree}; tau=0.)
	S = [SplitList(t) for t in (t1,t2, tn...)]
	Snew = resolve!(S, [t1, t2, tn...])
	for (t, s) in zip((t1,t2,tn...), S)
		resolve!(t, s, conflict=:fail, usemask=false, tau=tau)
	end

	return Snew
end




###############################################################################################################
############################## Resolve trees using inferred compatible clades #################################
###############################################################################################################

"""
	resolve!(t1::Tree, t2::Tree, MCCs; tau = 0.)

Resolve `t1` using `t2` and inversely using the list of MCCs.
New branches have a length `tau`.
Return the list of resolved splits in each tree.
"""
function resolve!(t1::Tree, t2::Tree, MCCs; tau = 0.)
	resolvable_splits = TreeKnit.new_splits(MCCs, t1, t2)
	resolve!(t1, resolvable_splits[1]; conflict=:fail, tau)
	resolve!(t2, resolvable_splits[2]; conflict=:fail, tau)

	return resolvable_splits
end




###############################################################################################################
########################### Resolve multiple trees in a list in a compatible way ##############################
###############################################################################################################



###############################################################################################################
############################## Resolve trees in a compatible way #################################
###############################################################################################################

## This was not working great and it's removed from current TreeKnit. Keeping it here for reference.

# """
# 	resolve_from_mccs!(infer_mccs::Function, trees::Vararg{Tree}; verbose=false, kwargs...)
# 	resolve_from_mccs!(infer_mccs::Function, trees::Dict; verbose=false, kwargs...)

# Resolve `trees` using pairwise MCC **inference**.
# The `infer_mccs` function must take a `Dict{<:Any, Tree}` as input.
# `kwargs...` are fed to `infer_mccs`.
# Return the set of compatible splits introduced.

# This is done in four steps:
# 1. Compute MCCs for pairs of trees while resolving them
# 2. Find all new splits in each tree that are introduced by the found MCCs
# 3. Select only the splits compatible with all others
# 4. Introduce those in trees

# These four steps are **iterated**, while new compatible splits are found.
# When no new compatible split is found, compute final MCCs without resolving.
# """
# function resolve_from_mccs!(infer_mccs::Function, trees::Dict{<:Any, <:Tree}; verbose=false)
# 	maxit = 10
# 	#
# 	verbose && println("It. 1/$(maxit)")
# 	newsplits = _resolve_from_mccs!(infer_mccs, trees, verbose=verbose)
# 	#
# 	nS = newsplits
# 	it = 1
# 	while !mapreduce(isempty, *, values(nS), init=true) && it < maxit
# 		verbose && println("It. $(it+1)/$(maxit)")
# 		nS = _resolve_from_mccs!(infer_mccs, trees, verbose=verbose)
# 		for (s,S) in nS
# 			append!(newsplits[s].splits, S.splits)
# 		end
# 		it += 1
# 	end
# 	return newsplits
# end

# function _resolve_from_mccs!(infer_mccs::Function, trees::Dict{<:Any, <:Tree}; verbose=false)
# 	verbose && println("--- Resolving with MCCs... ---\n")
# 	# 1. `MCCs[u,v]` corresponds to `trees[u]` and `trees[v]`
# 	MCCs = _computeMCCs(infer_mccs, trees)
# 	# 2. `resolvable_splits[s]` is an array of `SplitList` objects that can be introduced in `trees[s]` from each other tree.
# 	resolvable_splits = TreeKnit.new_splits(trees, MCCs)
# 	if verbose
# 		for (s,t) in trees
# 			Y = union(resolvable_splits[s]...)
# 			println("Tree $s: $(length(Y)) resolvable splits (potentially incompatible).")
# 		end
# 		println()
# 	end
# 	# 3. `cmpt_splits[s]` is a `SplitList` object containing a selection of splits
# 	# we can choose between several methods for this choice.
# 	# cmpt_splits = TreeKnit.max_clique_splits(resolvable_splits; verbose)
# 	cmpt_splits = TreeKnit.compatible_splits(resolvable_splits)
# 	if verbose
# 		for (s,t) in trees
# 			println("Tree $s: $(length(cmpt_splits[s])) compatible resolvable splits.")
# 		end
# 		println()
# 	end
# 	# 4.
# 	for (s,S) in cmpt_splits
# 		resolve!(trees[s], S, conflict=:fail)
# 		TreeTools.check_tree(trees[s])
# 	end
# 	#
# 	return cmpt_splits
# end


# function resolve_from_mccs!(
# 	infer_mccs::Function, trees::Vararg{Tree};
# 	verbose=false, kwargs...
# )
# 	resolve_from_mccs!(
# 		infer_mccs, Dict(i=>t for (i,t) in enumerate(trees));
# 		verbose=verbose, kwargs...
# 	)
# end


# """
#     compatible_splits(rS::Dict{<:Any, Array{SplitList,1}})

# Find splits in `rS` that are compatible with all other splits. Return a dictionary with the same keys as `rS`.
# """
# function compatible_splits(rS::Dict{<:Any, Array{SplitList,1}})
#     cs = Dict{Any, SplitList}() # Compatible splits (output)
#     for (seg, aS) in rS # aS: new splits for segment seg from all other trees
#         for (i,S) in enumerate(aS) # S: new splits in seg from tree i
#             if !haskey(cs, seg)
#                 cs[seg] = SplitList(S.leaves)
#             end
#             for s in S
#                 # Check if s is compatible with aS[j] for all `j!=i`
#                 cmpt = true
#                 for j in Iterators.filter(!=(i), 1:length(aS))
#                     if !iscompatible(s, aS[j], usemask=false)
#                         cmpt = false
#                         break
#                     end
#                 end
#                 cmpt && push!(cs[seg].splits, s)
#             end
#         end
#         cs[seg] = unique(cs[seg])
#     end
#     return cs
# end

# function max_clique_splits(nS::Dict; verbose=false)
# 	cS = Dict{Any, SplitList}()
# 	for (seg, S) in nS
# 		verbose && println(seg)
# 		cS[seg] = max_clique_splits(S; verbose)
# 	end
# 	return cS
# end
# function max_clique_splits(nS; verbose=false)
# 	S = union(nS...)
# 	verbose && println("Finding max clique among $(length(S)) splits.")
# 	g = build_compat_graph(S)
# 	all_cliques = LightGraphs.maximal_cliques(g)
# 	verbose && println("Found $(length(all_cliques)) cliques.")
# 	if isempty(all_cliques)
# 		return S
# 	else
# 		max_clique = sort(all_cliques, by=x->length(x), rev=true)[1]
# 		# Delete splits not in the max clique
# 		todel = findall(!in(max_clique), 1:length(S))
# 		deleteat!(S.splits, todel)
# 		return S
# 	end
# end

# function build_compat_graph(S::SplitList)
# 	L = length(S)
# 	g = Graph(L)
# 	for i in 1:L, j in (i+1):L
# 		if arecompatible(S[i], S[j])
# 			add_edge!(g, i, j)
# 		end
# 	end
# 	return g
# end
