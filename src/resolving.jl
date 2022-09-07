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
				if usemask
					roots = TreeTools.blca([t.lleaves[x] for x in leaves(S,i)]...)
				else
					roots = TreeTools.blca([t.lleaves[x] for x in leaves(S,i)]...)
				end
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
	resolve!(S1::SplitList, t1::Tree, S2::SplitList)

Add splits of `S2` in `S1` if they resolve `t1`.
"""
function resolve!(S1new, S1::SplitList, t1::Tree, S2::SplitList)
	c = 0
	for s2 in S2
		r1 = t1.lleaves[S2.leaves[s2.dat[1]]]
		for i in 2:length(s2.dat)
			r1 = lca(r1, t1.lleaves[S2.leaves[s2.dat[i]]])
		end
		#r1 = lca(t1, S2.leaves[s2.dat]) # Ancestor of nodes in s2 in t1
		s1 = S1.splitmap[r1.label]
		if s1 != s2 && !in(s2, S1new) && arecompatible(s1, s2)
			# Consider the set of splits just below r1 that are subsplits of s2
			# If I join those, I should get exactly s2
			# Otherwise, can't use s2 to resolve r1
			stmp = Split(0)
			for n in r1.child
				if n.isleaf
					i = findfirst(==(n.label), S2.leaves)
					if in(i, s2.dat)
						TreeTools.joinsplits!(stmp, Split([i]))
					end
				else
					if TreeTools.is_sub_split(S1.splitmap[n.label], s2)
						TreeTools.joinsplits!(stmp, S1.splitmap[n.label])
					end
				end
			end

			if stmp == s2
				push!(S1.splits, s2)
				push!(S1new.splits, s2)
				c += 1
			end
		end
	end

	return c
end

function resolve!(S1::SplitList, S2::SplitList, t1::Tree, t2::Tree)
	nit = 0
	nitmax = 20
	flag = true
	S1new = SplitList(S1.leaves)
	S2new = SplitList(S2.leaves)
	while flag && nit < nitmax
		flag = false
		c = resolve!(S1new, S1, t1, S2)
		c != 0 && (flag = true)
		c = resolve!(S2new, S2, t2, S1)
		c != 0 && (flag = true)
		nit += 1
	end

	if nit == nitmax
		@warn "Maximum number of iterations reached"
	end

	return [S1new, S2new]
end

###############################################################################################################
################################### Resolve trees using eachother splits ######################################
###############################################################################################################

"""
	resolve!(t1::Tree, t2::Tree; tau=0.)

Resolve `t1` using splits of `t2` and inversely. Every split of `t2` a tree that is compatible with `t1` is introduced in `t1` with branch length `tau` (and inversely). Return new splits in each tree.
"""
function resolve!(t1::Tree, t2::Tree; tau=0.)
	S = [SplitList(t) for t in (t1,t2)]
	Snew = resolve!(S[1], S[2], t1, t2)
	for (t, s) in zip((t1,t2), S)
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
