###############################################################################################################
########################################## Basic resolve functions ############################################
###############################################################################################################

"""
	resolve!(t::Tree, S::SplitList; conflict=:fail, usemask=false, tau=0.)

Add splits in `S` to `t` by introducing internal nodes. New nodes are assigned a time `tau` (`0` by default).
If `conflict != :ignore`, will fail if a split `s` in `S` is not compatible with `t`. Otherwise, silently skip the conflicting splits.
"""
function resolve!(t::Tree{T}, S::SplitList; conflict=:fail, usemask=false, tau=0.) where T
	# Label for created nodes
	label_i = parse(Int64, TreeTools.create_label(t, "RESOLVED")[10:end])
	#
	tsplits = SplitList(t)
	for (i,s) in enumerate(S)
		if !in(s, tsplits; usemask)
			if iscompatible(s, tsplits; usemask)
				if usemask
					roots = TreeTools.blca([t.lleaves[x] for x in S.leaves[S[i].dat .* S.mask]]...)
				else
					roots = TreeTools.blca([t.lleaves[x] for x in S.leaves[S[i].dat]]...)
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
				@error "Conflicting splits"
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
	resolved_trees = deepcopy(trees)
	for (s,S) in splits
		resolve!(resolved_trees[s], S; kwargs...)
	end
	return resolved_trees
end

#"""
#	resolve_ref!(Sref::SplitList, S::Vararg{SplitList}, usemask=false)

#Add new and compatible splits of `S` to `Sref`. If `usemask`, masks are used to determine compatibility. Return the number of added splits.
#**Note**: the order of `S` matters if its elements contain incompatible splits!
#"""
#function resolve_ref!(Sref::SplitList, S::Vararg{SplitList}; usemask=false)
#	c = 0
#	for s in S
#		for x in s
#			if !in(x, Sref; usemask) && iscompatible(x, Sref; usemask)
#				push!(Sref.splits, x)
#				c += 1
#			end
#		end
#	end
#	return c
#end

#"""
#	resolve!(S::Vararg{SplitList})

#Resolve each element of `S` using other elements by calling `resolve!(S[i],S)` for all `i` until no new split can be added. If `usemask`, masks are used to determine compatibility.
#"""
#function resolve!(S::Vararg{SplitList}; usemask=false)
#	nit = 0
#	nitmax = 20
#	flag = true
#	while flag && nit < nitmax
#		flag = false
#		for s in S
#			c = resolve_ref!(s, S...; usemask)
#			if c != 0
#				flag = true
#			end
#		end
#		nit += 1
#	end
#	if nit == nitmax
#		@warn "Maximum number of iterations reached"
#	end
#	nothing
#end


#####
#####
"""
	resolve!(S1::SplitList, t1::Tree, S2::SplitList)

Add splits of `S2` in `S1` if they resolve `t1`.
"""
function resolve!(S1new, S1::SplitList, t1::Tree, S2::SplitList)
	c = 0
	for s2 in S2
		r1 = lca(t1, S2.leaves[s2.dat]) # Ancestor of nodes in s2 in t1
		s1 = S1.splitmap[r1.label]
		if s1 != s2 && !in(s2, S1new) && arecompatible(s1, s2)
			# Consider the set of splits just below r1 that are subsplits of s2
			# If I join those, I should get exactly s2
			# Otherwise, can't use s2 to resolve r1
			stmp = Split(length(s2))
			for n in r1.child
				if n.isleaf
					i = findfirst(==(n.label), S2.leaves)
					if s2.dat[i]
						stmp.dat[i] = true
					end
				else
					if TreeTools.is_sub_split(S1.splitmap[n.label], s2)
						TreeTools.joinsplits!(stmp, S1.splitmap[n.label])
					end
				end
			end
			if stmp == s2
				println("add split $(S2.leaves[s2.dat])")
				push!(S1.splits, s2)
				push!(S1new, s2)
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
	S1new = []
	S2new = []
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

	return nothing
end

###############################################################################################################
################################### Resolve trees using eachother splits ######################################
###############################################################################################################

"""
	resolve!(t1::Tree, t2::Tree; tau=0.)

Resolve `t1` using splits of `t2` and inversely. Every split of `t2` a tree that is compatible with `t1` is introduced in `t1` (and inversely). Return new splits in each tree.
"""
function resolve!(t1::Tree, t2::Tree; tau=0.)
	S = [SplitList(t) for t in (t1,t2)]
	tsplits_a = deepcopy(S) # for comparison in return value
	#resolve!(S...; usemask=false)
	resolve!(S[1], S[2], t1, t2)
	for (t,s) in zip((t1,t2), S)
		resolve!(t, s, conflict=:fail, usemask=false, tau=tau)
	end

	return [setdiff(S[i], tsplits_a[i], :none) for i in eachindex(S)]
end



###############################################################################################################
############################## Resolve trees using inferred compatible clades #################################
###############################################################################################################

"""
	resolve_from_mccs!(infer_mccs::Function, trees::Dict; verbose=false, kwargs...)

Resolve `trees` using pairwise MCC inference. The `infer_mccs` function must take a `Dict{<:Any, Tree}` as input. `kwargs...` are fed to `infer_mccs`. Return the set of compatible splits introduced.

This is done in four steps:
1. Compute MCCs for pairs of trees while resolving them
2. Find all new splits in each tree that are introduced by the found MCCs
3. Select only the splits compatible with all others
4. Introduce those in trees

These four steps are **iterated**, while new compatible splits are found.
When no new compatible split is found, compute final MCCs without resolving.
"""
function resolve_from_mccs!(infer_mccs::Function, trees::Dict{<:Any, <:Tree}; verbose=false)
	newsplits = _resolve_from_mccs!(infer_mccs, trees, verbose=verbose)
	nS = newsplits
	maxit = 10
	it = 0
	while !mapreduce(isempty, *, values(nS), init=true) && it < maxit
		verbose && println("It. $(it+1)/$(maxit)")
		nS = _resolve_from_mccs!(infer_mccs, trees, verbose=verbose)
		for (s,S) in nS
			append!(newsplits[s].splits, S.splits)
		end
		it += 1
	end
	return newsplits
end

function _resolve_from_mccs!(infer_mccs::Function, trees::Dict{<:Any, <:Tree}; verbose=false)
	verbose && println("--- Resolving with MCCs... ---\n")
	# 1. `MCCs[u,v]` gives corresponds to `trees[u]` and `trees[v]`
	MCCs = _computeMCCs(infer_mccs, trees)
	# 2. `resolvable_splits[s]` is an array of `SplitList` objects that can be introduced in `trees[s]` from each other tree.
	resolvable_splits = RecombTools.new_splits(trees, MCCs)
	if verbose
		for (s,t) in trees
			Y = unique(vcat([x.splits for x in resolvable_splits[s]]...))
			println("Tree $s: $(length(Y)) resolvable splits (potentially incompatible).")
		end
		println()
	end
	# 3. `cmpt_splits[s]` is a `SplitList` object containing only splits compatible with all others
	cmpt_splits = RecombTools.compatible_splits(resolvable_splits)
	if verbose
		for (s,t) in trees
			println("Tree $s: $(length(cmpt_splits[s])) compatible resolvable splits.")
		end
		println()
	end
	# 4.
	for (s,S) in cmpt_splits
		resolve!(trees[s], S, conflict=:fail)
		TreeTools.check_tree(trees[s])
	end
	#
	return cmpt_splits
end

"""
	resolve_from_mccs!(infer_mccs::Function, trees::Vararg{Tree}; verbose=false, kwargs...)
"""
resolve_from_mccs!(infer_mccs::Function, trees::Vararg{Tree}; verbose=false, kwargs...) = resolve_from_mccs!(infer_mccs, Dict(i=>t for (i,t) in enumerate(trees)), verbose=verbose; kwargs...)


"""
    compatible_splits(rS::Dict{<:Any, Array{SplitList,1}})

Find splits in `rS` that are compatible with all other splits. Return a dictionary with the same keys as `rS`.
"""
function compatible_splits(rS::Dict{<:Any, Array{SplitList,1}})
    cs = Dict{Any, SplitList}() # Compatible splits (output)
    for (seg, aS) in rS # aS: new splits for segment seg from all other trees
        for (i,S) in enumerate(aS) # S: new splits in seg from tree i
            if !haskey(cs, seg)
                cs[seg] = SplitList(S.leaves)
            end
            for s in S
                # Check if s is compatible with aS[j] for all `j!=i`
                cmpt = true
                for j in Iterators.filter(!=(i), 1:length(aS))
                    if !iscompatible(s, aS[j], usemask=false)
                        cmpt = false
                        break
                    end
                end
                cmpt && push!(cs[seg].splits, s)
            end
        end
        cs[seg] = unique(cs[seg])
    end
    return cs
end



###############################################################################################################
############################# Resolve polytomies using cross-mapped mutations #################################
###############################################################################################################
"""
	resolve_crossmapped_muts!(trees::Dict; cmkey=:cmmuts)
"""
function resolve_crossmapped_muts!(trees::Dict; cmkey=:cmmuts, infer_ancestral=true)
	infer_ancestral && ancestral_sequences!(trees, self=false, crossmapped=true)
	for (sref, tref) in trees
		Snew = resolve_crossmapped_muts!(tref, cmkey=cmkey)
	end
end
"""
	resolve_crossmapped_muts(t::Tree, refseg=nothing; cmkey=:cmmuts)

Resolve polytomies in `t` using crossmapped mutations from other segments. For each tree node `n`, if `isnothing(refseg)`, consider all other segments in `n.data.dat[cmkey]`. Otherwise, consider only mutations in `n.data.dat[cmkey][s]` for `s` in `refseg`.
"""
function resolve_crossmapped_muts(t::Tree, refseg=nothing; cmkey=:cmmuts)
	Sref = SplitList(t)
	Snew = SplitList(copy(Sref.leaves))
	for n in values(t.lnodes)
		if length(n.child) > 2
			polyS = resolve_polytomy(n, refseg; cmkey=cmkey)
			_add_splits!(Snew, polyS, Sref)
		end
	end
	return Snew
end
resolve_crossmapped_muts!(t::Tree, refseg=nothing; cmkey=:cmmuts) = resolve!(t, resolve_crossmapped_muts(t, refseg, cmkey=cmkey))

function _add_splits!(Snew, polyS, Sref)
	for ps in polyS
		snew = TreeTools.Split(length(Snew.leaves))
		for (i,x) in enumerate(ps.dat)
			if x # Polytomy-level leaf `i` is in the split
				sref = get(Sref.splitmap, polyS.leaves[i]) do
					snew.dat[findfirst(==(polyS.leaves[i]), Sref.leaves)] = true
					nothing
				end # Corresponding split at the tree level - if it's a leaf, it won't be there, hence the get(f,dict,key) do syntax
				!isnothing(sref) && TreeTools.joinsplits!(snew, sref)
			end
		end
		push!(Snew.splits, snew)
	end
end

"""
	resolve_polytomy(a::TreeNode, refseg=nothing; cmkey=:cmmuts)

Attempt to resolve polytomy with root `a` using mutations from other segments. For all children of `c` of `a`, search for mutations in other segments `s` in `c.data.dat[cmkey][s]`. If `isnothing(refseg)`, consider all other segments. Otherwise, consider the ones in `refseg`.
Return a `TreeTools.SplitList` object containing new splits. Does not modify the tree.

## Warning
Returned splits are based on nodes below `a`.
"""
function resolve_polytomy(a::TreeNode, refseg=nothing; cmkey=:cmmuts)
	# Count mutations
	mutcount = Dict{Any, Int64}()
	for c in a.child
		for (s,muts) in c.data.dat[cmkey], m in muts
			if isnothing(refseg) || _inrefseg(s,refseg)
				mutcount[s,m] = get(mutcount, (s,m), 0) + 1
			end
		end
	end
	# Build character table from mutations that appear more than once
	chartab, labels, charid = get_cmmut_chartab(a, mutcount, cmkey)
	#
	return CompatibilityTree.max_compatibility_splits(chartab, labels)
end
_inrefseg(s, refseg::AbstractArray) = in(s,refseg)
_inrefseg(s, refseg::AbstractString) = (s == refseg)

function get_cmmut_chartab(a::TreeNode, mutcount, cmkey)
	# Map of mutation to id in chartable
	charid = Dict()
	j = 1
	for (m,c) in mutcount
		if c > 1
			charid[m] = j
			j += 1
		end
	end
	nsegregating = length(charid)
	#
	chartab = zeros(Bool, length(a.child) + 1, nsegregating)
	labels = Array{String,1}(undef, length(a.child) + 1)
	for (i,c) in enumerate(a.child)
		labels[i] = c.label
		for (s,muts) in c.data.dat[cmkey], m in muts
			if haskey(charid, (s,m))
				chartab[i,charid[s,m]] = true
			end
		end
	end
	labels[end] = a.label
	return chartab, labels, charid
end
