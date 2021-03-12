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
function resolve_from_mccs!(infer_mccs::Function, trees::Dict{<:Any, Tree}; verbose=false, kwargs...)
	verbose && println("--- Resolving with MCCs... ---\n")
	f(t) = inter_mccs(t; kwargs...)
	# 1. `MCCs[u,v]` gives corresponds to `trees[u]` and `trees[v]` 
	MCCs = computeMCCs(f, trees)
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
		TreeTools.resolve!(trees[s], S, conflict=:fail)
		TreeTools.check_tree(trees[s])
	end
	# 
	return cmpt_splits
end

"""
	resolve_from_mccs!(infer_mccs::Function, trees::Vararg{Tree}; verbose=false, kwargs...)
"""
resolve_from_mccs!(infer_mccs::Function, trees::Vararg{Tree}; verbose=false, kwargs...) = resolve_from_mccs!(infer_mccs, Dict(i=>t for (i,t) in enumerate(trees)), verbose=verbose; kwargs...)