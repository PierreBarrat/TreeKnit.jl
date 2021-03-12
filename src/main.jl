using SplitGraph

"""
	computeMCCs(trees::Dict{<:Any, Tree}, oa::OptArgs; preresolve=true)

Compute pairwise MCCs for `trees` by calling function `SplitGraph.runopt(oa,t1,t2)` on pairs of trees. Return MCCs and resolved splits. 
About `preresolve`: 
- If `true`, a first pass of MCC computation is made and trees are resolved using the results, keeping only compatible splits. A second pass of MCC computation is then made without resolving. 
- Else, only the first pass is performed. For more than two trees, this may find MCCs that introduce incompatible splits in case of poorly resolved input trees. 
"""
function computeMCCs(trees::Dict{<:Any, Tree}, oa::OptArgs; preresolve=true)
	if resolve
		return computeMCCs_preresolve(trees, oa)
	else
		return computeMCCs_dynresolve(trees, oa)
	end
end

function computeMCCs_preresolve(trees::Dict{<:Any, Tree}, oa::OptArgs)
	oac = @set oa.resolve = true
	resolved_splits = resolve_from_mccs!((t1,t2) -> runopt(oac,t1,t2), trees)
	oac = @set oa.resolve = false
	MCCs = _computeMCCs((t1,t2) -> runopt(oac,t1,t2), trees)
	return MCCs, resolved_splits
end

function computeMCCs_dynresolve(trees::Dict{<:Any, Tree}, oa::OptArgs)
	MCCs = _computeMCCs((t1,t2) -> runopt(oa,t1,t2), trees)
	resolved_splits = newsplits(trees, MCCs)
	return MCCs, resolved_splits
end

"""
	_computeMCCs(f::Function, trees::Dict{<:Any, Tree})

Compute MCCs for each pair of trees in `trees` by calling `f(t1,t2)`. Return a dictionary indexed with pairs of keys of `trees`. 
For simplicity, output for a pair of the same key is `[String[]]` (instead of not being indexed at all)
"""
function _computeMCCs(f::Function, trees::Dict{<:Any, Tree})
	MCCs = Dict()
	segments = sort(collect(keys(trees)))
	for i in 1:length(trees), j in (i+1):length(trees)
		s1 = segments[i]
		s2 = segments[j]
		MCCs[s1,s2] = f(trees[s1], trees[s2])
		MCCs[s2,s1] = MCCs[s1,s2]
	end
	for s in segments
		MCCs[s,s] = [String[]]
	end
	return MCCs
end