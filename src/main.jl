export computeMCCs!

"""
	computeMCCs(trees::Dict{<:Any, <:Tree}, oa::OptArgs=OptArgs(); preresolve=true, naive=false)
"""
function computeMCCs(trees::Dict{<:Any, <:Tree}, oa::OptArgs=OptArgs(); preresolve=true, naive=false)
	ct = deepcopy(trees)
	computeMCCs!(ct, oa, preresolve=preresolve, naive=naive)
end

"""
	computeMCCs!(trees::Dict{<:Any, <:Tree}, oa::OptArgs; preresolve=true, naive=false)

Compute pairwise MCCs for `trees` by calling function `SplitGraph.runopt(oa,t1,t2)` on pairs of trees. Return MCCs and resolved splits. 
About `preresolve`: 
- If `true`, a first pass of MCC computation is made and trees are resolved using the results, keeping only compatible splits. A second pass of MCC computation is then made without resolving. 
- Else, only the first pass is performed. For more than two trees, this may find MCCs that introduce incompatible splits in case of poorly resolved input trees. 
"""
function computeMCCs!(trees::Dict{<:Any, <:Tree}, oa::OptArgs=OptArgs(); preresolve=true, naive=false)
	if naive
		return computeMCCs_naive!(trees)
	end
	if preresolve
		return computeMCCs_preresolve!(trees, oa)
	else
		return computeMCCs_dynresolve!(trees, oa)
	end
end

function computeMCCs_naive!(trees::Dict{<:Any, <:Tree})
	rS = resolve!(values(trees)...)
	resolved_splits = Dict(k=>rS[i] for (i,k) in enumerate(keys(trees)))
	MCCs = _computeMCCs(ts -> RecombTools.maximal_coherent_clades(collect(values(ts))...), trees)
	return MCCs, resolved_splits	
end

function computeMCCs_preresolve!(trees::Dict{<:Any, <:Tree}, oa::OptArgs)
	oac = @set oa.resolve = true
	oac = @set oa.output = :mccs
	resolved_splits = resolve_from_mccs!(ts -> runopt(oac,ts), trees)
	oac = @set oa.resolve = false
	MCCs = _computeMCCs(ts -> runopt(oac,ts), trees)
	return MCCs, resolved_splits
end

function computeMCCs_dynresolve!(trees::Dict{<:Any, <:Tree}, oa::OptArgs)
	MCCs = _computeMCCs(ts -> runopt(oa,ts), trees)
	resolved_splits = new_splits(trees, MCCs)
	resolved_splits = Dict(k=>TreeTools.cat(aS...) for (k,aS) in resolved_splits) # new_splits returns a Dict{T, Array{SplitList}}
	for k in keys(trees)
		resolve!(trees[k], resolved_splits[k])
	end
	return MCCs, resolved_splits
end

"""
	_computeMCCs(f::Function, trees::Dict{<:Any, <:Tree})

Compute MCCs for each pair of trees in `trees` by calling `f(t1,t2)`. Return a dictionary indexed with pairs of keys of `trees`. 
For simplicity, output for a pair of the same key is `[String[]]` (instead of not being indexed at all)
"""
function _computeMCCs(f::Function, trees::Dict{<:Any, <:Tree})
	MCCs = Dict()
	segments = collect(keys(trees))
	for i in 1:length(trees), j in (i+1):length(trees)
		s1 = segments[i]
		s2 = segments[j]
		MCCs[s1,s2] = f(Dict(s1=>trees[s1], s2=>trees[s2]))
		MCCs[s2,s1] = MCCs[s1,s2]
	end
	for s in segments
		MCCs[s,s] = [String[]]
	end
	return MCCs
end