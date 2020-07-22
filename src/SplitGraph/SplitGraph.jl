module SplitGraph

using TreeTools
using RecombTools

export opttrees!

include("objects.jl")
include("tools.jl")
include("energy.jl")


"""
	opttrees!(t... ; γ=1.05, Trange=0.5:-0.01:0.05, M = 1000)

Return a list of MCCs for input trees. 
Output:
1. 
"""
opttrees!(t... ; γ=1.05, Trange=0.5:-0.01:0.05, M = 1000) = opttrees!(γ, Trange, M, t...)
function opttrees!(γ, Trange, M, t::Vararg{Tree})
	treelist = collect(t)
	mcc = maximal_coherent_clades(treelist)
	mcc_names = name_mcc_clades!(treelist, mcc)
	for (i,t) in enumerate(treelist)
		treelist[i] = reduce_to_mcc(t, mcc)
	end
	g = trees2graph(treelist)
	oconf, E, F = sa_opt(g, γ=γ, Trange=Trange, M=M)
	converged = compute_energy(oconf[1],g)==0
	# Note
	# If !converged, then oconf does not represent an overarching MCC, but rather a collection of MCCs that do not include the root node
	# In this case, my MCCs should be 
	# 1. [mcc_names[x] for x in g.labels[.!conf]], that is all the MCCs that are *not* in oconf
	# 2. The MCCs that I find by keeping only leaves that are in oconf. This could very well be *less* than sum(oconf), since some mismatches may have been removed. 
	# This means that I should start at high γ, find an oconf such that sum(oconf)!=L, and this way I can confidently remove mccs corresponding to 1.
	# I can then iterate this process, *updating* my definition of MCCs as I go. 
	# This requires a meta-optmization that recomputes MCCs every time.  
	return [[mcc_names[x] for x in g.labels[.!conf]] for conf in oconf], [g.labels[.!conf] for conf in oconf], E, F, converged
end

"""
	pruneconf!(clades, trees::Vararg{Tree})

Prune `clades` from `trees...`. 
"""
function pruneconf!(clades, trees::Vararg{Tree})
	for t in trees
		for st in clades
			TreeTools.prunesubtree!(t, st, clade_only=true)
		end
	end
end
pruneconf!(trees, clades) = pruneconf!(clades, trees...)
"""
	pruneconf!(trees, mcc_names, mcc_conf)

Prune MCCs `mcc_names[x]` for all `x` in `mcc_conf` from trees `t...`. 
"""
pruneconf!(trees, mcc_names, mcc_conf) = pruneconf!([mcc_names[x] for x in mcc_conf], trees...)


end