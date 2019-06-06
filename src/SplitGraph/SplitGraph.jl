module SplitGraph

using TreeTools
using RecombTools

export prunetrees 

include("objects.jl")
include("tools.jl")
include("energy.jl")


"""
"""
prunetrees(t... ; γ=1.05, Trange=0.5:-0.01:0.05, M = 1000) = prunetrees(γ, Trange, M, t...)
function prunetrees(γ, Trange, M, t::Vararg{Tree})
	treelist = collect(t)
	mcc = maximal_coherent_clades(treelist)
	mcc_names = name_mcc_clades!(treelist, mcc)
	for (i,t) in enumerate(treelist)
		treelist[i] = reduce_to_mcc(t, mcc)
	end
	g = trees2graph(treelist)
	oconf, E, F = sa_opt(g, γ=γ, Trange=Trange, M=M)
	return [mcc_names[x] for x in g.labels[.!oconf]],E,F
end

end