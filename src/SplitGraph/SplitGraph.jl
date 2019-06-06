module SplitGraph

using TreeTools
using RecombTools

export opttrees 

include("objects.jl")
include("tools.jl")
include("energy.jl")


"""
"""
opttrees(t... ; γ=1.05, Trange=0.5:-0.01:0.05, M = 1000) = opttrees(γ, Trange, M, t...)
function opttrees(γ, Trange, M, t::Vararg{Tree})
	treelist = collect(t)
	mcc = maximal_coherent_clades(treelist)
	mcc_names = name_mcc_clades!(treelist, mcc)
	for (i,t) in enumerate(treelist)
		treelist[i] = reduce_to_mcc(t, mcc)
	end
	g = trees2graph(treelist)
	oconf, E, F = sa_opt(g, γ=γ, Trange=Trange, M=M)
	converged = compute_energy(oconf,g)==0
	return [mcc_names[x] for x in g.labels[.!oconf]], g.labels[.!oconf],E,F,converged
end

end