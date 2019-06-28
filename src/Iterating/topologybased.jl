function prunetrees_sa(segtrees ; maxsize = 50, verbose=true)
	st = deepcopy(segtrees)

	# Initialization
	verbose && println("Resolving trees...")
	_resolve_trees!(st)
	MCC_init = maximal_coherent_clades(collect(values(st)))
	verbose && println("Found $(length(MCC_init)) MCCs of average size $(mean([length(x) for x in MCC_init]))")
	mcc_names = name_mcc_clades!(collect(values(st)), MCC_init)
	if length(MCC_init) == 1
		println("Found only one MCC, nothing to prune.")
		return segtrees, [MCC_init[1]]
	end
	supra = unique([Set(x) for x in RecombTools.supraMCC(values(st), MCC_init)])
	verbose && println("Found $(length(supra)) common clades from MCCs")

	#
	mcclist = []
	for s in supra
		subtrees = make_subtrees(st, s)
		locMCC = maximal_coherent_clades(collect(values(subtrees)))
		if length(locMCC) < maxsize
			out = SplitGraph.opttrees(collect(values(subtrees))..., γ=0.95, Trange = .5:-0.01:0.05, M=2500);
			if out[end]
				append!(mcclist, unique(cat(out[1]..., dims=1)))
			else
				verbose && @warn "Simulated annealing did not converge. Mismatches remain."
			end
		end
	end

	return st, unique(mcclist)
end


"""
"""
function make_subtrees(segtrees, s)
	subtrees = Dict{String, Tree}()
	for (k,t) in segtrees
	    r = lca([t.lnodes[x] for x in s])
	    if !r.isroot
	        subtrees[k] = node2tree(prunenode(r)[1])
	    else
	        subtrees[k] = deepcopy(segtrees[k])
	    end
	end	
	return subtrees
end

