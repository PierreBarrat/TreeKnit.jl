module SplitGraph

using TreeTools
using RecombTools
using StatsBase
using DataFrames
using SpecialFunctions
using Parameters


export runopt, OptArgs

include("objects.jl")
include("tools.jl")
include("energy.jl")
include("likelihood.jl")

let verbose::Bool = false, vverbose::Bool = false
	global v() = verbose
	global vv() = vverbose
	global set_verbose(v) = (verbose = v)
	global set_vverbose(v) = (vverbose = v)
end

"""
		runopt(t::Vararg{Tree}; kwargs...)
		runopt(oa::OptArgs, t::Vararg{Tree})

Run optimization at constant γ. See `?Optargs` for arguments. 
"""
runopt(t::Vararg{Tree}; kwargs...) = runopt(OptArgs(;kwargs...), t...)

function runopt(oa::OptArgs, t::Vararg{Tree})
	# 
	ot = deepcopy(collect(t))
	resolve_trees!(ot...)
	# 
	iMCCs = maximal_coherent_clades(ot)
	Einit = count_mismatches(ot...)
	n0 = length(first(ot).lleaves)
	df = DataFrame(nleaves=Int64[n0],
		nMCCs=length(iMCCs),
		γ=Any[oa.γ], M=Any[missing], 
		Efinal=Any[Einit], Ffinal=Any[Einit],
		removedMCCs=Any[missing], all_removedMCCs=Any[missing], remainingMCCs=Any[iMCCs])
	MCCs = []
	Evals = Any[]
	Fvals = Any[]
	set_verbose(oa.verbose)
	set_vverbose(oa.vv)
	#
	for i in 1:oa.itmax
		flag = :init
		v() && println("\n --- \nIteration $i/$(oa.itmax) - $(df.nleaves[end]) leaves remaining")

		# Optimization
		n = length(first(ot).lleaves)
		M = getM(n, oa.Md)
		mccs, Efinal, Ffinal, E, F = opttrees!(ot..., γ=oa.γ, M=M, Trange=oa.Trange, likelihood_sort=oa.likelihood_sort)
		length(mccs) != 0 ? (flag = :found) : (v() && println("No solution found in current iteration"))
		append!(MCCs, mccs)

		# Checks
		!prod([check_tree(t) for t in ot]) && @error "Problem in one of the trees"
		
		# Actions if a final solution is found
		if flag == :found 
			if sum(length(m) for m in mccs) != length(ot[1].lleaves) # Found mccs do not cover all leaves
				pruneconf!(mccs, ot...)
				if complete_mccs!(MCCs, ot)
					v() && println("$flag: all mccs have been found")
					resolve_trees!(ot...)
					rMCCs = maximal_coherent_clades(ot)
					update_df!(df, length(ot[1].lleaves), length(rMCCs), oa.γ, M, Efinal, Ffinal, mccs, MCCs, rMCCs)
					break
				end
			else # Found mccs cover all leaves
				v() && println("$flag: all mccs have been found")
				resolve_trees!(ot...)
				rMCCs = maximal_coherent_clades(ot)
				update_df!(df, length(ot[1].lleaves), length(rMCCs), oa.γ, M, Efinal, Ffinal, mccs, MCCs, rMCCs)
				break
			end
		end

		# Action if a temporary solution is found (pruneconf! is called above)
		resolve_trees!(ot...)
		rMCCs = maximal_coherent_clades(ot)
		update_df!(df, length(ot[1].lleaves), length(rMCCs), oa.γ, M, Efinal, Ffinal, mccs, MCCs, rMCCs)
		push!(Evals, E)
		push!(Fvals, F)

		# Actions if no solution is found
		if i == oa.itmax || flag != :found
			v() && println("Maximum number of iterations reached or no solution found - stop")
			complete_mccs!(MCCs, ot, force=true)
			break
		end
	end
	return RecombTools.sort_mccs(MCCs), df, ot, Evals, Fvals
end


"""
	complete_mccs!(MCCs, ot; force=false)

Complete `MCCs` if `maximal_coherent_clades(ot)` is of length `1`, and return `true`. Otherwise (unless `force`), return `false`.
"""
function complete_mccs!(MCCs, ot; force=false)
	final_mccs = maximal_coherent_clades(ot)
	if length(final_mccs) != 1
		force && append!(MCCs, final_mccs)
		return false
	else
		append!(MCCs, final_mccs)
		return true
	end
end


"""
	opttrees!(t... ; γ=1.05, Trange=0.5:-0.01:0.05, M = 1000)

Return a list of MCCs for input trees. 
Output:
1. 
"""
opttrees!(t... ; γ=1.05, μ=ones(Float64, length(t)), Trange=reverse(0.01:0.05:1.1), M = 1000, likelihood_sort=true) = opttrees!(γ, Trange, M, μ, t...; likelihood_sort=likelihood_sort)
function opttrees!(γ, Trange, M, μ, t::Vararg{Tree}; likelihood_sort=true)
	# length(μ) != length(collect(t)) && @error "`μ` and `trees` do not have the same length."
	#
	treelist = collect(t)
	mcc = maximal_coherent_clades(treelist)
	if length(mcc) == 1
		return mcc, 0, 0., Int64[], Float64[]
	end
	mcc_names = name_mcc_clades!(treelist, mcc)
	for (i,t) in enumerate(treelist)
		treelist[i] = reduce_to_mcc(t, mcc)
	end
	g = trees2graph(treelist)

	# SA - Optimization
	oconfs, E, F = sa_opt(g, γ=γ, Trange=Trange, M=M)
	# Computing likelihoods
	if length(oconfs) != 1
		v() && println("Sorting $(length(oconfs)) topologically equivalent configurations.")
		vv() && println("Configurations\n", oconfs)
		vv() && println(g.labels)
		oconf = sortconf(oconfs, treelist, g, μ, mcc_names, likelihood_sort, false)
	else
		oconf = oconfs[1]
	end

	return [mcc_names[x] for x in g.labels[.!oconf]], compute_energy(oconf,g), compute_F(oconf, g, γ), E, F
end


function sortconf(oconfs, trees, g::Graph, μ, mcc_names, likelihood_sort, E_sort)
	if E_sort # Only considering configurations of lowest energies
		E = [compute_energy(conf,g) for conf in oconfs]
		Emin = minimum(E)
		oconfs_ = oconfs[findall(x->x==Emin, E)]
		v() && println("Removing ", length(oconfs) - length(oconfs_), " configurations using energy.")
	else # Removing configurations where nothing is removed
		oconfs_ = oconfs[findall(c->sum(c)<length(c), oconfs)]
	end
	# Sorting the remaining configurations using Poisson likelihood ratios
	if length(oconfs_) == 1 
		return oconfs_[1]
	elseif !likelihood_sort
		v() && println("Picking a random configuration among remaining ones")
		return rand(oconfs_)
	else
		v() && println("Comparing $(length(oconfs_)) configurations using likelihood")
		L = Float64[]
		for conf in oconfs_
			vv() && println("## Looking at configuration $conf with energy $(compute_energy(conf,g))")
			push!(L, conf_likelihood(conf, g, μ, trees, mode=:time, v=vv()))
			vv() && println()
		end
		# L = [conf_likelihood(conf, g, μ, trees) for conf in oconfs_]
		vv() && println("Confs: ", [[mcc_names[x] for x in g.labels[.!conf]] for conf in oconfs_])
		v() && println("Likelihoods: ", L)
		Lmax = maximum(L)
		oconfs_ = oconfs[findall(==(Lmax), L)]
		if length(oconfs) != 1 # Final sort by energy if more than one most likely conf
			E = [compute_energy(conf,g) for conf in oconfs_]
			Emin = minimum(E)
			oconfs_ = oconfs[findall(==(Emin), E)]
		end
		return rand(oconfs_)
	end
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
		TreeTools.remove_internal_singletons!(t, ptau=true)
	end
end
"""
	pruneconf!(trees, mcc_names, mcc_conf)

Prune MCCs `mcc_names[x]` for all `x` in `mcc_conf` from trees `t...`. 
"""
pruneconf!(trees, mcc_names, mcc_conf) = pruneconf!([mcc_names[x] for x in mcc_conf], trees...)

function pruneconf_guidetrees!(clades, trees::Vararg{Tree})
	for t in trees
		for c in clades
			TreeTools.prunenode!(t, c, propagate=true)
		end
		TreeTools.remove_internal_singletons!(t, ptau=true)
	end
end


update_df!(df::DataFrame, nleaves::Int64, nMCCs::Int64, γ, M, Efinal, Ffinal, rMCCs, arMCCs, remainingMCCs) = push!(df, 
	Dict(:nleaves=>nleaves, :nMCCs=>nMCCs, :γ=>γ, :M=>M,
		:Efinal=>Efinal, :Ffinal=>Ffinal, :removedMCCs=>rMCCs, :all_removedMCCs=>arMCCs, :remainingMCCs=>remainingMCCs))




end # module