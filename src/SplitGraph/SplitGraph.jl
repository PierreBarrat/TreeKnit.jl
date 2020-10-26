module SplitGraph

using TreeTools
using RecombTools
using StatsBase
using DataFrames
using SpecialFunctions


export opttrees!

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
	runopt(t::Vararg{Tree}; 
	γ=1.1, M=200, 
	itmax=50, stop_when_stuck=true, 
	Trange=reverse(0.01:0.2:1.1), verbose=false)

	runopt(γ::Real, M::Int64, t::Vararg{Tree}; 
	stop_when_stuck = true, 
	itmax=50, 
	Trange=reverse(0.01:0.2:1.1), 
	verbose=false)

Run optimization at constant γ
"""
runopt(t::Vararg{Tree}; γ=1.1, M=2000, 
	itmax=50, Mmax = 10_000, 
	Trange=reverse(0.01:0.2:1.1), verbose=false) = runopt(γ, M, t..., 
	itmax=itmax, Mmax=Mmax, 
	Trange=Trange, verbose=verbose)

function runopt(γ::Real, M::Int64, t::Vararg{Tree}; 
	Mmax = 100_000, 
	itmax=50, 
	Trange=reverse(0.01:0.2:1.1), 
	verbose=false)
	ot = deepcopy(collect(t))
	df = DataFrame(nleaves=Int64[length(ot[1].lleaves)],
		nMCCs=length(maximal_coherent_clades(ot)),
		γ=Any[missing], M=Any[missing], 
		Efinal=Any[count_mismatches(t...)], Ffinal=Any[count_mismatches(t...)],
		removedMCCs=Any[missing])
	MCCs = []
	Evals = Any[]
	Fvals = Any[]
	set_verbose(verbose)
	#
	for i in 1:itmax
		flag = :init
		v() && println("\n --- \nIteration $i/$(itmax) - $(df.nleaves[end]) leaves remaining")
		# Optimization
		mccs, Efinal, Ffinal, E, F = opttrees!(ot..., γ=γ, M=M, Trange=Trange)
		if length(mccs) != 0
			flag = :found
		end

		# Checks
		!prod([check_tree(t) for t in ot]) && @error "Problem in one of the trees"

		# Actions if a solution is found
		if flag == :found && Efinal == 0. # There is the chance that the algorithm converged
			v() && println("$flag: temporary solution found - `E == 0.`")
			append!(MCCs, mccs)
			if sum(length(m) for m in mccs) != length(ot[1].lleaves)
				pruneconf!(mccs, ot...)
				if complete_mccs!(MCCs, ot)
					v() && println("$flag: all mccs have been found")
					update_df!(df, ot[1], γ, M, E, F, Efinal, Ffinal, mccs)
					break
				end
			else # Found mccs cover all leaves
				v() && println("$flag: all mccs have been found")
				update_df!(df, ot[1], γ, M, E, F, Efinal, Ffinal, mccs)
				break
			end
		elseif flag == :found && Efinal != 0.
			v() && println("$flag: temporary solution found - `E != 0.`")
			pruneconf!(mccs, ot...)
			append!(MCCs, mccs)
		end

		# Updating df based on found solution
		nMCCs = length(maximal_coherent_clades(ot))
		update_df!(df, length(ot[1].lleaves), nMCCs, γ, M, Efinal, Ffinal, mccs)
		push!(Evals, E)
		push!(Fvals, F)

		# Actions if no solution is found
		if i == itmax
			v() && println("Maximum number of iterations reached - stop")
			complete_mccs!(MCCs, ot, force=true)
		elseif flag != :found && M <= Mmax
			M *= 2
		elseif flag !=:found && M > Mmax
			v() && println("Maximum M reached without converging - stop")
			complete_mccs!(MCCs, ot, force=true)
			break
		end
		
	end
	return RecombTools.sort_mccs(MCCs), df, Evals, Fvals
end


"""
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
opttrees!(t... ; γ=1.05, μ=ones(Float64, length(t)), Trange=reverse(0.01:0.2:1.1), M = 200) = opttrees!(γ, Trange, M, μ, t...)
function opttrees!(γ, Trange, M, μ, t::Vararg{Tree})
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
		oconf = sortconf(oconfs, treelist, g, μ, mcc_names)
	else
		oconf = oconfs[1]
	end

	return [mcc_names[x] for x in g.labels[.!oconf]], compute_energy(oconf,g), compute_F(oconf, g, γ), E, F

	# return [[mcc_names[x] for x in g.labels[.!conf]] for conf in oconfs], [g.labels[.!conf] for conf in oconfs], [compute_energy(conf,g) for conf in oconfs], [compute_F(conf, g, γ) for conf in oconfs], E, F
end


function sortconf(oconfs, trees, g::Graph, μ, mcc_names)
	# Only considering configurations of lowest energies
	E = [compute_energy(conf,g) for conf in oconfs]
	Emin = minimum(E)
	oconfs_ = oconfs[findall(x->x==Emin, E)]
	v() && println("Removing ", length(oconfs) - length(oconfs_), " configurations using energy.")
	# Sorting the remaining configurations using Poisson likelihood ratios
	if length(oconfs_) == 1
		return oconfs_[1]
	else
		v() && println("Comparing $(length(oconfs_)) configurations using likelihood")
		L = Float64[]
		for conf in oconfs_
			vv() && println("## Looking at configuration $conf")
			push!(L, conf_likelihood(conf, g, μ, trees, mode=:time, v=vv()))
			vv() && println()
		end
		# L = [conf_likelihood(conf, g, μ, trees) for conf in oconfs_]
		vv() && println("Confs: ", [[mcc_names[x] for x in g.labels[.!conf]] for conf in oconfs_])
		v() && println("Likelihoods: ", L)
		return oconfs[findmax(L)[2]]
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


update_df!(df::DataFrame, nleaves::Int64, nMCCs::Int64, γ, M, Efinal, Ffinal, rMCCs) = push!(df, 
	Dict(:nleaves=>nleaves, :nMCCs=>nMCCs, :γ=>γ, :M=>M,
		:Efinal=>Efinal, :Ffinal=>Ffinal, :removedMCCs=>rMCCs))




end # module