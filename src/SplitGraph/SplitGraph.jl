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
	runopt(t::Vararg{Tree}; γ=3., 
		M=500, itmax=50, Mmax = 1_000, Trange=reverse(0.01:0.02:1.1), 
		verbose=false, vv=false, 
		guidetrees=(), likelihood_sort=true) 

		runopt(t::Vararg{Tree}, oa::OptArgs)
		runopt(t::Vararg{Tree}) 

Run optimization at constant γ
"""
function runopt(t::Vararg{Tree}; γ=3., 
	M=ceil(Int64, length(first(t).lleaves)/50), itmax=10, Mmax = M, Trange=reverse(0.001:0.01:1.1), 
	verbose=false, vv=false, 
	guidetrees=(), likelihood_sort=true) 
	runopt(OptArgs(γ=γ, M=M, itmax=itmax, Mmax=Mmax, Trange=Trange, 
		verbose=verbose, vv=vv, guidetrees=guidetrees, likelihood_sort=likelihood_sort), 
	t...)
end

function runopt(oa::OptArgs, t::Vararg{Tree})
	#
	ot = deepcopy(collect(t))
	gt = deepcopy(collect(oa.guidetrees))
	resolve_trees!(vcat(ot,gt)...)

	nMCC = maximal_coherent_clades(ot)
	Einit = count_mismatches(ot...)
	df = DataFrame(nleaves=Int64[length(ot[1].lleaves)],
		nMCCs=length(nMCC),
		γ=Any[missing], M=Any[missing], 
		Efinal=Any[Einit], Ffinal=Any[Einit],
		removedMCCs=Any[missing], all_removedMCCs=Any[missing], remainingMCCs=Any[nMCC])
	MCCs = []
	Evals = Any[]
	Fvals = Any[]
	set_verbose(oa.verbose)
	set_vverbose(oa.vv)
	M = oa.M
	#
	for i in 1:oa.itmax
		flag = :init
		v() && println("\n --- \nIteration $i/$(oa.itmax) - $(df.nleaves[end]) leaves remaining")
		# Optimization
		mccs, Efinal, Ffinal, E, F = opttrees!(ot..., γ=oa.γ, M=ceil(Int64, length(first(ot).lleaves)/50), Trange=oa.Trange, likelihood_sort=oa.likelihood_sort)
		if length(mccs) != 0
			flag = :found
		else
			v() && println("No solution found in current iteration")
		end

		# Checks
		!prod([check_tree(t) for t in ot]) && @error "Problem in one of the trees"
		
		# Actions if a solution is found
		if flag == :found && Efinal == 0. # There is the chance that the algorithm converged
			v() && println("$flag: temporary solution found - `E == 0.`")
			append!(MCCs, mccs)
			if sum(length(m) for m in mccs) != length(ot[1].lleaves)
				pruneconf!(mccs, ot...)
				pruneconf_guidetrees!(mccs, gt...)
				if complete_mccs!(MCCs, ot)
					v() && println("$flag: all mccs have been found")
					resolve_trees!(vcat(ot,gt)...)
					nMCCs = maximal_coherent_clades(ot)
					update_df!(df, length(ot[1].lleaves), length(nMCCs), oa.γ, M, Efinal, Ffinal, mccs, MCCs, nMCCs)
					break
				end
			else # Found mccs cover all leaves
				v() && println("$flag: all mccs have been found")
				resolve_trees!(vcat(ot,gt)...)
				nMCCs = maximal_coherent_clades(ot)
				update_df!(df, length(ot[1].lleaves), length(nMCCs), oa.γ, M, Efinal, Ffinal, mccs, MCCs, nMCCs)
				break
			end
		elseif flag == :found && Efinal != 0.
			v() && println("$flag: temporary solution found - `E != 0.`")
			pruneconf!(mccs, ot...)
			pruneconf_guidetrees!(mccs, gt...)
			append!(MCCs, mccs)
		end

		# Updating df based on found solution -- will also update if no solution was found
		resolve_trees!(vcat(ot,gt)...)
		nMCCs = maximal_coherent_clades(ot)
		update_df!(df, length(ot[1].lleaves), length(nMCCs), oa.γ, M, Efinal, Ffinal, mccs, MCCs, nMCCs)
		push!(Evals, E)
		push!(Fvals, F)

		# Actions if no solution is found
		if i == oa.itmax
			v() && println("Maximum number of iterations reached - stop")
			complete_mccs!(MCCs, ot, force=true)
		elseif flag != :found && M <= oa.Mmax
			M *= 2
		elseif flag !=:found && M > oa.Mmax
			v() && println("Maximum M reached without converging - stop")
			complete_mccs!(MCCs, ot, force=true)
			break
		end
	end
	return RecombTools.sort_mccs(MCCs), df, ot, Evals, Fvals
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

	# return [[mcc_names[x] for x in g.labels[.!conf]] for conf in oconfs], [g.labels[.!conf] for conf in oconfs], [compute_energy(conf,g) for conf in oconfs], [compute_F(conf, g, γ) for conf in oconfs], E, F
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
		return oconfs_[findmax(L)[2]]
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