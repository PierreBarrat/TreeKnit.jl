module SplitGraph

using TreeTools
using RecombTools
using StatsBase
using DataFrames

export opttrees!

include("objects.jl")
include("tools.jl")
include("energy.jl")

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
runopt(t::Vararg{Tree}; γ=1.1, M=200, 
	itmax=50, Mmax = 10_000, 
	Trange=reverse(0.01:0.2:1.1), verbose=false) = runopt(γ, M, t..., 
	itmax=round(Int64, length(first(t).lleaves)/4), 
	Trange=Trange, verbose=verbose)

function runopt(γ::Real, M::Int64, t::Vararg{Tree}; 
	Mmax = 2_000, 
	itmax=50, 
	Trange=reverse(0.01:0.2:1.1), 
	verbose=false)
	ot = deepcopy(collect(t))
	df = DataFrame(nleaves=Int64[length(ot[1].lleaves)],
		γ=Any[missing], M=Any[missing], 
		E=Any[missing], F=Any[missing])
	MCCs = []
	#
	for i in 1:itmax
		flag = :init
		verbose && println("\n --- \nIteration $i/$(itmax) - $(df.nleaves[end]) leaves remaining")
		# Optimization
		oconf, tmp, Evals, Fvals = opttrees!(ot..., γ=γ, M=M)
		if maximum([length(x) for x in oconf]) != 0
			flag = :found
		end
		# Sorting configurations
		sortedconfs = sortconf(oconf)
		mccs = sortedconfs[1][1]
		E = Evals[sortedconfs[2]][1]
		F = Fvals[sortedconfs[2]][1]

		# Checks
		!prod([check_tree(t) for t in ot]) && @error "Problem in one of the trees"

		# Actions if a solution is found
		if flag == :found && E == 0. # There is the chance that the algorithm converged
			verbose && println("$flag: temporary solution found")
			append!(MCCs, mccs)
			if sum(length(m) for m in mccs) != length(ot[1].lleaves)
				pruneconf!(mccs, ot...)
				if complete_mccs!(MCCs, ot)
					verbose && println("$flag: all mccs have been found")
					break
				end
			end
		elseif flag == :found && E != 0.
			verbose && println("$flag: temporary solution found")
			pruneconf!(mccs, ot...)
			append!(MCCs, mccs)
		end

		# Updating df based on found solution
		update_df!(df, ot[1], γ, M, E, F)

		# Actions if no solution is found
		if i == itmax
			verbose && println("Maximum number of iterations reached - stop")
			complete_mccs!(MCCs, ot, force=true)
		elseif flag != :found && M <= Mmax
			M *= 2
		elseif flag !=:found && M > Mmax
			verbose && println("Maximum M reached without converging - stop")
			complete_mccs!(MCCs, ot, force=true)
			break
		end
		
	end
	return RecombTools.sort_mccs(MCCs), df
end

"""
	runopt(t::Vararg{Tree}; γinit=5.1, dγ=0.5, γmin=1., M=200, itmax=50, Trange=reverse(0.01:0.2:1.1), verbose=false)
	runopt(γinit::Real, dγ::Real, γmin::Real, M::Int64, t::Vararg{Tree}; itmax=50, Trange=reverse(0.01:0.2:1.1), verbose=false)
"""
# runopt(t::Vararg{Tree}; γinit=4.1, dγ=0.25, γmin=1., M=200, 
# 	itmax=50, stop_when_stuck=true, 
# 	Trange=reverse(0.01:0.2:1.1), verbose=false) = runopt(γinit, dγ, γmin, M, t..., 
# 	stop_when_stuck=stop_when_stuck, itmax=itmax, 
# 	Trange=Trange, verbose=verbose)

function runopt(γinit::Real, dγ::Real, γmin::Real, M::Int64, t::Vararg{Tree}; 
	stop_when_stuck = true, 
	itmax=50, 
	Trange=reverse(0.01:0.2:1.1), 
	verbose=false)
	ot = deepcopy(collect(t))
	γi = γinit; δγ = dγ; oM = M; 
	# 
	df = DataFrame(nleaves=Int64[length(ot[1].lleaves)],
		γf=Any[missing], γi=Any[missing], γmin=Any[missing], dγ=Any[missing], M=Any[missing], 
		E=Any[missing], F=Any[missing])
	MCCs = []
	# 
	for i in 1:itmax
		verbose && println("\n --- \nIteration $i/$(itmax) - $(df.nleaves[end]) leaves remaining")
		# γrun
		oconf, Evals, Fvals, γf, flag = γrun!(γi, δγ, γmin, oM, ot...)
		verbose && println("Final γ=$γf")
		# Sorting found configurations - the sort is a bit arbitrary
		# `E` (resp. `F`) is the energy of the top conf, `sortedconfs[2]` is a permutation
		sortedconfs = sortconf(oconf)
		mccs = sortedconfs[1][1]
		E = Evals[sortedconfs[2]][1]
		F = Fvals[sortedconfs[2]][1]
		update_df!(df, ot[1], γf, γi, δγ, γmin, oM, E, F)
		# verbose && println("Energy of best conf E=$E")
		# verbose && println("Fpotential of best conf F=$F")


		# Checks
		!prod([check_tree(t) for t in ot]) && @error "Problem in one of the trees"

		# Actions depending on γrun output flag
		if flag == :stuck && stop_when_stuck
			final_mccs = maximal_coherent_clades(ot)
			append!(MCCs, final_mccs)
			break
		elseif flag == :stuck && !stop_when_stuck
			verbose && println("$flag: changing parameters")
			γi, δγ, γmin, oM = update_γrun1(γi, δγ, γmin, oM)
			#
		elseif flag == :maxit
			verbose && println("$flag: changing parameters")
			γi, δγ, γmin, oM = update_γrun2(γi, δγ, γmin, oM)
			#
		elseif flag == :found && E == 0. # There is the chance that the algorithm converged
			verbose && println("$flag: temporary solution found")
			append!(MCCs, mccs)
			if sum(length(m) for m in mccs) != length(ot[1].lleaves)
				pruneconf!(mccs, ot...)
				if complete_mccs!(MCCs, ot)
					verbose && println("$flag: all mccs have been found")
					break
				end
			end
			γi = min(γinit, γf + 8*δγ); δγ = min(δγ*2, 0.5); # Resetting parameters
			#
		elseif flag == :found && E != 0.
			verbose && println("$flag: temporary solution found")
			pruneconf!(mccs, ot...)
			append!(MCCs, mccs)
			γi = min(γinit, γf + 8*δγ); δγ = min(δγ*2, 0.5); # Resetting parameters
			#
		else
			error("Unknown flag value $flag")
		end
		#
	end
	return RecombTools.sort_mccs(MCCs), df
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
	γrun!(γi, dγ, γmin, M, t::Vararg{Tree}; verbose=false, itmax = 50)

Search a configuration (*i.e.* CCs) that reduces the potential `F(γ)` below its initial value, for diminishing values of `γ`. 
`γ` is initialized at `γi` and reduced by `dγ` at every iteration. 
The first solution found is returned.   
Return values: CCs found, energy of final confs, potential F of final confs,final `γ` value, and `flag`:  
- `:stuck`: `γ<γmin` and no solution found. Might want to increase `M` (duration of annealing) or reduce `γmin`. 
- `:maxit`: Maximum number of iterations reached.
- `:found`: At least one CC found. Energy of best configuration is not 0. 
- **not used** `:converged`: Energy of best configuration is 0 (*i.e.* full decomposition in MCCs reached).

## Note
Having a solution for a value `γ` means that the program found a CC (consistent clade) that reduces the energy (or number of split mismatches) by at least `γ`, **without** counting its own removal. In other words, removing this CC reduces the energy by `γ+1`. 
"""
function γrun!(γi, dγ, γmin, M, t::Vararg{Tree}; verbose=false, itmax = 50)
    out = [[[]], [], [], []]
    γ = γi
    flag = :unfinished
    verbose && println("γrun starting at γ=$(γ)...")
    it = 1
    while flag == :unfinished
        if γ < γmin # Stuck
        	flag = :stuck
        	break
        elseif it > itmax
        	flag = :maxit
        	break
        elseif maximum([length(x) for x in out[1]]) != 0 # Found a non empty solution 
        	flag = :found
        	break
        end

    	verbose && print("γ=$(γ)							\r")
        out = opttrees!(t..., γ=γ, M=M)
        γ -= dγ
        it += 1
    end
    verbose && println("\n Done")
    return out[1], out[3], out[4], γ, flag
end
γrun!(t::Vararg{Tree}; γi=5.6, dγ=0.25, γmin=1., 
	M=1000, verbose=false, itmax=50) = γrun!(γi, dγ, γmin, M, t..., verbose=verbose, itmax=itmax)


"""
	opttrees!(t... ; γ=1.05, Trange=0.5:-0.01:0.05, M = 1000)

Return a list of MCCs for input trees. 
Output:
1. 
"""
opttrees!(t... ; γ=1.05, Trange=reverse(0.01:0.2:1.1), M = 200) = opttrees!(γ, Trange, M, t...)
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
	return [[mcc_names[x] for x in g.labels[.!conf]] for conf in oconf], [g.labels[.!conf] for conf in oconf], [compute_energy(conf,g) for conf in oconf], [compute_F(conf, g, γ) for conf in oconf], E, F
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

"""
	sortconf(cladelist)

Sort configurations by selecting the one with the most clades. In case of equality, the one with the largest clades on average is picked. 
"""
function sortconf(cladelist)
	clt(x,y) = begin
	    if isempty(x)
	    	return false
	    elseif isempty(y)
	    	return true
	    elseif length(x) == length(y)
	    	return mean(length(a) for a in x) > mean(length(a) for a in y)
	    else
	    	return length(x) < length(y)
	    end
	end
	return sort(cladelist, lt=clt), sortperm(cladelist, lt=clt)
end

"""
"""
update_df!(df::DataFrame, t::Tree, γf, γi, dγ, γmin, M, E, F) = push!(df, Dict(:nleaves=>length(t.lleaves), 
												:γf=>γf, :γi=>γi, :γmin=>γmin, :dγ=>dγ, :M=>M,
												:E=>E, :F=>F))

update_df!(df::DataFrame, t::Tree, γ, M, E, F) = push!(df, 
	Dict(:nleaves=>length(t.lleaves), :γ=>γ, :M=>M, :E=>E, :F=>F))

"""
	update_γrun1(γ, dγ, γmin, M)

Update parameters of `γrun` for the case where the previous run was stuck with `γ<γmin`. 
"""
function update_γrun1(γi, dγ, γmin, M)
	return γi - dγ, max(0.25, dγ/2), γmin - 0.25, min(2*M, 10_000)
end
"""
	update_γrun2(γi, dγ, γmin, M)

Update parameters of `γrun` for the case where the previous run reached the maximum number of iterations
"""
function update_γrun2(γi, dγ, γmin, M)
	return γi - dγ, dγ, γmin, M
end



end