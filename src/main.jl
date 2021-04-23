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

Compute pairwise MCCs for `trees` by calling function `runopt(oa,t1,t2)` on pairs of trees. Return MCCs and resolved splits. 
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
	MCCs = _computeMCCs(ts -> RecombTools.naive_mccs(collect(values(ts))...), trees)
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


###############################################################################################################
########################################## Main inference function ############################################
###############################################################################################################

"""
		runopt(t1::Tree, t2::Tree; kwargs...)
		runopt(oa::OptArgs, t1::Tree, t2::Tree)
		runopt(oa::OptArgs, trees::Dict{<:Any,<:Tree})

Run optimization at constant γ. See `?Optargs` for arguments. 
"""
runopt(t1::Tree, t2::Tree; kwargs...) = runopt(OptArgs(;kwargs...), t1, t2)

function runopt(oa::OptArgs, t1::Tree, t2::Tree)
	runopt(oa, Dict(1=>t1, 2=>t2))
end
function runopt(oa::OptArgs, trees::Dict)
	# 
	ot = deepcopy(collect(values(trees)))
	oa.resolve && resolve!(ot...)
	# 
	iMCCs = naive_mccs(ot)
	Einit = SplitGraph.count_mismatches(ot...)
	n0 = length(first(ot).lleaves)
	df = DataFrame(nleaves=Int64[n0],
		nMCCs=length(iMCCs),
		γ=Any[oa.γ], M=Any[missing], 
		Efinal=Any[Einit], Ffinal=Any[Einit],
		removedMCCs=Any[missing], all_removedMCCs=Any[missing], remainingMCCs=Any[iMCCs])
	MCCs = []
	Evals = Any[]
	Fvals = Any[]
	SplitGraph.set_verbose(oa.verbose)
	SplitGraph.set_vverbose(oa.vv)
	SplitGraph.set_resolve(oa.resolve)
	#
	for i in 1:oa.itmax
		flag = :init
		oa.verbose && println("\n --- \nIteration $i/$(oa.itmax) - $(df.nleaves[end]) leaves remaining")

		# Prune suspicious branches (if mut. cross-mapping)
		if oa.crossmap_prune
			mccs = RecombTools.prune_suspicious_mccs!(Dict(s=>t for (s,t) in zip(collect(keys(trees)), ot)), :suspicious_muts) # Need to give keys in a sensible way 
			if !isempty(mccs)
				rMCCs = naive_mccs(ot)
				append!(MCCs, mccs)
				update_df!(df, length(ot[1].lleaves), length(rMCCs), oa.γ, missing, missing, missing, mccs, MCCs, rMCCs)
			end
		end
		## Loop above
		## Should be a function
		## Should take care of doing ancestral reconstruction

		# Optimization
		n = length(first(ot).lleaves)
		M = getM(n, oa.Md)
		mccs, Efinal, Ffinal, E, F, lk = SplitGraph.opttrees!(ot..., γ=oa.γ, seq_lengths = oa.seq_lengths,
			M=M, Trange=oa.Trange, likelihood_sort=oa.likelihood_sort, resolve=oa.resolve, sa_rep = oa.sa_rep)
		length(mccs) != 0 ? (flag = :found) : (oa.verbose && println("No solution found in current iteration"))
		append!(MCCs, mccs)

		# Checks
		!prod([check_tree(t) for t in ot]) && @error "Problem in one of the trees"
		
		# Actions if a final solution is found
		if flag == :found 
			if sum(length(m) for m in mccs) != length(ot[1].lleaves) # Found mccs do not cover all leaves
				pruneconf!(mccs, ot...)
				if complete_mccs!(MCCs, ot)
					oa.verbose && println("$flag: all mccs have been found")
					oa.resolve && resolve!(ot...)
					rMCCs = naive_mccs(ot)
					update_df!(df, length(ot[1].lleaves), length(rMCCs), oa.γ, M, Efinal, Ffinal, mccs, MCCs, rMCCs)
					break
				end
			else # Found mccs cover all leaves
				oa.verbose && println("$flag: all mccs have been found")
				oa.resolve && resolve!(ot...)
				rMCCs = naive_mccs(ot)
				update_df!(df, length(ot[1].lleaves), length(rMCCs), oa.γ, M, Efinal, Ffinal, mccs, MCCs, rMCCs)
				break
			end
		end

		# Action if a temporary solution is found (pruneconf! is called above)
		oa.resolve && resolve!(ot...)
		rMCCs = naive_mccs(ot)
		update_df!(df, length(ot[1].lleaves), length(rMCCs), oa.γ, M, Efinal, Ffinal, mccs, MCCs, rMCCs)
		push!(Evals, E)
		push!(Fvals, F)

		# Actions if no solution is found
		if i == oa.itmax || flag != :found
			oa.verbose && println("Maximum number of iterations reached or no solution found - stop")
			complete_mccs!(MCCs, ot, force=true)
			break
		end
	end
	if oa.output == :all
		return RecombTools.sort_mccs(MCCs), df, ot, Evals, Fvals
	elseif oa.output == :mccs
		return RecombTools.sort_mccs(MCCs)
	elseif oa.output == :mccs_df
		return RecombTools.sort_mccs(MCCs), df
	else
		error("Unknown `oa.output` $(oa.output)")
	end
end


"""
	complete_mccs!(MCCs, ot; force=false)

Complete `MCCs` if `naive_mccs(ot)` is of length `1`, and return `true`. Otherwise (unless `force`), return `false`.
"""
function complete_mccs!(MCCs, ot; force=false)
	final_mccs = naive_mccs(ot)
	if length(final_mccs) != 1
		force && append!(MCCs, final_mccs)
		return false
	else
		append!(MCCs, final_mccs)
		return true
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



update_df!(df::DataFrame, nleaves::Int64, nMCCs::Int64, γ, M, Efinal, Ffinal, rMCCs, arMCCs, remainingMCCs) = push!(df, 
	Dict(:nleaves=>nleaves, :nMCCs=>nMCCs, :γ=>γ, :M=>M,
		:Efinal=>Efinal, :Ffinal=>Ffinal, :removedMCCs=>rMCCs, :all_removedMCCs=>arMCCs, :remainingMCCs=>remainingMCCs))
