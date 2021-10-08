function computeMCCs(
	t1::Tree, t2::Tree, oa::OptArgs = OptArgs();
	preresolve = false, naive = false, seqlengths = [1,1],
)
	computeMCCs(Dict(1=>t1, 2=>t2), oa; preresolve, naive, seqlengths)[1,2]
end

"""
	computeMCCs(
		t1::Tree, t2::Tree, oa::OptArgs = OptArgs();
		preresolve = false, naive = false, seqlengths = [1,1],
	)
	computeMCCs(
		trees::Dict, oa::OptArgs=OptArgs();
		preresolve = false, naive = false, seqlengths = Dict(s=>1 for s in keys(trees)),
	)

Compute pairwise MCCs for trees. Return MCCs and resolved splits. The `computeMCCs!`
version resolves the input trees with newly found splits.

# Inputs
### `oa::OptArgs`
Controls parameters of the MCC inference (unless `naive=true`). See `?OptArgs` for details.

### `preresolve = false`
- If `true`, a first pass of MCC computation is made and trees are resolved using the
  results, keeping only compatible splits if more than two trees are given as input. A
  second pass of MCC computation is then made without resolving.
- Else, only the first pass is performed. For more than two trees, this may find MCCs that
  introduce incompatible splits in case of poorly resolved input trees.

In general, this should be set to `true` if more than two trees are used, and to `false`
  for only two trees (for speed).

### `naive = false`
- If `true`, use a naive estimation for MCCs, *i.e.* find all clades that have an exactly
  matching topology in all trees.
- Else, use a pseudo-parsimonious method based (mostly) on topology. The method
  `runopt(oa,t1,t2)` is called on every pair of trees.
"""
function computeMCCs(
	trees::Dict{<:Any, <:Tree}, oa::OptArgs=OptArgs();
	preresolve = false, naive = false, seqlengths = Dict(s=>1 for s in keys(trees)),
)
	ct = Dict(k=>copy(t) for (k,t) in trees)
	return computeMCCs!(ct, oa; preresolve, naive, seqlengths)
end
"""
	computeMCCs!(
		trees::Dict, oa::OptArgs=OptArgs();
		preresolve = false, naive = false, seqlengths = Dict(s=>1 for s in keys(trees)),
	)

See `computeMCCs`.
"""
function computeMCCs!(
	trees::Dict, oa::OptArgs=OptArgs();
	preresolve=false, naive=false, seqlengths = Dict(s=>1 for s in keys(trees)),
)
	if naive
		if preresolve
			return computeMCCs_naive_preresolve!(trees, oa.resolve)
		else
			return computeMCCs_naive_dynresolve!(trees, oa.resolve)
		end
	else
		oac = @set oa.seq_lengths = seqlengths
		if preresolve
			return computeMCCs_preresolve!(trees, oac)
		else
			return computeMCCs_dynresolve(trees, oac)
		end
	end
end

function computeMCCs_naive_dynresolve!(trees::Dict, resolve)
	function naive_inf!(trees::Dict, resolve)
		resolve && resolve!(values(trees)...)
		RecombTools.naive_mccs(values(trees)...)
	end
	return _computeMCCs(ts -> naive_inf!(ts, resolve), trees)
end

function computeMCCs_naive_preresolve!(trees::Dict, resolve)
	function naive_inf!(trees::Dict, resolve)
		resolve && resolve!(values(trees)...)
		RecombTools.naive_mccs(values(trees)...)
	end
	resolve_from_mccs!(ts -> naive_inf!(ts, resolve), trees)
	MCCs = _computeMCCs(ts -> naive_inf!(ts, false), trees)
	return MCCs
end

function computeMCCs_preresolve!(trees::Dict, oa::OptArgs)
	# First pass: compute MCCs while resolving, and keep only introduced splits that
	# are compatible with all trees
	# oac = @set oa.resolve = true
	oac = @set oa.output = :mccs
	oac = @set oac.verbose = false
	resolve_from_mccs!(ts -> runopt(oac,ts), trees; verbose=oa.verbose)

	# Second pass: compute MCCs with pre-resolved trees.
	oac = @set oa.resolve = false
	MCCs = _computeMCCs(ts -> runopt(oac,ts), trees)
	return MCCs
end

function computeMCCs_dynresolve(trees::Dict, oa::OptArgs)
	return _computeMCCs(ts -> runopt(oa,ts), trees)
end

"""
	_computeMCCs(f::Function, trees::Dict)

Compute MCCs for each pair of trees in `trees` by calling `f(t1,t2)`.
Return a dictionary indexed with pairs of keys of `trees`.
For simplicity, output for a pair of the same key is `[String[]]`
  (instead of not being indexed at all)
"""
function _computeMCCs(f::Function, trees::Dict)
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

Run optimization at constant γ. See `?Optargs` for arguments. In the first form, keyword
  arguments are given to `OptArgs`.
"""
runopt(t1::Tree, t2::Tree; kwargs...) = runopt(OptArgs(;kwargs...), t1, t2)

function runopt(oa::OptArgs, t1::Tree, t2::Tree)
	runopt(oa, Dict(1=>t1, 2=>t2))
end
function runopt(oa::OptArgs, trees::Dict)

	#
	datatype = TreeTools.EmptyData

	# Copying input trees for optimization
	ot = Dict(k=>copy(t, datatype) for (k,t) in trees)
	oa.resolve && resolve!(values(ot)...)

	# Writing details to a dataframe object
	use_df_log = (oa.output != :mccs)
	if use_df_log
		iMCCs = naive_mccs(values(ot)...)
		Einit = SplitGraph.count_mismatches(values(ot)...)
		n0 = length(first(values(ot)).lleaves)
		dflog = DataFrame(nleaves=Int64[n0],
			nMCCs=length(iMCCs),
			method=Any[:init],
			γ=Any[missing],
			M=Any[missing],
			newMCCs=Any[[]],
			AllFinalMCCs=Any[[]],
			RemainingConsistentClades=Any[iMCCs],
			Ffinal=Any[Einit],
			Efinal=Any[Einit]
		)
	else
		dflog = nothing
	end
	MCCs = [] # All final MCCs found up to now

	# Misc.
	SplitGraph.set_verbose(oa.verbose)
	SplitGraph.set_vverbose(oa.vv)
	SplitGraph.set_resolve(oa.resolve)

	it = 1
	while true
		flag = :init
		oa.verbose && println("\n --- \nIteration $it/$(oa.itmax) - $(length(leaves(first(values(ot))))) leaves remaining")

		# Topology based inference
		oa.verbose && println("\n## Running optimization to find MCCs...")
		!share_labels(values(ot)...) && error("Trees do not share leaves")
		M = getM(length(first(values(ot)).lleaves), oa.Md)
		mccs, Efinal, Ffinal, lk = SplitGraph.opttrees(
			values(ot)...;
			γ=oa.γ, seq_lengths = [oa.seq_lengths[x] for x in keys(ot)], M=M, Trange=oa.Trange,
			likelihood_sort=oa.likelihood_sort, resolve=oa.resolve, sa_rep = oa.sa_rep
		)
		!isempty(mccs) && append!(MCCs, mccs)
		oa.verbose && println("Found $(length(mccs)) new mccs.")


		# Stopping condition
		oa.verbose && println("\n## Proceeding based on newly found MCCs...")
		flag, rMCCs = stop_conditions!(MCCs, mccs, oa, it, values(ot)...; hardstop=true)

		# Checks
		!prod([check_tree(t) for t in values(ot)]) && @error "Problem in a tree"

		# Log results
		if use_df_log
			update_df!(
				dflog, length(first(values(ot)).lleaves), length(rMCCs), oa.γ, M, Efinal, Ffinal,
				mccs, MCCs, rMCCs, :topology_optimization
			)
		end
		(flag == :stop) && break
		it += 1
	end

	# Output
	if oa.output == :all
		return RecombTools.sort_mccs(MCCs), dflog, values(ot)
	elseif oa.output == :mccs
		return RecombTools.sort_mccs(MCCs)
	elseif oa.output == :mccs_df
		return RecombTools.sort_mccs(MCCs), dflog
	else
		error("Unknown `output` field: $(oa.output). See `?OptArgs` for allowed values.")
	end
end

function stop_conditions!(previous_mccs, new_mccs, oa, it, trees... ; hardstop=true)
	# If no solution was found
	if length(new_mccs) == 0
		oa.verbose && print("No solution found... ")
		remaining_mccs = naive_mccs(trees)
		## If hardstop or maxit,  stop
		if hardstop || it > oa.itmax
			append!(previous_mccs, remaining_mccs)
			oa.verbose && println("Stopping.")
			return :stop, []
		else
		## Otherwise, continue
			oa.verbose && println("Continuing.")
			return :next, remaining_mccs # (can't call this after because of `pruneconf!` below)
		end
	end

	# If some new MCCs were found
	## If they cover all leaves of the tree: a final decomposition has been found.
	if sum(length(m) for m in new_mccs) == length(first(trees).lleaves)
		oa.verbose && println("Found mccs cover all leaves: final decomposition found. Stopping.")
		return :stop, []
	end
	## If they do not cover all leaves of the trees, remove them from the trees
	oa.verbose && println("Found mccs do not cover all leaves. Pruning them from trees. ")
	pruneconf!(new_mccs, trees...)
	oa.resolve && resolve!(trees...)
	remaining_mccs = naive_mccs(trees)
	### If they new trees have an obvious solution, append it to new_mccs and stop
	if length(remaining_mccs) == 1
		append!(new_mccs, remaining_mccs)
		append!(previous_mccs, remaining_mccs)
		oa.verbose && println("Resulting trees are compatible: final decomposition found. Stopping.")
		return :stop, []
	### If we reached the maximum number of iterations, stop
	elseif it > oa.itmax
		append!(previous_mccs, remaining_mccs)
		oa.verbose && println("Maximum number of iterations reached. Stopping.")
		return :stop, []
	### Otherwise, continue
	else
		oa.verbose && println("Resulting trees have incompatibilities. Continuing.")
		return :next, remaining_mccs
	end
end

function getM(L, Md)
	if L < 10 && Md >= 1
		return 10
	else
		return _getM(L, Md)
	end
	tmp
end
_getM(n,Md) = ceil(Int, n/Md)

"""
	pruneconf!(clades, trees::Vararg{Tree})

Prune `clades` from `trees...`.
"""
function pruneconf!(clades, trees::Vararg{Tree})
	for t in trees
		for st in clades
			TreeTools.prunesubtree!(t, st, clade_only=false)
		end
		TreeTools.remove_internal_singletons!(t, ptau=true)
	end
end
"""
	pruneconf!(trees, mcc_names, mcc_conf)

Prune MCCs `mcc_names[x]` for all `x` in `mcc_conf` from trees `t...`.
"""
pruneconf!(trees, mcc_names, mcc_conf) = pruneconf!([mcc_names[x] for x in mcc_conf], trees...)



function update_df!(df::DataFrame, nleaves::Int64, nMCCs::Int64, γ, M, Efinal, Ffinal, rMCCs, arMCCs, remainingMCCs, method)
	push!(df,
		Dict(:nleaves=>nleaves,
		:nMCCs=>nMCCs,
		:γ=>γ,
		:M=>M,
		:Efinal=>Efinal,
		:Ffinal=>Ffinal,
		:newMCCs=>copy(rMCCs),
		:AllFinalMCCs=>copy(arMCCs),
		:RemainingConsistentClades=>copy(remainingMCCs),
		:method=>method)
		)
end
