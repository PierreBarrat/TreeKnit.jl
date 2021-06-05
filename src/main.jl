"""
	computeMCCs(
		trees::Dict{<:Any, <:Tree}, oa::OptArgs=OptArgs();
		preresolve=true, naive=false
	)

Compute pairwise MCCs for `trees`. Return MCCs and resolved splits. The `computeMCCs!`
version resolved the input trees with newly found splits.

# Inputs
### `oa::OptArgs`
Controls parameters of the MCC inference (unless `naive=true`). See `?OptArgs` for details.

### `preresolve = true`
- If `true`, a first pass of MCC computation is made and trees are resolved using the
  results, keeping only compatible splits if more than two trees are given as input. A
  second pass of MCC computation is then made without resolving.
- Else, only the first pass is performed. For more than two trees, this may find MCCs that
  introduce incompatible splits in case of poorly resolved input trees.

### `naive = false`
- If `true`, use a naive estimation for MCCs, *i.e.* find all clades that have an exactly
  matching topology in all trees.
- Else, use a pseudo-parsimonious method based (mostly) on topology. The method
  `runopt(oa,t1,t2)` is called on every pair of trees.
  by calling function `runopt(oa,t1,t2)` on pairs of trees. Return MCCs and resolved splits.

### `output = :mccs`
Control the type of output. If `:mccs`, only MCCs are returned. If `:all`, splits added
	in trees during the inference process are returned as well.
"""
function computeMCCs(
	trees::Dict{<:Any, <:Tree}, oa::OptArgs=OptArgs();
	preresolve = true, naive = false, output = :mccs,
)
	ct = Dict(k=>copy(t) for (k,t) in trees)
	return computeMCCs!(ct, oa; preresolve, naive, output)
end
"""
	computeMCCs!(
		trees::Dict{<:Any, <:Tree}, oa::OptArgs=OptArgs();
		preresolve = true, naive = false, output = :mccs
	)

See `computeMCCs`.
"""
function computeMCCs!(
	trees::Dict{<:Any, <:Tree}, oa::OptArgs=OptArgs();
	preresolve = true, naive = false, output = :mccs
)
	mccs, splits = _computeMCCs!(trees, oa; preresolve, naive)
	if output == :mccs
		return mccs
	elseif output == :all
		return mccs, splits
	else
		@warn "Possible values for `output`: `:mccs` and `:all`. Using `:all`"
		return mccs, splits
	end
end
function _computeMCCs!(
	trees::Dict{<:Any, <:Tree}, oa::OptArgs=OptArgs();
	preresolve=true, naive=false
)
	if oa.crossmap_resolve
		resolve_crossmapped_muts!(trees)
	end
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
	# First pass: compute MCCs while resolving, and keep only introduced splits that
	# are compatible with all trees
	oac = @set oa.resolve = true
	oac = @set oa.output = :mccs
	resolved_splits = resolve_from_mccs!(ts -> runopt(oac,ts), trees)
	# Second pass: compute MCCs with pre-resolved trees.
	oac = @set oa.resolve = false
	MCCs = _computeMCCs(ts -> runopt(oac,ts), trees)
	return MCCs, resolved_splits
end

function computeMCCs_dynresolve!(trees::Dict{<:Any, <:Tree}, oa::OptArgs)
	MCCs = _computeMCCs(ts -> runopt(oa,ts), trees)
	resolved_splits = new_splits(trees, MCCs)
	# new_splits returns a Dict{T, Array{SplitList}}
	resolved_splits = Dict(k=>TreeTools.cat(aS...) for (k,aS) in resolved_splits)
	for k in keys(trees)
		resolve!(trees[k], resolved_splits[k])
	end
	return MCCs, resolved_splits
end

"""
	_computeMCCs(f::Function, trees::Dict{<:Any, <:Tree})

Compute MCCs for each pair of trees in `trees` by calling `f(t1,t2)`.
Return a dictionary indexed with pairs of keys of `trees`.
For simplicity, output for a pair of the same key is `[String[]]`
  (instead of not being indexed at all)
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

Run optimization at constant γ. See `?Optargs` for arguments. In the first form, keyword
  arguments are given to `OptArgs`.
"""
runopt(t1::Tree, t2::Tree; kwargs...) = runopt(OptArgs(;kwargs...), t1, t2)

function runopt(oa::OptArgs, t1::Tree, t2::Tree)
	runopt(oa, Dict(1=>t1, 2=>t2))
end
function runopt(oa::OptArgs, trees::Dict)

	#
	datatype = oa.crossmap_prune ? TreeTools.MiscData : TreeTools.EmptyData

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
			Efinal=Any[Einit],
			Ffinal=Any[Einit],
		)
		Evals = Any[]
		Fvals = Any[]
	else
		dflog = nothing
	end
	MCCs = [] # All final MCCs found up to now

	# Misc.
	SplitGraph.set_verbose(oa.verbose)
	SplitGraph.set_vverbose(oa.vv)
	SplitGraph.set_resolve(oa.resolve)

	# If using crossmapped mutations, we also need to compute self mutations for reference
	if oa.crossmap_prune
		ancestral_sequences!(ot, self=true, crossmapped=false)
	end
	for i in 1:oa.itmax
		flag = :init
		oa.verbose && println("\n --- \nIteration $i/$(oa.itmax) - $(length(leaves(first(values(ot))))) leaves remaining")

		# Find MCCs by pruning suspicious branches
		if oa.crossmap_prune
			oa.verbose && println("\n## Using cross-mapped mutations to cut branches...")
			flag = crossmap_prune!(ot, MCCs, dflog, oa)
			(flag == :stop) && break
		end

		# Topology based inference
		oa.verbose && println("\n## Running optimization to find MCCs...")
		!share_labels(values(ot)...) && error("Trees do not share leaves")
		M = getM(length(first(values(ot)).lleaves), oa.Md)
		mccs, Efinal, Ffinal, E, F, lk = SplitGraph.opttrees(
			values(ot)...;
			γ=oa.γ, seq_lengths = oa.seq_lengths, M=M, Trange=oa.Trange,
			likelihood_sort=oa.likelihood_sort, resolve=oa.resolve, sa_rep = oa.sa_rep
		)
		!isempty(mccs) && append!(MCCs, mccs)
		oa.verbose && println("Found $(length(mccs)) new mccs.")


		# Stopping condition
		oa.verbose && println("\n## Proceeding based on newly found MCCs...")
		flag, rMCCs = stop_conditions!(MCCs, mccs, oa, values(ot)...; hardstop=true)

		# Checks
		!prod([check_tree(t) for t in values(ot)]) && @error "Problem in a tree"

		# Log results
		use_df_log && update_df!(
			dflog, length(first(values(ot)).lleaves), length(rMCCs), oa.γ, M, Efinal, Ffinal,
			mccs, MCCs, rMCCs, :topology_optimization
		)
		(flag == :stop) && break
	end

	# Output
	if oa.output == :all
		return RecombTools.sort_mccs(MCCs), dflog, values(ot), Evals, Fvals
	elseif oa.output == :mccs
		return RecombTools.sort_mccs(MCCs)
	elseif oa.output == :mccs_df
		return RecombTools.sort_mccs(MCCs), dflog
	else
		error("Unknown `oa.output` $(oa.output)")
	end
end


function crossmap_prune!(ot, MCCs, logdf, oa)
	!share_labels(values(ot)...) && error("Trees do not share leaves")
	mccs = crossmap_prune(ot, oa.suspmut_threshold) # Makes a copy of the trees
	!isempty(mccs) && append!(MCCs, mccs)
	oa.verbose && println("Found $(length(mccs)) new mccs.")

	# Stop condition
	oa.verbose && println("\n## Proceeding based on newly found MCCs...")
	flag, rMCCs = stop_conditions!(MCCs, mccs, oa, values(ot)...; hardstop=false)

	# Checks
	!prod([check_tree(t) for t in values(ot)]) && @error "Problem in a tree"

	# Log results
	!isnothing(dflog) && update_df!(
		logdf, length(first(values(ot)).lleaves), length(rMCCs), missing, missing,
		missing, missing, mccs, MCCs, rMCCs, :crossmapped_mutations
	)
	#
	return flag
end

"""
"""
function crossmap_prune(trees, suspmut_threshold)
	ot = deepcopy(trees)
	# Self ancestral states can be computed once and for all at the start,
	# but we need to introduce mutations above new resolved internal nodes
	for t in values(ot), n in values(t.lnodes)
		!haskey(n.data.dat, :selfmuts) && (n.data.dat[:selfmuts] = Array{TreeTools.Mutation,1}(undef, 0))
	end
	# Cross-mapped mutations
	ancestral_sequences!(ot, self=false, crossmapped=true)
	suspicious_branches!(ot)
	#
	mccs = prune_suspicious_mccs!(deepcopy(ot), :suspicious_muts, suspmut_threshold)
end

function stop_conditions!(previous_mccs, new_mccs, oa, trees... ; hardstop=true)
	# If no solution was found
	if length(new_mccs) == 0
		oa.verbose && print("No solution found... ")
		remaining_mccs = naive_mccs(trees)
		## If hardstop,  stop
		if hardstop
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
	### Otherwise, continue
	else
		oa.verbose && println("Resulting trees have incompatibilities. Continuing.")
		return :next, remaining_mccs
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
