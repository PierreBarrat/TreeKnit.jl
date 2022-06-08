using Combinatorics

"""
	inferARG(
		t1::Tree, t2::Tree, oa::OptArgs = OptArgs();
		naive = false, seqlengths = [1,1],
	)
"""
function inferARG(
	t1::Tree, t2::Tree, oa::OptArgs = OptArgs();
	naive = false,
)
	MCCs = computeMCCs(t1, t2, oa)
	arg = SRG.arg_from_trees(t1, t2, MCCs)[1]
	return arg[1]
end

"""
	computeMCCs(
		t1::Tree, t2::Tree, oa::OptArgs = OptArgs();
		naive = false, seqlengths = [1,1],
	)

Compute pairwise MCCs for trees. Return MCCs and resolved splits. The `computeMCCs!`
version resolves the input trees with newly found splits.

# Inputs
### `oa::OptArgs`
Controls parameters of the MCC inference (unless `naive=true`). See `?OptArgs` for details.

In general, this should be set to `true` if more than two trees are used, and to `false`
  for only two trees (for speed).

### `naive = false`
- If `true`, use a naive estimation for MCCs, *i.e.* find all clades that have an exactly
  matching topology in all trees.
- Else, use a pseudo-parsimonious method based (mostly) on topology. The method
  `runopt(oa,t1,t2)` is called on every pair of trees.
"""
function computeMCCs(
	t1::Tree, t2::Tree, oa::OptArgs=OptArgs();
	naive=false
)
	if naive
		return naive_mccs(t1, t2)
	else
		# oac = @set oa.seq_lengths = seqlengths
		return runopt(oa, t1, t2; output = :mccs)
	end
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

function runopt(oa::OptArgs, t1::Tree, t2::Tree, tn::Vararg{Tree}; output = :mccs)
	# Copying input trees for optimization
	ot = [t for t in (t1, t2, tn...)];
	ot, copy_leaves = prepare_copies!(ot);

	# Resolve
	oa.resolve && resolve!(ot...);

	for tree_pair in Combinatorics.combinations(1:length(ot), 2)
		iMCCs = naive_mccs([ot[tree_pair[1]], ot[tree_pair[2]]], copy_leaves[tree_pair]);
		oa.verbose && @info "Initial state: $(length(iMCCs)) naive MCCs for tree pair: $(tree_pair)"
	end

	MCCs = [] # All final MCCs found up to now
	it = 1
	mcc_names = Dict()
	while true
		flag = :init
		oa.verbose && @info "--- Iteration $it (max. $(oa.itmax))"
		for copy in keys(copy_leaves)
			oa.verbose && @info "- $(length(copy_leaves[copy])) leaves remaining in tree pair $copy ---\n"
		end
		# Topology based inference
		oa.verbose && @info "Running optimization to find MCCs..."
		# TO DO: @assert share_labels(ot...) "Trees do not share leaves"
		M = [Int(ceil(length(t.lleaves) * oa.nMCMC / length(oa.Trange))) for t in ot]
		mccs, Efinal, Ffinal, lk = SplitGraph.opttrees(
			ot, copy_leaves; mcc_names=mcc_names,
			γ=oa.γ, seq_lengths = oa.seq_lengths, M=M[1], Trange=oa.Trange,
			likelihood_sort=oa.likelihood_sort, resolve=oa.resolve, sa_rep = oa.sa_rep, oa.verbose
		)
		!isempty(mccs) && append!(MCCs, mccs)
		oa.verbose && @info "Found $(length(mccs)) new mccs."

		# Stopping condition
		oa.verbose && @info "Proceeding based on newly found MCCs..."
		flag, rMCCs = stop_conditions!(MCCs, mccs, oa, it, ot[1], ot[2]; copy_leaves=copy_leaves, hardstop=true)
		#=
			Note on variables at this point
		- mccs: MCCs removed after simulated annealing (`opttrees` step)
		- rMCCs: remaining naive MCCs after removing the ones found in the `opttrees` step.
		- MCCs: all the MCCs found up to now (mccs + all the ones in the previous iters)
		=#

		# Checks
		print("copy leaves")
		print(copy_leaves)
		print(ot[1])
		print(ot[2])
		@assert prod([check_tree(t) for t in ot]) "Problem in a tree during opt."

		(flag == :stop) && break
		it += 1
	end

	if output == :all
		return TreeKnit.sort_mccs(MCCs), ot...
	else
		return sort_mccs(MCCs)
	end
end

function stop_conditions!(previous_mccs, new_mccs, oa, it, trees...; copy_leaves=nothing, hardstop=true)
	# If no solution was found
	if length(new_mccs) == 0
		# oa.verbose && print("No solution found... ")
		oa.verbose && @info "No solution found... "
		remaining_mccs = naive_mccs(trees...)
		## If hardstop or maxit,  stop
		if hardstop || it > oa.itmax
			append!(previous_mccs, remaining_mccs)
			# oa.verbose && println("Stopping.")
			oa.verbose && @info "Stopping\n"
			return :stop, []
		else
		## Otherwise, continue
			# oa.verbose && println("Continuing.")
			oa.verbose && @info "Continuing\n"
			return :next, remaining_mccs # (can't call this after because of `pruneconf!` below)
		end
	end

	# If some new MCCs were found
	## If they cover all leaves of the tree: a final decomposition has been found.
	if sum(length(m) for m in new_mccs) == length(first(trees).lleaves)
		# oa.verbose && println("Found mccs cover all leaves: final decomposition found. Stopping.")
		oa.verbose && @info "Found mccs cover all leaves: final decomposition found. Stopping.\n"
		return :stop, []
	end
	## If they do not cover all leaves of the trees, remove them from the trees
	# oa.verbose && println("Found mccs do not cover all leaves. Pruning them from trees. ")
	oa.verbose && @info "Found mccs do not cover all leaves. Pruning them from trees."
	pruneconf!(new_mccs, trees...)
	copy_leaves[[1,2]] = Set(keys(trees[1].lleaves))
	oa.resolve && resolve!(trees...)
	remaining_mccs = naive_mccs(trees...)
	### If they new trees have an obvious solution, append it to new_mccs and stop
	if length(remaining_mccs) == 1
		append!(new_mccs, remaining_mccs)
		append!(previous_mccs, remaining_mccs)
		# oa.verbose && println("Resulting trees are compatible: final decomposition found. Stopping.")
		oa.verbose && @info "Resulting trees are compatible: final decomposition found. Stopping.\n"
		return :stop, []
	### If we reached the maximum number of iterations, stop
	elseif it > oa.itmax
		append!(previous_mccs, remaining_mccs)
		# oa.verbose && println("Maximum number of iterations reached. Stopping.")
		oa.verbose && @info "Maximum number of iterations reached. Stopping.\n"
		return :stop, []
	### Otherwise, continue
	else
		# oa.verbose && println("Resulting trees have incompatibilities. Continuing.")
		oa.verbose && @info "Resulting trees have incompatibilities ($(length(remaining_mccs)) naive mccs left). Continuing.\n\n"
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


"""
	prepare_copies!(trees::Vector{Tree{T}}) where T

Prepare trees for finding MCCs, check trees share labels at start, give each node a copy mask,
each tree has l-1 copies, for tree i, n.dat["copy"][j] lets us know if there is this node is still 
included in the copy. Additionally, return the current leaves of each copy pair in a dictionary
the key is the tree pair.
"""
function prepare_copies!(trees::Vector{Tree{T}}) where T
	if !isa(trees[1], Tree{TreeTools.MiscData})
		trees = [convert(Tree{TreeTools.MiscData}, t) for t in trees]
	end
	# Check that trees share leaves
    sh = mapreduce(t->share_labels(t[1],t[2]), *, zip(trees[1:end-1], trees[2:end]))
    !sh && error("Can only be used on trees that share leaf nodes.")

	##give each node in each tree a copy mask, at the beginning every node has (l-1) copies 
	##but give a l -dimensional list for faster look-up when comparing pairs
	l = length(trees)
	for t in trees
		for n in values(t.lnodes)
			n.data = TreeTools.MiscData(Dict("copy"=>[1 for i in 1:l]))
		end
	end

	## return the leaves for each tree copy pair
	copy_leaves = Dict(c=>deepcopy(Set(keys(trees[1].lleaves))) for c in Combinatorics.combinations(1:l, 2))
	return trees, copy_leaves
end


# function update_df!(df, nleaves::Int64, nMCCs::Int64, γ, M, Efinal, Ffinal, rMCCs, arMCCs, remainingMCCs, method)
# 	push!(df,
# 		Dict(:nleaves=>nleaves,
# 		:nMCCs=>nMCCs,
# 		:γ=>γ,
# 		:M=>M,
# 		:Efinal=>Efinal,
# 		:Ffinal=>Ffinal,
# 		:newMCCs=>copy(rMCCs),
# 		:AllFinalMCCs=>copy(arMCCs),
# 		:RemainingConsistentClades=>copy(remainingMCCs),
# 		:method=>method)
# 		)
# end
