"""
	run_standard_treeknit!(
		trees::AbstractVector{<:Tree}, oa::OptArgs;
		infer_mccs_function = TreeKnit.runopt
	)

Subfunction to sequentially compute MCCs for each pair of trees in `trees`. Iterate over
all pairs in the same order a total of `oa.rounds` times.

## Details

Let `l_t = length(trees)`. Then for each round and for `1 <= i < j <= l_t`:

- calculate the `MCC_ij` between `trees[i]` and `trees[j]` using the standard TreeKnit
  procedure.

- resolve `trees[i]` and `trees[j]` using `MCC_ij`, unless this is the last round.

## Rounds

Going over all pairs of trees once is called a *round*. MultiTreeKnit performs several
rounds for reasons of consistency between inferred MCCs and tree topology.

The minimal number of rounds should be 2, but the effect of more than 2 rounds has not been
tested. During the final round, the trees are not resolved anymore. This can be changed
by setting `oa.final_no_resolve=false`, but it is not recommended for more than two trees.
"""
function run_standard_treeknit!(
	trees::AbstractVector{<:Tree}, oa::OptArgs;
	infer_mccs_function = TreeKnit.runopt
)
	@assert length(oa.seq_lengths) == length(trees) || length(oa.seq_lengths) == 2
	l_t = length(trees)
	pair_MCCs = MCC_set(l_t, [t.label for t in trees])
	for r in 1:oa.rounds
		@logmsg LogLevel(-1) "ROUND: $r\n\n"
		for (i, j) in _tree_pairs(l_t)
			poa = _pair_options(oa, i, j, r)
			@logmsg LogLevel(-1) """
			Infering MCCs for trees: $(trees[j].label) and $(trees[i].label)\n\n
			"""
			mccs = infer_mccs_function(poa, trees[i], trees[j]; output = :mccs)
			add!(pair_MCCs, mccs, (i, j))
			if poa.resolve
				TreeKnit.resolve!(trees[i], trees[j], mccs; poa.strict)
			end
			@logmsg LogLevel(-1) "found MCCs for trees: $(trees[j].label) and $(trees[i].label)"
			r == oa.rounds && _sort_pair!(trees, i, j, mccs, poa.strict)
		end
	end
	return pair_MCCs
end

_tree_pairs(K) = [(i, j) for i in 1:(K-1) for j in (i+1):K]
_infer_pair(f, oa, t1, t2) = f(oa, t1, t2; output = :mccs)

"""
	_pair_options(oa, i, j, round)

Copy of `oa` for inferring MCCs of trees `i` and `j` in round `round`: sequence lengths of
the two trees, and no resolution in the final round if `oa.final_no_resolve`.
The input `oa` is not modified.
"""
function _pair_options(oa::OptArgs, i, j, round)
	poa = deepcopy(oa)
	if length(oa.seq_lengths) > 2
		poa.seq_lengths = [oa.seq_lengths[i], oa.seq_lengths[j]]
	end
	if oa.final_no_resolve && round == oa.rounds
		poa.resolve = false
		poa.strict = false
	end
	return poa
end

# Ladderize the first tree and sort polytomies so that MCCs of trees `i` and `j` face each other.
function _sort_pair!(trees, i, j, mccs, strict)
	i == 1 && TreeTools.ladderize!(trees[i]) # only the first tree is ladderized
	TreeKnit.sort_polytomies!(trees[i], trees[j], mccs; strict)
	@logmsg LogLevel(-1) "ladderized and sorted trees: $(trees[j].label) and $(trees[i].label)"
end

"""
	run_parallel_treeknit!(trees, oa::OptArgs; infer_mccs_function = TreeKnit.runopt)

Parallel version of `run_standard_treeknit!`, with the same result. In rounds without
resolution, MCCs of all tree pairs are inferred in parallel. In rounds that resolve trees,
each pair modifies trees used by later pairs, so pairs are processed sequentially.
"""
function run_parallel_treeknit!(
	trees::AbstractVector{<:Tree}, oa::OptArgs;
	infer_mccs_function = TreeKnit.runopt
)
	@assert length(oa.seq_lengths) == length(trees) || length(oa.seq_lengths) == 2
	l_t = length(trees)
	pairs = _tree_pairs(l_t)
	pair_MCCs = MCC_set(l_t, [t.label for t in trees])
	for r in 1:oa.rounds
		resolving = any(((i, j),) -> _pair_options(oa, i, j, r).resolve, pairs)
		if resolving
			for (i, j) in pairs
				poa = _pair_options(oa, i, j, r)
				mccs = infer_mccs_function(poa, trees[i], trees[j]; output = :mccs)
				add!(pair_MCCs, mccs, (i, j))
				TreeKnit.resolve!(trees[i], trees[j], mccs; poa.strict)
				r == oa.rounds && _sort_pair!(trees, i, j, mccs, poa.strict)
			end
		else
			tasks = map(pairs) do (i, j)
				poa = _pair_options(oa, i, j, r)
				Dagger.@spawn _infer_pair(infer_mccs_function, poa, trees[i], trees[j])
			end
			for ((i, j), task) in zip(pairs, tasks)
				add!(pair_MCCs, fetch(task), (i, j))
			end
			if r == oa.rounds
				for (i, j) in pairs
					_sort_pair!(trees, i, j, get(pair_MCCs, (i, j)), _pair_options(oa, i, j, r).strict)
				end
			end
		end
	end
	return pair_MCCs
end

"""
	run_treeknit!(trees, oa::OptArgs; naive=false) where T
	run_treeknit!(trees; kwargs...)

Main TreeKnit run function. 

Compute MCCs of all tree pairs in tree list `trees`. The order in which pairs are processed
is the natural one given the order of trees: `(1,2), (1,3) ... (K-1, K)`.

In the second form, the keyword arguments are passed to `OptArgs`.

## Parameters:
- `pre_resolve=true`: input trees are resolved with each other prior to MCC computation.
- `resolve=true`: input trees are resolved in each pair-wise MCC computation, resolved trees are used as 
   input trees for the next pair.
- `naive`: return naive MCCs of all tree pairs.
- `strict`: apply conservative resolution. If an MCC implies a coalescence event occured, but the order of reassortment and
   coalescence is ambiguous, more than one split could be added to the tree. In such an event `strict` resolve does not add a
   split, however `liberal = (strict==false)` resolution would choose one such order of events and add that respective split. 
- `parallel`: parallelize MCC computation of tree pairs as much as possible.
"""
function run_treeknit!(trees::AbstractVector{<:Tree}, oa::OptArgs; naive=false)

    oa.pre_resolve && resolve!(trees...)

    # function used for inferring the MCCs
    # for every pair of trees (i,j), will be called as
    # `f(oa, trees[i], trees[j]; output = :mccs)`
    f = if naive
		(oa, t1, t2; kwargs...) -> TreeKnit.naive_mccs(t1, t2)
	else
		TreeKnit.runopt
	end
    return if oa.parallel == true
        run_parallel_treeknit!(trees, oa; infer_mccs_function=f)
    else
        run_standard_treeknit!(trees, oa; infer_mccs_function=f)
    end
end

function run_treeknit!(trees; naive=false, kwargs...)
	return run_treeknit!(trees, OptArgs(length(trees); kwargs...); naive)
end

function run_treeknit!(t1::Tree, t2::Tree; naive=false, kwargs...)
    return run_treeknit!([t1,t2], OptArgs(;kwargs...); naive)
end
run_treeknit!(t1::Tree, t2::Tree, oa::OptArgs) = run_treeknit!([t1,t2], oa)

run_treeknit(trees, oa::OptArgs) = run_treeknit!(map(copy, trees), oa)
run_treeknit(t1::Tree, t2::Tree, oa::OptArgs) = run_treeknit!(copy(t1), copy(t2), oa)


###############################################################################################################
####################################### ARG inference - for 1 tree pair #######################################
###############################################################################################################

"""
	inferARG(
		t1::Tree, t2::Tree, oa::OptArgs = OptArgs()
	)
"""
function inferARG(
	t1::Tree, t2::Tree, oa::OptArgs = OptArgs()
)
	oa = deepcopy(oa)
	oa.resolve = false
	oa.strict = false #liberal resolve is needed for reconstructing ARGs
	trees = [copy(t1), copy(t2)]
	MCCs = run_standard_treeknit!(trees, oa)
	return SRG.arg_from_trees(trees[1], trees[2], get(MCCs, 1, 2))[1]
end


###############################################################################################################
############################### Main inference function - for 1 tree pair #####################################
###############################################################################################################

"""
		runopt(t1::Tree, t2::Tree; kwargs...)
		runopt(oa::OptArgs, t1::Tree, t2::Tree)

Run optimization at constant γ. See `?Optargs` for arguments.
In the first form, keyword arguments are given to `OptArgs`.
"""
function runopt(t1::Tree, t2::Tree; kwargs...)
	runopt(OptArgs(2;kwargs...), t1, t2)
end

function runopt(oa::OptArgs, t1::Tree, t2::Tree; output = :mccs)
	# Copying input trees for optimization
	ot1 = copy(t1)
	ot2 = copy(t2)

	# Resolve
	oa.resolve && resolve!(ot1, ot2)

	# Initial state
	iMCCs = naive_mccs(ot1, ot2)
	@logmsg LogLevel(-1) "Initial state: $(length(iMCCs)) naive MCCs"
	MCCs = [] # All final MCCs found up to now
	it = 1
	while true
		flag = :init
		@logmsg LogLevel(-1) "--- Iteration $it (max. $(oa.itmax)) - $(length(leaves(ot1))) leaves remaining ---\n"

		# Topology based inference
		@logmsg LogLevel(-1) "Running optimization to find MCCs..."
		@logmsg LogLevel(-1) "Cooling schedule: $(oa.cooling_schedule) / Temperature values: $(length(oa.Trange)) / Total of MCMC steps: $(length(ot1.lleaves) * oa.nMCMC)"
		@assert share_labels(ot1, ot2) "Trees do not share leaves"
		M = Int(ceil(length(ot1.lleaves) * oa.nMCMC / length(oa.Trange)))
		mccs, Efinal, Ffinal, lk = SplitGraph.opttrees(
			ot1, ot2;
			oa.γ,
			oa.seq_lengths,
			M,
			oa.Trange,
			oa.likelihood_sort,
			oa.resolve,
			oa.sa_rep,
		)
		!isempty(mccs) && append!(MCCs, mccs)
		@logmsg LogLevel(-1) "Found $(length(mccs)) new mccs."


		# Stopping condition
		@logmsg LogLevel(-1) "Proceeding based on newly found MCCs..."
		flag, rMCCs = stop_conditions!(MCCs, mccs, oa, it, ot1, ot2; hardstop=true)

		#=
			Note on variables at this point
		- mccs: MCCs removed after simulated annealing (`opttrees` step)
		- rMCCs: remaining naive MCCs after removing the ones found in the `opttrees` step.
		- MCCs: all the MCCs found up to now (mccs + all the ones in the previous iters)
		=#

		# Checks
		@assert prod([check_tree(t) for t in (ot1, ot2)]) "Problem in a tree during opt."

		(flag == :stop) && break
		it += 1
	end

	if output == :all
		return TreeKnit.sort_mccs(MCCs), ot1, ot2
	else
		return sort_mccs(MCCs)
	end
end

function stop_conditions!(
	previous_mccs, new_mccs, oa, it, trees... ;
	hardstop=true
)
	# If no solution was found
	if length(new_mccs) == 0
		@logmsg LogLevel(-1) "No solution found... "
		remaining_mccs = naive_mccs(trees)
		## If hardstop or maxit,  stop
		if hardstop || it > oa.itmax
			append!(previous_mccs, remaining_mccs)
			@logmsg LogLevel(-1) "Stopping\n"
			return :stop, []
		else
		## Otherwise, continue
			@logmsg LogLevel(-1) "Continuing\n"
			return :next, remaining_mccs # (can't call this after because of `pruneconf!` below)
		end
	end

	# If some new MCCs were found
	## If they cover all leaves of the tree: a final decomposition has been found.
	if sum(length(m) for m in new_mccs) == length(first(trees).lleaves)
		@logmsg LogLevel(-1) "Found mccs cover all leaves: final decomposition found. Stopping.\n"
		return :stop, []
	end
	## If they do not cover all leaves of the trees, remove them from the trees
	@logmsg LogLevel(-1) "Found mccs do not cover all leaves. Pruning them from trees."
	@logmsg LogLevel(-2) "Pruned MCCs and trees: " new_mccs

	pruneconf!(new_mccs, trees...)
	oa.resolve && resolve!(trees...)
	remaining_mccs = naive_mccs(trees)
	### If they new trees have an obvious solution, append it to new_mccs and stop
	if length(remaining_mccs) == 1
		append!(new_mccs, remaining_mccs)
		append!(previous_mccs, remaining_mccs)
		@logmsg LogLevel(-1) "Resulting trees are compatible: final decomposition found. Stopping.\n"
		return :stop, []
	### If we reached the maximum number of iterations, stop
	elseif it > oa.itmax
		append!(previous_mccs, remaining_mccs)
		@logmsg LogLevel(-1) "Maximum number of iterations reached. Stopping.\n"
		return :stop, []
	### Otherwise, continue
	else
		@logmsg LogLevel(-1) "Resulting trees have incompatibilities ($(length(remaining_mccs)) naive mccs left). Continuing.\n\n"
		return :next, remaining_mccs
	end
end

function getM(L, Md)
	if L < 10 && Md >= 1
		return 10
	else
		return _getM(L, Md)
	end
end
_getM(n,Md) = ceil(Int, n/Md)

"""
	pruneconf!(clades, trees::Vararg{Tree})

Prune `clades` from `trees...`.
"""
function pruneconf!(clades, trees::Vararg{Tree})
	for t in trees
		for st in clades
			TreeTools.prunesubtree!(t, st, clade_only=false, create_leaf=true)
		end
		TreeTools.remove_internal_singletons!(t, delete_time=false)
	end
end
"""
	pruneconf!(trees, mcc_names, mcc_conf)

Prune MCCs `mcc_names[x]` for all `x` in `mcc_conf` from trees `t...`.
"""
function pruneconf!(trees, mcc_names, mcc_conf)
	pruneconf!([mcc_names[x] for x in mcc_conf], trees...)
end
