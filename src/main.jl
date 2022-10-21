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

Compute pairwise MCCs for two trees. Return MCCs and resolved splits.

# Inputs
### `oa::OptArgs`
Controls parameters of the MCC inference (unless `naive=true`). See `?OptArgs` for details.


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
  arguments are given to `OptArgs`. If `constraint` is given (in the form of an MCC where nodes that should 
  be together are in the same cluster) this will be used as `shared_branch`-constraint while performing 
  simulated annealing, lowering the likelihood that nodes that should be in the same MCC are split from each other.
"""
runopt(t1::Tree, t2::Tree, constraint::Union{Nothing, Vector{Vector{String}}}; kwargs...) = runopt(OptArgs(;kwargs...), t1, t2, constraint)
runopt(t1::Tree, t2::Tree; kwargs...) = runopt(OptArgs(;kwargs...), t1, t2, nothing)
runopt(oa::OptArgs, t1::Tree, t2::Tree; output = :mccs) = runopt(oa, t1, t2, nothing; output=output)

function runopt(oa::OptArgs, t1::Tree, t2::Tree, constraint::Union{Nothing, Vector{Vector{String}}}; output = :mccs)
	# Copying input trees for optimization
	ot1 = copy(convert(Tree{TreeTools.MiscData}, t1))
	ot2 = copy(convert(Tree{TreeTools.MiscData}, t2))

	# Resolve
	oa.resolve && resolve!(ot1, ot2)

	format_constraint!(constraint, t1)
	mark_shared_branches!(constraint, ot1, ot2)

	iMCCs = naive_mccs(ot1, ot2)
	oa.verbose && @info "Initial state: $(length(iMCCs)) naive MCCs"

	MCCs = [] # All final MCCs found up to now

	it = 1
	while true
		flag = :init
		oa.verbose && @info "--- Iteration $it (max. $(oa.itmax)) - $(length(leaves(ot1))) leaves remaining ---\n"

		# Topology based inference
		oa.verbose && @info "Running optimization to find MCCs..."
		oa.verbose && @info "Cooling schedule: $(oa.cooling_schedule) / Temperature values: $(length(oa.Trange)) / Total of MCMC steps: $(length(ot1.lleaves) * oa.nMCMC)"
		@assert share_labels(ot1, ot2) "Trees do not share leaves"
		M = Int(ceil(length(ot1.lleaves) * oa.nMCMC / length(oa.Trange)))
		mccs, Efinal, Ffinal, lk = SplitGraph.opttrees(
			ot1, ot2;
			γ=oa.γ, seq_lengths = oa.seq_lengths, M=M, Trange=oa.Trange, likelihood_sort=oa.likelihood_sort, 
			resolve=oa.resolve, sa_rep=oa.sa_rep, consistent =oa.consistent, constraint_cost=oa.constraint_cost, oa.verbose
		)
		!isempty(mccs) && append!(MCCs, mccs)
		oa.verbose && @info "Found $(length(mccs)) new mccs."


		# Stopping condition
		oa.verbose && @info "Proceeding based on newly found MCCs..."
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

function stop_conditions!(previous_mccs, new_mccs, oa, it, trees... ; hardstop=true)
	# If no solution was found
	if length(new_mccs) == 0
		# oa.verbose && print("No solution found... ")
		oa.verbose && @info "No solution found... "
		remaining_mccs = naive_mccs(trees)
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
	oa.vv && @info "Pruned MCCs and trees: " new_mccs
	oa.vv && show(trees)
	pruneconf!(new_mccs, trees...)
	oa.resolve && resolve!(trees...)
	remaining_mccs = naive_mccs(trees)
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
function pruneconf!(trees, mcc_names, mcc_conf)
	pruneconf!([mcc_names[x] for x in mcc_conf], trees...)
end
