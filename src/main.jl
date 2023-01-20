"""
	run_standard_treeknit!(trees::AbstractVector{Tree{T}}, oa::OptArgs, func_::Function)

Subfunction to sequentially compute MCCs for each pair of trees in `trees`. Iterate over
all pairs in the same order a total of `oa.rounds` times.

## Details

Let `l_t = length(trees)`. Then for each round and for `1 <= i < j <= l_t`:

- calculate the `MCC_ij` between `trees[i]` and `trees[j]` using the standard TreeKnit
  procedure.

- resolve `trees[i]` and `trees[j]` using `MCC_ij`, unless this is the last round.

### Rounds

Going over all pairs of trees once is called a *round*. MultiTreeKnit performs several
rounds for reasons of consistency between inferred MCCs and tree topology.

The minimal number of rounds should be 2, but the effect of more than 2 rounds has not been
tested. During the final round, the trees are not resolved anymore. This can be changed
by setting `oa.final_no_resolve=false`, but it is not recommended for more than two trees.
"""
function run_standard_treeknit!(trees::AbstractVector{Tree{T}}, oa::OptArgs; func_=TreeKnit.runopt) where T
    l_t = length(trees)
    pair_MCCs = MCC_set(l_t, [t.label for t in trees])
    for r in 1:oa.rounds
        oa.verbose && @info "ROUND:"*string(r)
        for i in 1:(l_t-1), j in (i+1):l_t

            # do we resolve for this round?
            if oa.final_no_resolve && r==oa.rounds
                oa.resolve = false
                oa.strict = false
            end

            # compute MCCs for current tree pair and add it to set
			println("Infering MCCs for trees: "*trees[j].label*" and "*trees[i].label)
            mccs = func_(oa, trees[i], trees[j]; output = :mccs)
            add!(pair_MCCs, mccs, (i, j))

            if oa.resolve 
                rS = TreeKnit.resolve!(trees[i], trees[j], get(pair_MCCs, (j, i)); oa.strict)
            end
            oa.verbose && @info "found MCCs for trees: "*trees[j].label*" and "*trees[i].label
            if r==oa.rounds 
				i == 1 && TreeTools.ladderize!(trees[i]) # only the first tree is ladderized
                TreeKnit.sort_polytomies!(trees[i], trees[j], get(pair_MCCs, trees[i].label, trees[j].label); oa.strict)
                oa.verbose && @info "ladderized and sorted trees: "*trees[j].label*" and "*trees[i].label
            end
        end
    end
    return pair_MCCs
end

function run_step!(oa::OptArgs, tree1::Tree, tree2::Tree, func_::Function, r, pos)
    if oa.final_no_resolve && r==oa.rounds
        oa = fetch(oa)
        oa.resolve = false
        oa.strict = false
        MCC = func_(oa, tree1, tree2; output = :mccs)
    else
        MCC = func_(oa, tree1, tree2; output = :mccs)
        rS = TreeKnit.resolve!(tree1, tree2, MCC; oa.strict)
    end
    if r==oa.rounds 
        if pos ==1
            TreeTools.ladderize!(tree1)
        end
        TreeKnit.sort_polytomies!(tree1, tree2, MCC; oa.strict)
    end
    return MCC
end

"""
	run_parallel_treeknit!(trees::Vector{Tree{T}}, oa::OptArgs, func_::Function) where T

Parallelized version of `run_standard_treeknit!(trees, oa, func_)`.
"""
function run_parallel_treeknit!(trees::Vector{Tree{T}}, oa::OptArgs; func_=TreeKnit.runopt) where T 
    l_t = length(trees)
    parallel_MCCs = Dict()
    for r in 1:oa.rounds
        for i in 1:(l_t-1), j in (i+1):l_t
            parallel_MCCs[Set([i,j])] = Dagger.@spawn run_step!(oa, trees[i], trees[j], func_, r, i)
        end
    end
    pair_MCCs = MCC_set(l_t, [t.label for t in trees])
    for (key, mcc) in parallel_MCCs
        add!(pair_MCCs, fetch(mcc), Tuple(key))
    end
    return pair_MCCs
end

"""
	run_treeknit!(trees::Vector{Tree{T}}, oa::OptArgs; naive=false) where T

Main TreeKnit run function. 

Computes MCCs of all tree pairs in tree list `trees` using `TreeKnit.run_opt`. 

## Parameters:
- `pre_resolve=true`: input trees are resolved with each other prior to MCC computation.
- `resolve=true`: input trees are resolved in each pair-wise MCC computation, resolved trees are used as 
input trees for the next pair, the order is specified in Combinatorics.combinations(1:length(trees), 2). 
- `naive`: return naive MCCs of all tree pairs.
- `strict`: Apply conservative resolution. If an MCC implies a coalescence event occured, but the order of reassortment and
coalescence is ambiguous, more than one split could be added to the tree. In such an event `strict` resolve does not add a
split, however `liberal = (strict==false)` resolution would choose one such order of events and add that respective split. 
- `parallel`: Parallelize MCC computation of tree pairs as much as possible.
"""
function run_treeknit!(trees::Vector{Tree{T}}, oa::OptArgs; naive=false) where T 

    oa.pre_resolve && resolve!(trees...)

    if naive
		func_ = TreeKnit.naive_mccs
	else
		func_ = TreeKnit.runopt
	end

    if oa.parallel == true
        pair_MCCs = run_parallel_treeknit!(trees, oa; func_)
    else
        pair_MCCs = run_standard_treeknit!(trees, oa; func_)
    end

    return pair_MCCs
end

function run_treeknit!(trees::Vector{Tree{T}}; kwargs...) where T 
    return run_treeknit!(trees, OptArgs(length(trees);kwargs...))
end

function run_treeknit!(t1::Tree{T}, t2::Tree{T}; kwargs...) where T 
    return run_treeknit!([t1,t2], OptArgs(2;kwargs...))
end

function run_treeknit!(t1::Tree{T}, t2::Tree{T}, oa::OptArgs) where T 
    return run_treeknit!([t1,t2], oa)
end

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
	oa.resolve= false
	oa.strict = false #liberal resolve is needed for reconstructing ARGs
	MCCs = run_standard_treeknit!(trees, oa)
	arg = SRG.arg_from_trees(t1, t2, MCCs)[1]
	return arg[1]
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
			oa.γ,
			oa.seq_lengths,
			M,
			oa.Trange,
			oa.likelihood_sort,
			oa.resolve,
			oa.sa_rep,
			oa.verbose
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

function stop_conditions!(
	previous_mccs, new_mccs, oa, it, trees... ;
	hardstop=true
)
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
