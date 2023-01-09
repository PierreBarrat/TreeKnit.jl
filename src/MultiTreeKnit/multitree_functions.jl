"""
	compute_mcc_pairs!(trees, oa)

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
function compute_mcc_pairs!(trees::Vector{Tree{T}}, oa::OptArgs) where T 
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
            mccs = TreeKnit.runopt(oa, trees[i], trees[j]; output = :mccs)
            add!(pair_MCCs, mccs, (i, j))

            if oa.resolve 
                rS = TreeKnit.resolve!(trees[i], trees[j], get(pair_MCCs, (j, i)); oa.strict)
            end
            oa.verbose && @info "found MCCs for trees: "*trees[j].label*" and "*trees[i].label
            if r==oa.rounds 
                if i ==1 ##only the first tree should be ladderized
                    TreeTools.ladderize!(trees[i])
                end
                TreeKnit.sort_polytomies!(trees[i], trees[j], get(pair_MCCs, trees[i].label, trees[j].label); oa.strict)
                oa.verbose && @info "ladderized and sorted trees: "*trees[j].label*" and "*trees[i].label
            end
        end
    end
    return pair_MCCs
end

function run_step!(oa::OptArgs, tree1::Tree, tree2::Tree, r, pos)
    if oa.final_no_resolve && r==oa.rounds
        oa = fetch(oa)
        oa.resolve = false
        oa.strict = false
        MCC = TreeKnit.runopt(oa, tree1, tree2; output = :mccs)
    else
        MCC = TreeKnit.runopt(oa, tree1, tree2; output = :mccs)
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
	parallelized_compute_mccs!(trees::Vector{Tree{T}}, oa::OptArgs, strict) where T

Parallelized version of `compute_mcc_pairs!(trees, oa)`
"""
function parallelized_compute_mccs!(trees::Vector{Tree{T}}, oa::OptArgs) where T 
    l_t = length(trees)
    parallel_MCCs = Dict()
    for r in 1:oa.rounds
        for i in 1:(l_t-1), j in (i+1):l_t
            parallel_MCCs[Set([i,j])] = Dagger.@spawn run_step!(oa, trees[i], trees[j], r, i)
        end
    end
    pair_MCCs = MCC_set(l_t, [t.label for t in trees])
    for (key, mcc) in parallel_MCCs
        add!(pair_MCCs, fetch(mcc), Tuple(key))
    end
    return pair_MCCs
end

function compute_naive_mcc_pairs!(trees::Vector{Tree{T}}; oa::OptArgs) where T 
    l_t = length(trees)
    pair_MCCs = MCC_set(l_t, [t.label for t in trees])
    for r in 1:oa.rounds
        for i in 1:(l_t-1), j in (i+1):l_t
            if oa.final_no_resolve && r==oa.rounds
                oa.resolve = false
                oa.strict = false
            end
            TreeKnit.add!(pair_MCCs, TreeKnit.naive_mccs(trees[i], trees[j]))
            if oa.resolve
                rS = TreeKnit.resolve!(trees[i], trees[j], get(pair_MCCs, (j, i)); oa.strict)
            end
            if i==1
                TreeTools.ladderize!(trees[i])
            end
            TreeKnit.sort_polytomies!(trees[i], trees[j], get(pair_MCCs, trees[i].label, trees[j].label); strict)
        end
    end
    return pair_MCCs
end

"""
    get_infered_MCC_pairs!(trees::Vector{Tree{T}}, oa::OptArgs=OptArgs()) where T

Function to compute MCCs of all tree pairs in tree list `trees` using `TreeKnit.run_opt`. Resolved trees
from previous MCC calculations are used as the input trees for the next pair, the order is specified in
Combinatorics.combinations(1:length(trees), 2). 

## Parameters:
- `naive`: return naive MCCs of all tree pairs.
- `strict`: Apply conservative resolution. If a reassortment event happened in a polytomy that could be resolved by a MCC 
it is unclear if the reassortment or the coalescence happened first, in order to introduce such a split the order must be 
randomly assigned. If this is desired `strict` should be set to false.
- `oa.parallel`: Parallelize MCC computation of tree pairs as much as possible.
"""
function get_infered_MCC_pairs!(trees::Vector{Tree{T}}, oa::OptArgs; naive=false) where T 

    oa.pre_resolve && resolve!(trees...)

    if naive
        return compute_naive_mcc_pairs!(trees; oa)
    end

    if oa.parallel == true
        pair_MCCs = parallelized_compute_mccs!(trees, oa)
    else
        pair_MCCs = compute_mcc_pairs!(trees, oa)
    end

    return pair_MCCs
end

function get_infered_MCC_pairs!(trees::Vector{Tree{T}}; kwargs...) where T 
    return get_infered_MCC_pairs!(trees, OptArgs(;kwargs...))
end