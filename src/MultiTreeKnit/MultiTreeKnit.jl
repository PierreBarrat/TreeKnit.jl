using Dagger

"""
compute_mcc_pairs!(trees, oa)

subfunction to compute tree pair MCCs

"""
function compute_mcc_pairs!(trees::Vector{Tree{TreeTools.MiscData}}, oa::OptArgs; strict=true)
    l_t = length(trees)
    pair_MCCs = MCC_set(l_t, [t.label for t in trees])
    for r in 1:oa.rounds
        oa.verbose && @info "ROUND:"*string(r)
        for i in 1:(l_t-1)
            for j in (i+1):l_t
                joint_MCCs = nothing
                if oa.consistent && (i>1 || r>1) && l_t>1
                    if r>1
                        range = filter(e->e∉Set([i,j]), 1:l_t)
                    else
                        range = 1:(i-1)
                    end
                    for x in range
                        first = get(pair_MCCs, (i, x))
                        second = get(pair_MCCs, (j, x))
                        pre_joint_MCCs = TreeKnit.join_sets([first, second])
                        if !isnothing(joint_MCCs)
                            joint_MCCs = TreeKnit.join_sets([joint_MCCs, pre_joint_MCCs])
                        else
                            joint_MCCs = pre_joint_MCCs
                        end
                    end
                end
                TreeKnit.add!(pair_MCCs, TreeKnit.runopt(oa, trees[i], trees[j], joint_MCCs; output = :mccs), (i, j))
                if strict==false || r==oa.rounds
                    rS = TreeKnit.resolve!(trees[i], trees[j], get(pair_MCCs, (j, i)))
                    TreeTools.ladderize!(trees[i])
                    TreeKnit.sort_polytomies!(trees[i], trees[j], get(pair_MCCs, trees[i].label, trees[j].label))
                else
                    rS = TreeKnit.resolve_strict!(trees[i], trees[j], get(pair_MCCs, (j, i)))
                    TreeTools.ladderize!(trees[i])
                    TreeKnit.sort_polytomies_strict!(trees[i], trees[j], get(pair_MCCs, trees[i].label, trees[j].label))
                end
                oa.verbose && @info "found MCCs for trees: "*trees[j].label*" and "*trees[i].label
            end
        end
    end
    return pair_MCCs
end

function run_step!(oa::OptArgs, tree1::Tree, tree2::Tree, constraints, strict, r)
    if !isnothing(constraints)
        constraint = fetch(constraints[1])
        for mcc_con in constraints[2:end]
            constraint = TreeKnit.join_sets([constraint, fetch(mcc_con)])
        end
    else
        constraint = nothing
    end
    MCC = TreeKnit.runopt(oa, tree1, tree2, constraint; output = :mccs)
    if strict==false || r==oa.rounds
        rS = TreeKnit.resolve!(tree1, tree2, MCC)
        TreeTools.ladderize!(tree1)
        TreeKnit.sort_polytomies!(tree1, tree2, MCC)
    else
        rS = TreeKnit.resolve_strict!(tree1, tree2, MCC)
        TreeTools.ladderize!(tree1)
        TreeKnit.sort_polytomies_strict!(tree1, tree2, MCC)
    end
    oa.verbose && @info "found MCCs for trees: "*trees[j].label*" and "*trees[i].label
    return MCC
end

function parallelized_compute_mccs!(trees::Vector{Tree{TreeTools.MiscData}}, oa::OptArgs, strict)
    l_t = length(trees)
    parallel_MCCs = Dict()
    for r in 1:oa.rounds
        for i in 1:(l_t-1)
            for j in (i+1):l_t
                MCC_list = nothing
                if oa.consistent && (i>1 || r>1) && l_t >2
                    MCC_list = []
                    if r>1
                        range = filter(e->e∉Set([i,j]), 1:l_t)
                    else
                        range = 1:(i-1)
                    end
                    for x in range
                        append!(MCC_list, [parallel_MCCs[Set([i, x])], parallel_MCCs[Set([j,x])]])
                    end
                end
                parallel_MCCs[Set([i,j])] = Dagger.@spawn run_step!(oa, trees[i], trees[j], MCC_list, strict, r)
            end
        end
    end
    pair_MCCs = MCC_set(l_t, [t.label for t in trees])
    for (key, mcc) in parallel_MCCs
        TreeKnit.add!(pair_MCCs, fetch(mcc), Tuple(key))
    end
    return pair_MCCs
end
parallelized_compute_mccs!(trees::Vector{Tree{TreeTools.MiscData}}, oa::OptArgs) = parallelized_compute_mccs!(trees::Vector{Tree{TreeTools.MiscData}}, oa::OptArgs, true)


function compute_naive_mcc_pairs!(trees::Vector{Tree{TreeTools.MiscData}}; strict=true)
    l_t = length(trees)
    pair_MCCs = MCC_set(l_t, [t.label for t in trees])
    for i in 1:(l_t-1)
        for j in (i+1):l_t
            TreeKnit.add!(pair_MCCs, TreeKnit.naive_mccs(trees[i], trees[j]))
        end
        if strict==false
            rS = TreeKnit.resolve!(trees[i], trees[j], get(pair_MCCs, (j, i)))
            TreeTools.ladderize!(trees[i])
            TreeKnit.sort_polytomies!(trees[i], trees[j], get(pair_MCCs, trees[i].label, trees[j].label))
        else
            rS = TreeKnit.resolve_strict!(trees[i], trees[j], get(pair_MCCs, (j, i)))
            TreeTools.ladderize!(trees[i])
            TreeKnit.sort_polytomies_strict!(trees[i], trees[j], get(pair_MCCs, trees[i].label, trees[j].label))
        end
    end
    return pair_MCCs
end

"""
function get_infered_MCC_pairs!(trees::Vector{Tree{T}}, oa::OptArgs=OptArgs()) where T

function to get all MCCs of all tree pairs in tree list `trees` using TreeKnit,
resolved trees from previous MCC calculations are used as the tree input for the next pair, 
the order is specified in Combinatorics.combinations(1:length(trees), 2). If MCC pairs should be 
consistent with previous MCCs the `consistent` flag can be set to true. A `shared_branch_constraint` will then be 
computed from the previous MCC pairs that will discourage TreeKnit from removing nodes in an inconsistent manner.

For example if node `a` and node `b` are both in the same MCC clade for MCC12 and MCC13 they should also 
be together in MCC23, otherwise the MCC pairs are inconsistent.

"""
function get_infered_MCC_pairs!(trees::Vector{Tree{T}}, oa::OptArgs; strict=true, naive=false) where T

    if naive
        return compute_naive_mcc_pairs!(trees; strict)
    end

    trees = [convert(Tree{TreeTools.MiscData}, t) for t in trees]

    if oa.parallel == true
        pair_MCCs = parallelized_compute_mccs!(trees, oa, strict)
    else
        pair_MCCs = compute_mcc_pairs!(trees, oa; strict=strict)
    end

    return pair_MCCs
end

function get_infered_MCC_pairs!(trees::Vector{Tree{T}}; kwargs...) where T 
    return get_infered_MCC_pairs!(trees, OptArgs(;kwargs...))
end


###functions for calculating the constraints previously inferred MCCs put on later MCCs

"""
join_sets(input_sets::Vector{Vector{String}})

given as input a list of MCC lists this function joins these MCC lists (sets) by calculating their intersection recursively
and returning their intersection
"""
function join_sets(MCCs::Vector{Vector{Vector{String}}})
    sets = [union([Set([m... ]) for m in mcc]...) for mcc in MCCs]
    @assert union(sets...) == sets[1] ## make sure labels are the same in all trees
    mcc_map = Dict{String, Vector{Int}}()
    reverse_mcc_map = Dict{Vector{Int}, Set{String}}()
    for MCC in MCCs
        for (i,mcc) in enumerate(MCC)
            for node in mcc
                if haskey(mcc_map, node)
                    append!(mcc_map[node], i)
                else
                    mcc_map[node] = [i]
                end
            end
        end
    end
    for node in keys(mcc_map)
        if haskey(reverse_mcc_map, mcc_map[node])
            push!(reverse_mcc_map[mcc_map[node]], node)
        else
            reverse_mcc_map[mcc_map[node]] = Set([node])
        end
    end
    MCCs_new = Vector{String}[]
    for (num, mcc) in reverse_mcc_map
        append!(MCCs_new, [sort(collect(mcc))])
    end
    return TreeKnit.sort(MCCs_new; lt=TreeKnit.clt)
end

function join_sets_to_dict(MCCs::Vector{Vector{Vector{String}}})
    sets = [union([Set([m... ]) for m in mcc]...) for mcc in MCCs]
    @assert union(sets...) == sets[1] ## make sure labels are the same in all trees
    mcc_map = Dict{String, Vector{Int}}()
    reverse_mcc_map = Dict{Vector{Int}, Set{String}}()
    for MCC in MCCs
        for (i,mcc) in enumerate(MCC)
            for node in mcc
                if haskey(mcc_map, node)
                    append!(mcc_map[node], i)
                else
                    mcc_map[node] = [i]
                end
            end
        end
    end
    for node in keys(mcc_map)
        if haskey(reverse_mcc_map, mcc_map[node])
            push!(reverse_mcc_map[mcc_map[node]], node)
        else
            reverse_mcc_map[mcc_map[node]] = Set([node])
        end
    end
    MCCs_new = Dict{Int, Set{String}}()
    for (num, key) in enumerate(keys(reverse_mcc_map))
        MCCs_new[num] = reverse_mcc_map[key]
    end
    return MCCs_new
end
