using Dagger
include("GenerateTrees.jl")

"""
compute_mcc_pairs!(trees, oa)

subfunction to compute tree pair MCCs

"""
function compute_mcc_pairs!(trees::Vector{Tree{TreeTools.MiscData}}, oa::OptArgs; strict=false)
    l_t = length(trees)
    pair_MCCs = MCC_set(l_t, [t.label for t in trees])
    for r in 1:oa.rounds
        oa.verbose && @info "ROUND:"*string(r)
        for i in 1:(l_t-1)
            for j in (i+1):l_t
                joint_MCCs = nothing
                if oa.consistent && (i>1 || r>1)
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
                TreeKnit.add!(pair_MCCs, TreeKnit.runopt(TreeKnit.OptArgs(;constraint_cost=oa.constraint_cost), trees[i], trees[j], joint_MCCs; output = :mccs), (i, j))
                rS = TreeKnit.resolve!(trees[i], trees[j], get(pair_MCCs, (j, i)); strict=strict)
                TreeTools.ladderize!(trees[i])
                if strict == false
                    TreeKnit.sort_polytomies!(trees[i], trees[j], get(pair_MCCs, trees[i].label, trees[j].label))
                end
                oa.verbose && @info "found MCCs for trees: "*trees[j].label*" and "*trees[i].label
            end
        end
    end
    return pair_MCCs
end

function run_step!(oa::OptArgs, tree1::Tree, tree2::Tree, constraints, strict)
    if !isnothing(constraints)
        constraint = fetch(constraints[1])
        for mcc_con in constraints[2:end]
            constraint = TreeKnit.join_sets([constraint, fetch(mcc_con)])
        end
    else
        constraint = nothing
    end
    MCC = TreeKnit.runopt(oa, tree1, tree2, constraint; output = :mccs)
    rS = TreeKnit.resolve!(tree1, tree2, MCC; strict=strict)
    TreeTools.ladderize!(tree1)
    if strict==false
        TreeKnit.sort_polytomies!(tree1, tree2, MCC)
    end
    oa.verbose && @info "found MCCs for trees: "*trees[j].label*" and "*trees[i].label
    return MCC
end

function parallelized_compute_mccs!(trees::Vector{Tree{TreeTools.MiscData}}, oa::OptArgs; strict=false)
    l_t = length(trees)
    parallel_MCCs = Dict()
    for r in 1:oa.rounds
        for i in 1:(l_t-1)
            for j in (i+1):l_t
                MCC_list = nothing
                if oa.consistent && (i>1 || r>1)
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
                parallel_MCCs[Set([i,j])] = Dagger.@spawn run_step!(oa, trees[i], trees[j], MCC_list, strict)
            end
        end
    end
    pair_MCCs = MCC_set(l_t, [t.label for t in trees])
    for (key, mcc) in parallel_MCCs
        TreeKnit.add!(pair_MCCs, fetch(mcc), Tuple(key))
    end
    return pair_MCCs
end

"""
function get_infered_MCC_pairs!(trees::Vector{Tree{T}}, oa::OptArgs=OptArgs()) where T

function to get all MCCs of all tree pairs in tree list `trees` using TreeKnit,
resolved trees from previous MCC calculations are used as the tree input for the next pair, 
the order is specified in Combinatorics.combinations(1:length(trees), 2). If MCC pairs should be 
consistent with previous MCCs the `consistent` flag can be set to true. A `shared_branch_constraint` will then be 
computed from the previous MCC pairs to prevent nodes from being removed in an inconsistent manner.

For example if node `a` and node `b` are both in the same MCC clade for MCC12 and MCC13 they should also 
be together in MCC23, otherwise the MCC pairs are inconsistent. Note that inconsistent recombination 
events cannot be viewed together in `ARGPlot`.

"""

function get_infered_MCC_pairs!(trees::Vector{Tree{T}}, oa::OptArgs; strict=false) where T

    trees = [convert(Tree{TreeTools.MiscData}, t) for t in trees]

    if oa.parallel == true
        pair_MCCs = parallelized_compute_mccs!(trees, oa; strict=strict)
    else
        pair_MCCs = compute_mcc_pairs!(trees, oa; strict=strict)
    end

    if oa.force_consist
        oa.verbose && @info "Fix consistency"
        pair_MCCs = fix_consist_sets!(pair_MCCs, trees; topo=oa.force_topo_consist)
    end
    return pair_MCCs
end

function get_infered_MCC_pairs!(trees::Vector{Tree{T}}; kwargs...) where T 
    return get_infered_MCC_pairs!(trees, OptArgs(;kwargs...))
end


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

"""
infer_benchmark_MCCs(input_trees::Vector{Tree{T}}, tree_names::Union{Nothing,Vector{String}}; debug=false, order="resolution")

Given `input_trees` with names in `tree_names` order trees either using: `"input"`, `"resolution"` (from most to least resolved)
or the `"RF_distance"` (from most like all trees to least like all other trees using the Robinson-Fould distance) and then infer 
pairwise MCCs iteratively always using the resolved trees from the last pairwise MCC calculations when calculating later trees.
To find recombination events of all tree combinations the individual tree pairs are joined. The MCCs of all possible tree 
combinations are returned in a MCC_set object `MCC_dict` where the sorted list of tree name combinations is the key, 
as well as the list of resolved trees `input_trees` in input order. 
"""
function infer_benchmark_MCCs!(trees::Vector{Tree{TreeTools.MiscData}}; debug=false, consistent=true, order="input", rev=false)

    ##if desired first change the order of the trees
    order = get_tree_order(trees ;order=order, rev=rev)
    trees = trees[order]

    if debug
        ##print the input trees
        for i in range(1,length(trees))
            TreeTools.print_tree_ascii("", trees[i])
        end
    end

    MCC_dict = get_infered_MCC_pairs!(trees; consistent=consistent)

    for k in 3:MCC_dict.no_trees
        k_iters = Combinatorics.combinations(1:MCC_dict.no_trees, k)
        for combination in k_iters
            all_sub_combinations = Combinatorics.combinations(combination, k-1)
            all_sets = Vector{Vector{String}}[]
            #println("MCCs of all sub combinations: \n")
            for sub_combo in all_sub_combinations
                mccs = get(MCC_dict,Tuple([i for i in sub_combo]))
                push!(all_sets, mccs)
            end
            joint_sets = join_sets(all_sets)
            add!(MCC_dict, joint_sets, Tuple([i for i in combination]))
        end
    end

    if debug
        print(MCC_dict)
    end

    return MCC_dict
end

function infer_benchmark_MCCs(no_trees::Int64, lineage_number::Int64; debug=false, consistent=true, remove=false, order="resolution")
    trees, arg = get_trees(no_trees, lineage_number, remove=remove);
    return infer_benchmark_MCCs!([trees...], debug=debug, consistent=consistent, order=order)
end

"""
get_tree_order(trees ;order="resolution")

Reorder the list of input trees according the `resolution` index (more resolved trees
are assumed to have more information and should be used first), or the `RF-distance` index
(trees that are the most similar to all other trees should be used first), or as `input`.
"""
function get_tree_order(trees ;order="resolution", rev=false)
    no_trees = length(trees)
    if order=="RF_distance"
        ##start with trees that are most similar to all other trees
            RF_index = Float16[]
            for i in range(1, no_trees)
                rf = 0
                for j in range(1, no_trees)
                    if j!=i
                        rf += RF_distance(trees[i], trees[j])^2
                    end
                end
                push!(RF_index, rf/(no_trees-1))
            end
            permvec = sortperm(RF_index, rev=rev)
    else
        if order=="input"
            permvec = range(1, no_trees)
        else  ##start with most resolved trees
            resol_index = [resolution_value(t) for t in trees]
            permvec = sortperm(resol_index, rev=!rev)
        end
    end

    return permvec
end


