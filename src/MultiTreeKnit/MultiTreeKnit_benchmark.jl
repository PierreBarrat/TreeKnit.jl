include("GenerateTrees.jl")

"""
    calc_label_map(trees::Vector{Tree{T}}) where T

Create a label_map to look up which tree should be used when. 
Either an ordered list of `trees` or an ordered list of `tree_names` 
Strings are needed. This is the order the trees will be resolved in when 
computing the benchmark.
"""

function calc_label_map(trees::Vector{Tree{T}}) where T
    
    tree_names = [tree.label for tree in trees]
    label_map = calc_label_map(tree_names)
    return label_map
end

function calc_label_map(tree_names::Vector{String})
    label_map = Dict{Int, String}()
    for (i, l) in enumerate(tree_names)
        label_map[i] = l
    end
    return label_map
end

# function get_MCC_list(MCC_dict::Dict{Set{String}, Vector{Vector{String}}}(), ordered_trees::Vector{Tree{T}}) where T
#     MCC_list = Vector{Vector{String}}[]
#     label_map = calc_label_map(ordered_trees)
#     l_t = length(trees)
#     for i in 1:(l_t-1)
#         for j in i:l_t
#             append!(MCC_list, MCC_dict[Set([label_map[i], label_map[j]])])
#         end
#     end
#     return MCC_list
# end

"""
    get_infered_MCC_pairs(trees::Vector{Tree{T}}, store_trees::Bool; consistant = true) where T

    function to get all MCCs of all tree pairs in tree list `trees` using TreeKnit,
    resolved trees from previous MCC calculations are used as the tree input for the next pair, 
    the order is specified in Combinatorics.combinations(1:length(trees), 2). If MCC pairs should be 
    consistent with previous MCCs the `consistant` flag can be set to true. A `shared_branch_constraint` will then be 
    computed from the previous MCC pairs to prevent nodes from being removed in an inconsistent manner.

    For example if node `a` and node `b` are both in the same MCC clade for MCC12 and MCC13 they should also 
    be together in MCC23, otherwise the MCC pairs are inconsistent. Note that inconsistent recombination 
    events cannot be viewed together in `ARGPlot`.

"""
function get_infered_MCC_pairs!(trees::Vector{Tree{T}}; consistant = true, order="input", rev=false, constraint_cost=4., rounds=2, force=false, force_rounds=5) where T

    l_t = length(trees)

    ##if desired first change the order of the trees
    order = TreeKnit.get_tree_order(trees ;order=order, rev=rev)
    trees = trees[order]

    label_map = TreeKnit.calc_label_map(trees)
    pair_MCCs = Dict{Set{String}, Vector{Vector{String}}}()
    for r in 1:rounds
        for i in 1:(l_t-1)
            for j in (i+1):l_t
                joint_MCCs = nothing
                if consistant
                    for x in 1:l_t
                        if x ∉ Set([i, j]) && haskey(pair_MCCs, Set([label_map[x], label_map[i]])) && haskey(pair_MCCs, Set([label_map[x], label_map[j]]))
                            first = pair_MCCs[Set([label_map[i], label_map[x]])]
                            second = pair_MCCs[Set([label_map[j], label_map[x]])]
                            pre_joint_MCCs = TreeKnit.join_sets([first, second])
                            if !isnothing(joint_MCCs)
                                joint_MCCs = TreeKnit.join_sets([joint_MCCs, pre_joint_MCCs])
                            else
                                joint_MCCs = pre_joint_MCCs
                            end
                        end
                    end
                end
                pair_MCCs[Set([label_map[i], label_map[j]])] = TreeKnit.runopt(TreeKnit.OptArgs(;constraint=joint_MCCs, constraint_cost=constraint_cost), trees[i], trees[j]; output = :mccs)
                rS = TreeKnit.resolve!(trees[i], trees[j], pair_MCCs[Set([label_map[i], label_map[j]])])
            end
        end
    end
    if force
        rep = 0
        consti = consistency_rate(pair_MCCs, trees)
        while consti >0 && rep <force_rounds
            for i in 1:(l_t-1)
                for j in (i+1):l_t
                    for x in 1:l_t
                        if x ∉ Set([i, j]) 
                            first = pair_MCCs[Set([label_map[i], label_map[x]])]
                            second = pair_MCCs[Set([label_map[j], label_map[x]])]
                            third = pair_MCCs[Set([label_map[j], label_map[i]])]
                            fix_consist!([first, second, third], [trees[j], trees[i]])
                        end
                    end
                end
            end
            rep +=1
        end
        if sum(consti) >0
            print("Cannot find a consistent ARG")
        end
    end
    return pair_MCCs
end



"""
    print_MCCs(MCCs_list::Vector{Vector{String}}, MCC_combinations_pos_to_trees_list::Vector{Int64})

    prints out infered MCCs 
"""
function print_MCCs(MCCs_dict::Dict{Set{String}, Vector{Vector{String}}})
    for (key, value) in MCCs_dict
        println("")
        println(sort(collect(key)))
        println(value)
    end
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

"""
infer_benchmark_MCCs(input_trees::Vector{Tree{T}}, tree_names::Union{Nothing,Vector{String}}; debug=false, order="resolution")

Given `input_trees` with names in `tree_names` order trees either using: `"input"`, `"resolution"` (from most to least resolved)
or the `"RF_distance"` (from most like all trees to least like all other trees using the Robinson-Fould distance) and then infer 
pairwise MCCs iteratively always using the resolved trees from the last pairwise MCC calculations when calculating later trees.
To find recombination events of all tree combinations the individual tree pairs are joined. The MCCs of all possible tree 
combinations are returned in a dictionary `MCC_dict` where the sorted list of tree name combinations is the key, 
as well as the list of resolved trees `input_trees` in input order. 
"""
function infer_benchmark_MCCs!(trees::Vector{Tree{TreeTools.MiscData}}; debug=false, consistant=true, order="resolution", rev=false)
    
    ##if desired first change the order of the trees
    order = get_tree_order(trees ;order=order, rev=rev)
    trees = trees[order]

    label_map = calc_label_map(trees)
    no_trees = length(trees)

    if debug
        ##print the input trees
        for i in range(1,length(trees))
            TreeTools.print_tree_ascii("", trees[i])
        end
    end

    MCC_dict = get_infered_MCC_pairs!(trees, consistant=consistant)

    for k in 3:no_trees
        k_iters = Combinatorics.combinations(1:no_trees, k)
        for combination in k_iters
            all_sub_combinations = Combinatorics.combinations(combination, k-1)
            all_sets = Vector{Vector{String}}[]
            #println("MCCs of all sub combinations: \n")
            for sub_combo in all_sub_combinations
                mccs = MCC_dict[Set([label_map[i] for i in sub_combo])]
                push!(all_sets, mccs)
            end
            joint_sets = join_sets(all_sets)
            MCC_dict[Set([label_map[i] for i in combination])] = joint_sets
        end
    end

    if debug
        println("Found MCCs:")
        print_MCCs(MCC_dict)
    end

    return MCC_dict
end

function infer_benchmark_MCCs(no_trees::Int64, lineage_number::Int64; debug=false, consistant=true, remove=false, order="resolution")
    trees, arg = get_trees(no_trees, lineage_number, remove=remove);
    return infer_benchmark_MCCs!([trees...], debug=debug, consistant=consistant, order=order)
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


