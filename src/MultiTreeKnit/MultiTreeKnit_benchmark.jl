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

"""
    get_infered_MCC_pairs(trees::Vector{Tree{T}}, store_trees::Bool; consistant = true) where T

    function to get all MCCs of all tree pairs in tree list `trees` using TreeKnit,
    resolved trees from previous MCC calculations are used as the tree input for the next pair, 
    the order is specified in Combinatorics.combinations(1:length(trees), 2). If MCC pairs should be 
    consistent with previous MCCs the `consistant` flag can be set to true. A `mask` will then be 
    computed from the previous MCC pairs to prevent nodes from being removed in an inconsistent manner.

    For example if node `a` and node `b` are both in the same MCC clade for MCC12 and MCC13 they should also 
    be together in MCC23, otherwise the MCC pairs are inconsistent. Note that inconsistent recombination 
    events cannot be viewed together in `ARGPlot`.

"""
function get_infered_MCC_pairs!(trees::Vector{Tree{T}}; consistant = true, order="input", rev=false, constraint_cost=2.) where T

    l_t = length(trees)

    ##if desired first change the order of the trees
    order = get_tree_order(trees ;order=order, rev=rev)
    trees = trees[order]

    label_map = calc_label_map(trees)
    pair_MCCs = Dict{Set{String}, Vector{Vector{String}}}()
    for i in 1:(l_t-1)
        for j in (i+1):l_t
            if i>1 && consistant
                x= i
                y =j
                first = pair_MCCs[Set([label_map[x-1], label_map[x]])]
                second = pair_MCCs[Set([label_map[x-1], label_map[y]])]
                joint_MCCs = join_sets([first, second])
            else
                joint_MCCs = nothing
            end
            pair_MCCs[Set([label_map[i], label_map[j]])] = runopt(OptArgs(;constraint=joint_MCCs, constraint_cost=constraint_cost), trees[i], trees[j]; output = :mccs)
            rS = resolve!(trees[i], trees[j], pair_MCCs[Set([label_map[i], label_map[j]])])
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
function join_sets(input_sets::Vector{Vector{Vector{String}}})
    start_set = input_sets[1]
    for i in 2:length(input_sets)
        joint_sets = Vector{String}[]
        to_be_joint_set = [Set{String}(s2) for s2 in input_sets[i]]
        for s1 in start_set
            nodes = length(s1)
            while nodes>0
                for s2 in to_be_joint_set
                    joint_set = intersect(Set{String}(s1), s2)
                    if !isempty(joint_set)
                        s2 = setdiff(s2,joint_set)
                        nodes -= length(joint_set)
                        append!(joint_sets, [sort(collect(joint_set))])
                    end
                end
            end
        end
        start_set = joint_sets
    end
    return sort(start_set; lt=clt)
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
function infer_benchmark_MCCs!(trees::Vector{Tree{T}}; debug=false, consistant=true, order="resolution", rev=false) where T
    
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


