include("GenerateTrees.jl")

"""
    assign_pos_maps(no_trees::Int64) 
    ...
    Output
    - `MCC_combinations_pos_to_trees_list::Vector{Int64}`: ordered list of all tree combinations 
        (e.g. [[1,2], [1,3], [2,3], [1,2,3]] for the 3 tree case)
        from 2 pairs to all (no_trees) trees
    - `MCC_combinations_trees_to_pos_dict::Dict{Vector{Int64}, Int64}`: dictionary, key=tree combination list, 
        value= pos in MCC_combinations_pos_to_trees_list
    ...
"""
function assign_pos_maps(no_trees::Int64) 
    l = 0
        MCC_combinations_pos_to_trees_list = Vector{Int64}[]
        MCC_combinations_trees_to_pos_dict = Dict{Vector{Int64}, Int64}()
        for k in 2:no_trees
            iter_k = Combinatorics.combinations(1:no_trees, k)
            merge!(MCC_combinations_trees_to_pos_dict, Dict([(i, pos) for (i,pos) in zip(iter_k, (l + 1):(l+binomial(no_trees,k)))]))
            append!(MCC_combinations_pos_to_trees_list, collect(iter_k))
            l = length(MCC_combinations_pos_to_trees_list)
    end 
    return MCC_combinations_pos_to_trees_list, MCC_combinations_trees_to_pos_dict
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
function get_infered_MCC_pairs(trees::Vector{Tree{T}}, store_trees::Bool; consistant = true) where T
    l_t = length(trees)
    if consistant
        iter_2 = Combinatorics.combinations(1:l_t, 2)
        tree_pairs_to_pos_dict = Dict([(i, pos) for (i,pos) in zip(iter_2, 1:binomial(l_t,2))])
    end
    MCC_ordered_pairs = Vector{Vector{String}}[]
    for i in 1:(l_t-1)
        for j in (i+1):l_t
            if i>1 && consistant
                x= i
                y =j
                first = MCC_ordered_pairs[tree_pairs_to_pos_dict[[x-1, x]]]
                second = MCC_ordered_pairs[tree_pairs_to_pos_dict[[x-1, y]]]
                joint_MCCs = join_sets([first, second])
            else
                joint_MCCs = nothing
            end
            oa = OptArgs(;γ = 2., likelihood_sort = true, resolve = true,
                        nMCMC = 25, verbose=false)
            #mCCs = runopt(oa, trees[i], trees[j]; output = :mccs, constraint= joint_MCCs)
            mCCs = runopt(oa, trees[i], trees[j]; output = :mccs, constraint=joint_MCCs)
            rS = resolve!(trees[i], trees[j], mCCs)
            if store_trees
                push!(MCC_ordered_pairs, [mCCs, trees[i], trees[j]])
            else
                push!(MCC_ordered_pairs, mCCs)
            end
        end
    end
    return MCC_ordered_pairs
end



"""
    print_MCCs(MCCs_list::Vector{Vector{String}}, MCC_combinations_pos_to_trees_list::Vector{Int64})

    prints out infered MCCs 
"""
function print_MCCs(MCCs_list::Vector{Vector{Vector{String}}}, MCC_combinations_pos_to_trees_list::Vector{Vector{Int64}}, name_list::Vector{String})
    for i in 1:length(MCCs_list)
        println("")
        println(sort([name_list[MCC_combinations_pos_to_trees_list[i][k]] for k in range(1, length(MCC_combinations_pos_to_trees_list[i]))]))
        println(MCCs_list[i])
        println("")
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
function infer_benchmark_MCCs(input_trees::Vector{Tree{T}}, tree_names::Union{Nothing,Vector{String}}; debug=false, consistant=true, order="resolution") where T
    MCC_dict = Dict()
    no_trees = length(input_trees)
    if isnothing(tree_names)
        tree_names = [string(i) for i in range(1, no_trees)]
    end
    if debug
        ##print the input trees
        for i in range(1,length(input_trees))
            print(tree_names[i])
            TreeTools.print_tree_ascii("", input_trees[i])
        end
    end

    tree_order = get_tree_order(input_trees, order=order)
    trees = input_trees[tree_order]

    MCCsInfered = get_infered_MCC_pairs(trees, false, consistant=consistant)
    MCC_combinations_pos_to_trees_list, MCC_combinations_trees_to_pos_dict = assign_pos_maps(no_trees) 
    for i in range(1, length(MCCsInfered))
        MCC_dict[sort([tree_names[tree_order][j] for j in MCC_combinations_pos_to_trees_list[i]])] = MCCsInfered[i]
    end

    for k in 3:no_trees
        k_iters = Combinatorics.combinations(1:no_trees, k)
        for combination in k_iters
            all_sub_combinations = Combinatorics.combinations(combination, k-1)
            all_sets = Vector{Vector{String}}[]
            #println("MCCs of all sub combinations: \n")
            for sub_combo in all_sub_combinations
                pos = MCC_combinations_trees_to_pos_dict[sub_combo]
                mccs = MCCsInfered[pos]
                #println(mccs)
                push!(all_sets, mccs)
            end
            joint_sets = join_sets(all_sets)
            push!(MCCsInfered, joint_sets)
            MCC_dict[sort([tree_names[tree_order][j] for j in combination])] = joint_sets
        end
    end

    if debug
        tree_strings = Vector{String}()
        for i in range(1, length(trees))
            tree_string= "";
            tree_string = TreeTools.write_newick!(tree_string, trees[i].root)
            append!(tree_strings, [tree_string])
            write_newick("tree"*string(tree_names[tree_order][i])*".nwk", trees[i])
        end
        for i in range(1, no_trees-1)
            pair = sort([tree_names[tree_order][1], tree_names[tree_order][i+1]])
            write_mccs("MCCs"*string(pair[1])*","*string(pair[2])*".dat", MCCsInfered[i])
        end
        ARGPlot(tree_strings, MCCsInfered[1:(no_trees-1)], draw_connections=true, tree_names=nothing)
    end
    println("Found MCCs:")
    print_MCCs(MCCsInfered, MCC_combinations_pos_to_trees_list, tree_names[tree_order])

    return MCC_dict, input_trees
end

function infer_benchmark_MCCs(no_trees::Int64, lineage_number::Int64; debug=false, consistant=true, remove=false, order="resolution")
    trees, arg = get_trees(no_trees, lineage_number, remove=remove);
    return infer_benchmark_MCCs([trees...], nothing, debug=debug, consistant=consistant, order=order)
end

function infer_benchmark_MCCs(input_trees::Vector{Tree{T}}; debug=false, consistant=true, order="resolution") where T
    return infer_benchmark_MCCs(input_trees, nothing, debug=debug, consistant=consistant, order=order)
end

"""
get_tree_order(trees ;order="resolution")

Reorder the list of input trees according the `resolution` index (more resolved trees
are assumed to have more information and should be used first), or the `RF-distance` index
(trees that are the most similar to all other trees should be used first), or as `input`.
"""
function get_tree_order(trees ;order="resolution")
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
            permvec = sortperm(RF_index)
    else
        if order=="input"
            permvec = range(1, no_trees)
        else  ##start with most resolved trees
            resol_index = [resolution_value(t) for t in trees]
            permvec = sortperm(resol_index, rev=true)
        end
    end

    return permvec
end


