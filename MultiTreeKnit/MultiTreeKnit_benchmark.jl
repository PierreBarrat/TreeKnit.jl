using TreeKnit
using TreeTools
using Combinatorics

include("GenerateTrees.jl")
include("ARG_Plot_functions.jl")

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
    get_MCC_pairs(trees::Vector{Any}, store_trees::Bool)

    function to get all MCCs of all tree pairs in tree list `trees` using TreeKnit,
    resolved trees from previous MCC calculations are used as the tree input for the next pair, 
    the order is specified in Combinatorics.combinations(1:length(trees), 2)
"""
function get_infered_MCC_pairs(trees::Vector{Tree}, store_trees::Bool)
    MCC_ordered_pairs = Vector{Vector{String}}[]
    l_t = length(trees)
    for i in 1:(l_t-1)
        for j in (i+1):l_t
            oa = OptArgs(;γ = 2., likelihood_sort = true, resolve = true,
                        nMCMC = 25, verbose=true,)
            #mCCs= TreeKnit.runopt(oa, trees[i], trees[j]; output = :all)
            #trees[i] = t_i_r
            #trees[j] = t_j_r
            mCCs = TreeKnit.runopt(oa, trees[i], trees[j]; output = :mccs)
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
                        append!(joint_sets, Vector{Vector{String}}(sort([collect(joint_set)])))
                    end
                end
            end
        end
        start_set = joint_sets
    end
    sort!(start_set)
    return start_set
end

"""
    print_MCCs(MCCs_list::Vector{Vector{String}}, MCC_combinations_pos_to_trees_list::Vector{Int64})

    prints out infered MCCs 
"""
function print_MCCs(MCCs_list::Vector{Vector{Vector{String}}}, MCC_combinations_pos_to_trees_list::Vector{Vector{Int64}})
    for i in 1:length(MCCs_list)
        println("")
        println(MCC_combinations_pos_to_trees_list[i])
        println(MCCs_list[i])
        println("")
    end
end


"""
    infer_benchmark_MCCs(no_trees::Int64, lineage_number::Int64, debug=true)
"""
function infer_benchmark_MCCs(no_trees::Int64, lineage_number::Int64; debug=true)
    trees, MCCsReal = get_trees(no_trees, lineage_number; get_real_MCCs=true);

    if debug
        ##print the trees and the true MCCs
        for tree in trees
            TreeTools.print(tree)
        end
    end

    MCCsInfered = get_infered_MCC_pairs(trees, false)
    MCC_combinations_pos_to_trees_list, MCC_combinations_trees_to_pos_dict = assign_pos_maps(no_trees) 

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
            #println("MCCs of joint sets: \n")
            #println(joint_sets)
        end
    end

    tree_strings = Vector{String}()
    for tree in trees
        tree_string= "";
        tree_string = TreeTools.write_newick!(tree_string, tree.root)
        append!(tree_strings, [tree_string])
    end
    py"ARGPlot"(tree_strings, MCCsInfered[1:(no_trees-1)], draw_connections=true, tree_names=nothing)
    if debug
        println("Found MCCs:")
        print_MCCs(MCCsInfered, MCC_combinations_pos_to_trees_list)
    end


end

no_trees = 2
lineage_number = 6
infer_benchmark_MCCs(no_trees, lineage_number; debug=true)
println("done")


