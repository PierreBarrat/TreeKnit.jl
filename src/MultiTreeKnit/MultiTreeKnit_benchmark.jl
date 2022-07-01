using TestRecombTools
using Combinatorics
using Random

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
                        nMCMC = 25, verbose=false,)
            #mCCs= TreeKnit.runopt(oa, trees[i], trees[j]; output = :all)
            #trees[i] = t_i_r
            #trees[j] = t_j_r
            mCCs = runopt(oa, trees[i], trees[j]; output = :mccs)
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
function print_MCCs(MCCs_list::Vector{Vector{Vector{String}}}, MCC_combinations_pos_to_trees_list::Vector{Vector{Int64}})
    for i in 1:length(MCCs_list)
        println("")
        println(MCC_combinations_pos_to_trees_list[i])
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
    infer_benchmark_MCCs(no_trees::Int64, lineage_number::Int64, debug=true)
"""
function infer_benchmark_MCCs(no_trees::Int64, lineage_number::Int64; debug=true)
    trees, arg = get_trees(no_trees, lineage_number);

    if debug
        ##print the trees and the true MCCs
        for i in range(1,length(trees))
            TreeTools.print_tree_ascii("", trees[i])
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
    for i in range(1, length(trees))
        tree_string= "";
        tree_string = TreeTools.write_newick!(tree_string, trees[i].root)
        append!(tree_strings, [tree_string])
        if debug
            write_newick("tree"*string(i)*".nwk", trees[i])
        end
    end
    if debug
        for i in range(1, no_trees-1)
            write_mccs("MCCs1"*string(i+1)*".dat", MCCsInfered[i])
        end
    end
    ARGPlot(tree_strings, MCCsInfered[1:(no_trees-1)], draw_connections=true, tree_names=nothing)
    if debug
        println("Found MCCs:")
        print_MCCs(MCCsInfered, MCC_combinations_pos_to_trees_list)
    end

    return trees
end


"""
    is_MCC_subset(MCC1::Vector{Vector{String}}, MCC2::Vector{Vector{String}})

    function to check that every set in MCC1 is a subset of a set in MCC2
"""
function is_MCC_subset(MCC1::Vector{Vector{String}}, MCC2::Vector{Vector{String}})
    for mcc1 in MCC1
        subset = false
        for mcc2 in MCC2
            if issubset(Set{String}(mcc1), Set{String}(mcc2))
                subset = true
                break
            end
        end
        if !subset
            return false
        end
    end
    return true
end

"""
MCC degeneracy measure
"""
function is_degenerate(no_trees, MCCs_list)
    MCC_combinations_pos_to_trees_list, MCC_combinations_trees_to_pos_dict = assign_pos_maps(no_trees) 
    k_iters = Combinatorics.combinations(1:no_trees, 3)
    for combinations in k_iters
        pos1 = MCC_combinations_trees_to_pos_dict[sort([combinations[1],combinations[2]])]
        pos2 = MCC_combinations_trees_to_pos_dict[sort([combinations[1],combinations[3]])]
        pos3 = MCC_combinations_trees_to_pos_dict[sort([combinations[3],combinations[2]])]
        join_sets([MCCs_list[pos1], MCCs_list[pos2]])
        if (!is_MCC_subset(join_sets([MCCs_list[pos1], MCCs_list[pos2]]), MCCs_list[pos3]) ||
            !is_MCC_subset(join_sets([MCCs_list[pos1], MCCs_list[pos3]]), MCCs_list[pos2]) ||
            !is_MCC_subset(join_sets([MCCs_list[pos3], MCCs_list[pos2]]), MCCs_list[pos1]))
            return true
        end
    end
    return false
end

function check_MCCs(no_trees::Int64, lineage_number::Int64; debug=true)

    MCC_combinations_pos_to_trees_list, MCC_combinations_trees_to_pos_dict = assign_pos_maps(no_trees) 
    trees, arg = get_trees(no_trees, lineage_number, remove=true);
    rMCCs = get_real_MCCs(no_trees, arg)
    RF_order = get_tree_order(trees, order="RF_distance")
    res_order = get_tree_order(trees, order="resolution")
    iMCCs = get_infered_MCC_pairs(trees, false)
    RF_iMCCs = get_infered_MCC_pairs(trees[RF_order], false)
    res_iMCCs = get_infered_MCC_pairs(trees[res_order], false)
    deg = is_degenerate(no_trees, rMCCs[1:binomial(no_trees,2)])
    ideg = is_degenerate(no_trees, iMCCs)
    RF_ideg = is_degenerate(no_trees, RF_iMCCs)
    res_ideg = is_degenerate(no_trees, res_iMCCs)

    if debug
        ##print the trees and the true MCCs
        for i in range(1,length(trees))
            TreeTools.print_tree_ascii("", trees[i])
        end
        print(rMCCs[1:binomial(no_trees,2)])
        print("\n infered:")
        print(iMCCs)
        print("\n using resolution order:")
        print(res_iMCCs)
        print("\n using RF clustering:")
        print(RF_iMCCs)
    end

    rand_index_random = Float64[] #we would like the rand index to be close to 1
    rand_index_RF = Float64[]
    rand_index_res = Float64[]
    var_index_random = Float64[] # we would like the varinfo index to be close to 0
    var_index_RF = Float64[]
    var_index_res= Float64[]

    k_iters = collect(Combinatorics.combinations(1:no_trees, 2))
    for i in range(1, binomial(no_trees,2))
        comb = k_iters[i]
        push!(rand_index_random, TestRecombTools.rand_index_similarity(rMCCs[i], iMCCs[i]))
        push!(var_index_random, TestRecombTools.varinfo_similarity(rMCCs[i], iMCCs[i]))
        pos_res = MCC_combinations_trees_to_pos_dict[sort([res_order[comb[1]], res_order[comb[2]]])]
        pos_RF = MCC_combinations_trees_to_pos_dict[sort([RF_order[comb[1]], RF_order[comb[2]]])]
        push!(rand_index_res, TestRecombTools.rand_index_similarity(rMCCs[i], res_iMCCs[pos_res]))
        push!(var_index_res, TestRecombTools.varinfo_similarity(rMCCs[i], res_iMCCs[pos_res]))
        push!(rand_index_RF, TestRecombTools.rand_index_similarity(rMCCs[i], RF_iMCCs[pos_RF]))
        push!(var_index_RF, TestRecombTools.varinfo_similarity(rMCCs[i], RF_iMCCs[pos_RF]))
    end

    return [deg, ideg, RF_ideg, res_ideg], [rand_index_random, rand_index_RF, rand_index_res], [var_index_random, var_index_RF, var_index_res]
end

