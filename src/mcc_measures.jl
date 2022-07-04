"""
    resolution_value(t::Tree)

    Resolution measure, ranges from 0 to 1, where 1 is a fully resolved tree and 0 is a tree with no internal nodes. 
    In a fully resolved tree with n leaves, there are (n-1) internal nodes. 
"""
function resolution_value(t::Tree)
    if length(keys(t.lleaves))==1 # if tree only contains 1 node it is resolved by definition
        return 1
    else
        return (length(keys(t.lnodes)) - length(keys(t.lleaves)))/ (length(keys(t.lleaves)) -1)
    end
end

"""
    RF_distance(t1::Tree, t2::Tree)

    Compute the Robinson–Foulds metric, or distance between two trees defined as the number of partitions 
    implied by tree 1 and not tree 2 plus the number of partitions implied by tree 2 and not tree 1
"""
function RF_distance(t1::Tree, t2::Tree)
    s1 = SplitList(t1)
    s2 = SplitList(t2)
    return length(s1) + length(s2) - 2*length(TreeTools.intersect(s1, s2))
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
    MCC degeneracy measure, the 
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


