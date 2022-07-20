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
        if (!is_MCC_subset(join_sets([MCCs_list[pos1], MCCs_list[pos2]]), MCCs_list[pos3]) ||
            !is_MCC_subset(join_sets([MCCs_list[pos1], MCCs_list[pos3]]), MCCs_list[pos2]) ||
            !is_MCC_subset(join_sets([MCCs_list[pos3], MCCs_list[pos2]]), MCCs_list[pos1]))
            return true
        end
    end
    return false
end


function consistency_rate(MCCs, trees)
    l_t = length(trees)
    score = 0
    if binomial(l_t,2) != length(MCCs)
        return "Error MCC list does not align with number of trees"
    end
    MCC_combinations_pos_to_trees_list, MCC_combinations_trees_to_pos_dict = assign_pos_maps(l_t) 
    k_iters = Combinatorics.combinations(1:l_t, 3)
    for combinations in k_iters
        pos1 = MCC_combinations_trees_to_pos_dict[sort([combinations[1],combinations[2]])]
        pos2 = MCC_combinations_trees_to_pos_dict[sort([combinations[1],combinations[3]])]
        pos3 = MCC_combinations_trees_to_pos_dict[sort([combinations[3],combinations[2]])]
        s = consistency_rate(MCCs[pos1], MCCs[pos2], MCCs[pos3], trees[combinations])
        score += s
    end

	return score / binomial(l_t,3)
end

"""
    Calculate the average consistency of MCCs of a triplet of trees, each mcc triplet combination is evaluated
    and the consistency score is averaged. The consistency score for one triplet combination e.g. M12, M13 and M23
    is calculated as the total number of branches that are inconsistent (i.e. not grouped as would be expected)
    divided by the total number of branches. 

"""
function consistency_rate(M12, M13, M23, trees)
    consistency_rate = 0
    MCCs = [M12, M13, M23]
    k_iters = Combinatorics.combinations(1:3, 2)
	for (i, combinations) in enumerate(k_iters)
        last = filter(e->e∉combinations,1:3)
        order = append!(combinations, last)
        tree_order = [i, filter(e->e∉[i],1:3)...]
        score = consistent_mcc_triplets(MCCs[order], trees[tree_order])
        consistency_rate += score
    end

	return consistency_rate / 3
end


function consistent_mcc_triplets(MCCs, trees)
    s = 0
	Z = 0
    for n in nodes(trees[1])
        if isroot(n)
            continue
        end
        if !isnothing(TreeKnit.find_mcc_with_branch(n, MCCs[1])) && !isnothing(TreeKnit.find_mcc_with_branch(n, MCCs[2]))
            Z += 1
            # Must find a corresponding common branch in trees[2] (or trees[3])
            # We find it by taking the intesection of the clade of n with the MCC n
            # belongs to
            m12 = find_mcc_with_branch(n, MCCs[1])[2]
            m13 = find_mcc_with_branch(n, MCCs[2])[2]
            clade = [x.label for x in POTleaves(n)]
            n2 = lca(trees[2], intersect(clade, m12))
            n3 = lca(trees[3], intersect(clade, m13))
            if !(!isnothing(TreeKnit.find_mcc_with_branch(n2, MCCs[3])) || isroot(n2)) || !(!isnothing(TreeKnit.find_mcc_with_branch(n3, MCCs[3])) || isroot(n3))
                s += 1
            end
        end
    end
    return  Z == 0 ? 0.0 : s / Z
end

