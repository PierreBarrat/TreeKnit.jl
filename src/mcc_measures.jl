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
function is_degenerate(MCC_dict::Dict{Set{String}, Vector{Vector{String}}})
    tree_names = collect([union([key for key in keys(MCC_dict)]...)]...)
    label_map = calc_label_map(tree_names)
    no_trees = length(tree_names)
    k_iters = Combinatorics.combinations(1:no_trees, 3)
    for comb in k_iters
        MCC1 = MCC_dict[Set([label_map[comb[1]], label_map[comb[2]]])]
        MCC2 = MCC_dict[Set([label_map[comb[1]], label_map[comb[3]]])]
        MCC3 = MCC_dict[Set([label_map[comb[3]], label_map[comb[2]]])]
        if is_degenerate(MCC1, MCC2, MCC3)
            return true
        end
    end
    return false
end

function is_degenerate(MCC1::Vector{Vector{String}}, MCC2::Vector{Vector{String}}, MCC3::Vector{Vector{String}})
    if (!is_MCC_subset(join_sets([MCC1, MCC2]), MCC3) || 
        !is_MCC_subset(join_sets([MCC1, MCC3]), MCC2) ||
        !is_MCC_subset(join_sets([MCC3, MCC2]), MCC1))
        return true
    else
        return false
    end
end


function consistency_rate(MCC_dict::Dict{Set{String}, Vector{Vector{String}}}, trees::Vector{Tree{T}}) where T
    score = 0
    
    label_map = calc_label_map(trees)
    l_t = length(trees)
    k_iters = Combinatorics.combinations(1:l_t, 3)
    for comb in k_iters
        MCC1 = MCC_dict[Set([label_map[comb[1]], label_map[comb[2]]])]
        MCC2 = MCC_dict[Set([label_map[comb[1]], label_map[comb[3]]])]
        MCC3 = MCC_dict[Set([label_map[comb[3]], label_map[comb[2]]])]
        s = sum(consistency_rate(MCC1, MCC2, MCC3, trees[comb]))/3
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
function consistency_rate(M12::Vector{Vector{String}}, M13::Vector{Vector{String}}, M23::Vector{Vector{String}}, trees::Vector{Tree{T}}) where T
    
    @assert length(trees) == 3
    consistency_rate = []
    MCCs = [M12, M13, M23]
    k_iters = Combinatorics.combinations(1:3, 2)
	for (i, combinations) in enumerate(k_iters)
        last = filter(e->e∉combinations,1:3)
        order = append!(combinations, last)
        tree_order = [filter(e->e∉[i],1:3)...]
        score = consistent_mcc_triplets(MCCs[order], trees[tree_order])
        append!(consistency_rate, score)
    end

	return consistency_rate
end


function consistent_mcc_triplets(MCCs, trees; masked=false)
    constraint = TreeKnit.join_sets([MCCs[1], MCCs[2]])
    s = 0
	Z = 0
    for t in trees
        tree = copy(t)
        TreeKnit.mark_shared_branches!(constraint, tree)
        for n in nodes(tree)
            if isroot(n)
                continue
            end
            if n.data.dat["shared_branch_constraint"]==true ##should be in a MCC
                Z += 1
                if !(TreeKnit.is_branch_in_mccs(n, MCCs[3]))
                    s += 1
                end
            end
        end
    end
    @assert (masked && s!=0 && TreeKnit.is_MCC_subset(TreeKnit.join_sets([MCCs[1], MCCs[2]]), MCCs[3])) == false "Error: Should be consistent" 
    return  Z == 0 ? 0.0 : s / Z
end

