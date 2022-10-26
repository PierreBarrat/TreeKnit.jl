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

function is_MCC_subset_dict(MCC1::Dict{Int, Set{String}}, MCC2::Dict{String, Union{Int, Nothing}})
    for (keys, mcc1) in MCC1
        if length(Set([MCC2[m] for m in mcc1]))!=1
            return false
        end
    end
    return true
end

"""
is_degenerate(M::MCC_set)

Function to check if consistency holds for all MCC subtriplets in the MCC set.
This uses only the MCC sets themselves and does not take tree topology into consideration.
Returns a boolean value.
"""
function is_degenerate(M::MCC_set)
    k_iters = Combinatorics.combinations(1:M.no_trees, 3)
    for comb in k_iters
        MCC1 = get(M, comb[1], comb[2])
        MCC2 = get(M, comb[1], comb[3])
        MCC3 = get(M, comb[2], comb[3])
        if is_degenerate(MCC1, MCC2, MCC3)
            return true
        end
    end
    return false
end

function is_degenerate(MCC1::Vector{Vector{String}}, MCC2::Vector{Vector{String}}, MCC3::Vector{Vector{String}})
    if (!is_MCC_subset_dict(MTK.MCC_join_constraint([MCC1, MCC2]; dict=true), map_mccs_leaves(MCC3)) || 
        !is_MCC_subset_dict(MTK.MCC_join_constraint([MCC1, MCC3]; dict=true), map_mccs_leaves(MCC2)) ||
        !is_MCC_subset_dict(MTK.MCC_join_constraint([MCC3, MCC2]; dict=true), map_mccs_leaves(MCC1)))
        return true
    else
        return false
    end
end


function consistency_rate(M::MCC_set, trees::Vector{Tree{T}}) where T
    @assert M.order_trees == [t.label for t in trees]
    score = 0

    k_iters = Combinatorics.combinations(1:M.no_trees, 3)
    for comb in k_iters
        MCC1 = get(M, comb[1], comb[2])
        MCC2 = get(M, comb[1], comb[3])
        MCC3 = get(M, comb[2], comb[3])
        s = sum(consistency_rate(MCC1, MCC2, MCC3, trees[comb]))/3
        score += s
    end
    return score / binomial(M.no_trees,3)
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
    constraint = MTK.MCC_join_constraint([MCCs[1], MCCs[2]])
    s = 0
	Z = 0
    for t in trees
        tree = copy(t)
        MTK.mark_shared_branches!(constraint, tree)
        for n in nodes(tree)
            if isroot(n)
                continue
            end
            if n.data.dat["shared_branch"]==true ##should be in a MCC
                Z += 1
                if !(TreeKnit.is_branch_in_mccs(n, MCCs[3]))
                    s += 1
                end
            end
        end
    end
    @assert (masked && s!=0 && TreeKnit.is_MCC_subset(MTK.MCC_join_constraint([MCCs[1], MCCs[2]]), MCCs[3])) == false "Error: Should be consistent" 
    return  Z == 0 ? 0.0 : s / Z
end

"""
accuracy_shared_branches(tree, true_tree, MCC, rMCC)

Calculate the number of TP (true positive), FP, TN, and FN (false negative)
for shared branch inference. Label the true underlying tree `true_tree` (prior to branch removal)
using the real MCCs `MCC`, do the same with the infered, resolved tree `tree` (assuming TreeKnit received
an unresolved version of `true_tree` with removed branches and branches were then added to that tree using 
the inferred MCCs) and inferred MCCs `MCC`. For each branch that exists in both trees (use splitlist)
check if it has been correctly labeled as shared or not shared.
"""
function accuracy_shared_branches(tree, true_tree, MCC, rMCC)
    true_positive = 0
    false_positive = 0
    true_negative = 0
    false_negative = 0
    ##note the true tree should be fully resolved
    TreeKnit.mark_shared_branches!(MCC, tree)
    TreeKnit.mark_shared_branches!(rMCC, true_tree)
    true_splits = SplitList(true_tree)
    true_splits_dict = Dict()
    for n in nodes(true_tree)
        if !isleaf(n) && !isroot(n)
            split = true_splits.splitmap[n.label]
            true_splits_dict[split] = n
        end
    end
    s1 = SplitList(tree)
    for n in nodes(tree)
        node_in_true_tree = nothing
        if !isleaf(n) && !isroot(n)
            split = s1.splitmap[n.label]
            if split ∈ SplitList(true_tree)
                ##true split, branch exists in true tree
                node_in_true_tree = true_splits_dict[split]
            end
        end
        if isleaf(n)
            node_in_true_tree = true_tree.lnodes[n.label]
        end
        if !isnothing(node_in_true_tree)
            if n.data["shared_branch"] 
                if node_in_true_tree.data["shared_branch"]
                    true_positive += 1
                else
                    false_positive += 1
                end
            else
                if node_in_true_tree.data["shared_branch"]
                    false_negative += 1
                else
                    true_negative += 1
                end
            end
        end
    end
    return true_positive, false_positive, false_negative, true_negative
end


"""
is_full_topo_compatible(t1::Tree{T}, t2::Tree{T}, mcc::Vector{String}) where T

Check if the MCC `mcc` is fully compatible (up to resolution) with the topology of `t1` and `t2`.
After 1 standard round of TreeKnit this will always be the case, but for sequential TreeKnit on multiple trees,
where trees are additionally resolved during each pair-wise iteration this may not always be the case.
"""
function is_full_topo_compatible(t1::Tree{T}, t2::Tree{T}, mcc::Vector{String}) where T
	mask = [leaf in mcc for leaf in sort(collect(keys(t1.lleaves)))]
	S1 = TreeTools.unique(TreeTools.SplitList(t1, mask))
	mask = [leaf in mcc for leaf in sort(collect(keys(t2.lleaves)))]
	S2 = TreeTools.unique(TreeTools.SplitList(t2, mask))
    filtered_values = filter(i->mask[i], 1:length(mask))
	for s1 in S1.splits
        if length(filter(x->x in filtered_values, s1.dat))> 1
            found = false
            for s2 in S2.splits
                if TreeTools.isequal(s1, s2, mask)
                    found = true
                    break
                end
            end
            if !found
                return false
            end
        end
	end
	return true
end

"""
is_topo_compatible(t1::Tree{T}, t2::Tree{T}, mcc::Vector{String}) where T

Check if the MCC `mcc` is compatible (with additional - yet not performed- resolution) with the topology of `t1` and `t2`.
After 1 standard round of TreeKnit this will always be the case, but for sequential TreeKnit on multiple trees,
where trees are additionally resolved during each pair-wise iteration this may not always be the case.
"""
function is_topo_compatible(t1::Tree{T}, t2::Tree{T}, mcc::Vector{String}) where T
	mask = [leaf in mcc for leaf in sort(collect(keys(t1.lleaves)))]
	S1 = TreeTools.SplitList(t1, mask)
	mask = [leaf in mcc for leaf in sort(collect(keys(t2.lleaves)))]
	S2 = TreeTools.SplitList(t2, mask)
	for s1 in S1.splits
		for s2 in S2.splits
			if !arecompatible(s1, s2, S2.mask)
				return false
			end
		end
	end
	return true
end




