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

function is_MCC_subset_dict(MCC1::Dict{Int, Set{String}}, MCC2::Dict{String, Int})
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
    if (!is_MCC_subset_dict(join_sets_to_dict([MCC1, MCC2]), get_mcc_map(MCC3)) || 
        !is_MCC_subset_dict(join_sets_to_dict([MCC1, MCC3]), get_mcc_map(MCC2)) ||
        !is_MCC_subset_dict(join_sets_to_dict([MCC3, MCC2]), get_mcc_map(MCC1)))
        return true
    else
        return false
    end
end

function is_topologically_degenerate(M::MCC_set, trees::Vector{Tree{T}}) where T 
    k_iters = Combinatorics.combinations(1:M.no_trees, 3)
    for comb in k_iters
        MCC1 = get(M, comb[1], comb[2])
        MCC2 = get(M, comb[1], comb[3])
        MCC3 = get(M, comb[2], comb[3])
        if (is_topologically_degenerate(MCC1, MCC2, MCC3, trees[comb[2]], trees[comb[3]]) ||
            is_topologically_degenerate(MCC1, MCC3, MCC2, trees[comb[1]], trees[comb[3]]) ||
            is_topologically_degenerate(MCC2, MCC3, MCC1, trees[comb[1]], trees[comb[2]]))
            return true
        end
    end
    return false
end

function is_topologically_degenerate(MCC1::Vector{Vector{String}}, MCC2::Vector{Vector{String}}, MCC3::Vector{Vector{String}}, tree2::Tree{T}, tree3::Tree{T}) where T
    ##if mcc12 is a strictly smaller (not equal) subset of mcc13 this means that there must be a recombination event above mcc12 in mcc23 as well (and that this can be drawn onto tree 2)
    needed_splits_12 = Dict()
    needed_splits_13 = Dict()
    for mcc1 in MCC1
        for mcc2 in MCC2
            if issubset(Set{String}(mcc1), Set{String}(mcc2)) && !isempty(setdiff(Set{String}(mcc2), Set{String}(mcc1))) #check if mcc1 is a notequal subset of mcc2
                if haskey(needed_splits_12, mcc2)
                    append!(needed_splits_12[mcc2], [mcc1])
                else
                    needed_splits_12[mcc2] = [mcc1]
                end
            elseif issubset(Set{String}(mcc2), Set{String}(mcc1)) && !isempty(setdiff(Set{String}(mcc1), Set{String}(mcc2))) #check if mcc1 is a notequal subset of mcc2
                if haskey(needed_splits_13, mcc1)
                    append!(needed_splits_13[mcc1], [mcc2])
                else
                    needed_splits_13[mcc1] = [mcc2]
                end
            end
        end
    end
    if isempty(needed_splits_12) && isempty(needed_splits_13)
        return false
    end
    TreeKnit.mark_shared_branches!(MCC3, tree2)
    for (key, n) in needed_splits_12
        # split_needed = [!isroot(TreeTools.lca(tree2, s...)) && TreeTools.lca(tree2, s...).data.dat["shared_branch_constraint"]==true for s in n]
        # no_split_needed = count(x->x==true,split_needed)
        # if no_split_needed >1 && union([Set(l) for l in n]...) == Set(key)
        #     return true
        # elseif no_split_needed >0 && union([Set(l) for l in n]...) != Set(key)
        #     return true
        # end

        for s in n
            if !isroot(TreeTools.lca(tree2, s...)) && TreeTools.lca(tree2, s...).data.dat["shared_branch_constraint"]==true
            ##need to split
            return true
            end
        end

    end
    TreeKnit.mark_shared_branches!(MCC3, tree3)
    for (key, n) in needed_splits_13
        # split_needed = [!isroot(TreeTools.lca(tree3, s...)) && TreeTools.lca(tree3, s...).data.dat["shared_branch_constraint"]==true for s in n]
        # no_split_needed = count(x->x==true,split_needed)
        # if no_split_needed >1 && union([Set(l) for l in n]...) == Set(key)
        #     return true
        # elseif no_split_needed >0 && union([Set(l) for l in n]...) != Set(key)
        #     return true
        # end

        for s in n
            if !isroot(TreeTools.lca(tree3, s...)) && TreeTools.lca(tree3, s...).data.dat["shared_branch_constraint"]==true
            ##need to split
            return true
            end
        end
    end
    return false
end

function count_new_splits_for_topo_degeneracy(M::MCC_set, trees::Vector{Tree{T}}) where T 
    count = 0
    k_iters = Combinatorics.combinations(1:M.no_trees, 3)
    for comb in k_iters
        MCC1 = get(M, comb[1], comb[2])
        MCC2 = get(M, comb[1], comb[3])
        MCC3 = get(M, comb[2], comb[3])
        count += count_new_splits_for_topo_degeneracy(MCC1, MCC2, MCC3, trees[comb[2]], trees[comb[3]])
        count += count_new_splits_for_topo_degeneracy(MCC1, MCC3, MCC2, trees[comb[1]], trees[comb[3]])
        count += count_new_splits_for_topo_degeneracy(MCC2, MCC3, MCC1, trees[comb[1]], trees[comb[2]])
    end
    return count
end

function count_new_splits_for_topo_degeneracy(MCC1::Vector{Vector{String}}, MCC2::Vector{Vector{String}}, MCC3::Vector{Vector{String}}, tree2::Tree{T}, tree3::Tree{T}) where T
    ##if mcc12 is a strictly smaller (not equal) subset of mcc13 this means that there must be a recombination event above mcc12 in mcc23 as well (and that this can be drawn onto tree 2)
    needed_splits_12 = []
    needed_splits_13 = []
    for mcc1 in MCC1
        for mcc2 in MCC2
            if issubset(Set{String}(mcc1), Set{String}(mcc2)) && !isempty(setdiff(Set{String}(mcc2), Set{String}(mcc1))) #check if mcc1 is a notequal subset of mcc2
                append!(needed_splits_12, [mcc1])
            elseif issubset(Set{String}(mcc2), Set{String}(mcc1)) && !isempty(setdiff(Set{String}(mcc1), Set{String}(mcc2))) #check if mcc1 is a notequal subset of mcc2
                append!(needed_splits_13, [mcc2])
            end
        end
    end
    if isempty(needed_splits_12) && isempty(needed_splits_13)
        return 0
    end
    count = 0
    TreeKnit.mark_shared_branches!(MCC3, tree2)
    for n in needed_splits_12
        if !isroot(TreeTools.lca(tree2, n...)) && TreeTools.lca(tree2, n...).data.dat["shared_branch_constraint"]==true
            ##need to split
            count +=1 
        end
    end
    TreeKnit.mark_shared_branches!(MCC3, tree3)
    for n in needed_splits_13
        if !isroot(TreeTools.lca(tree3, n...)) && TreeTools.lca(tree3, n...).data.dat["shared_branch_constraint"]==true
            ##need to split
            count +=1 
        end
    end
    return count
end

function fix_topologically_degenerate(MCC1::Vector{Vector{String}}, MCC2::Vector{Vector{String}}, MCC3::Vector{Vector{String}}, tree2::Tree{T}, tree3::Tree{T}) where T
    ##if mcc12 is a strictly smaller (not equal) subset of mcc13 this means that there must be a recombination event above mcc12 in mcc23 as well (and that this can be drawn onto tree 2)
    needed_splits_12 = Dict()
    needed_splits_13 = Dict()
    for mcc1 in MCC1
        for mcc2 in MCC2
            if issubset(Set{String}(mcc1), Set{String}(mcc2)) && !isempty(setdiff(Set{String}(mcc2), Set{String}(mcc1))) #check if mcc1 is a notequal subset of mcc2
                if haskey(needed_splits_12, mcc2)
                    append!(needed_splits_12[mcc2], [mcc1])
                else
                    needed_splits_12[mcc2] = [mcc1]
                end
            elseif issubset(Set{String}(mcc2), Set{String}(mcc1)) && !isempty(setdiff(Set{String}(mcc1), Set{String}(mcc2))) #check if mcc1 is a notequal subset of mcc2
                if haskey(needed_splits_13, mcc1)
                    append!(needed_splits_13[mcc1], [mcc2])
                else
                    needed_splits_13[mcc1] = [mcc2]
                end
            end
        end
    end
    TreeKnit.mark_shared_branches!(MCC3, tree2)
    TreeKnit.mark_shared_branches!(MCC3, tree3)
    mcc_map = get_mcc_map(MCC3)
    next_i = length(MCC3) +1 
    for (key, n) in needed_splits_12
        for s in n
            if !isroot(TreeTools.lca(tree2, s...)) && TreeTools.lca(tree2, s...).data.dat["shared_branch_constraint"]==true
                node = TreeTools.lca(tree2, s...)
                ##there should be a split here but there is not one -> need to split
                mcc = node.data.dat["mcc"]
                for l in POT(node)
                    if l.data.dat["mcc"] == mcc
                        l.data.dat["mcc"] = next_i
                        if isleaf(l)
                            mcc_map[l.label] = next_i
                        end
                    end
                end
                next_i +=1 
            end
        end
        # split_needed = [!isroot(TreeTools.lca(tree2, s...)) && TreeTools.lca(tree2, s...).data.dat["shared_branch_constraint"]==true for s in n]
        # no_split_needed = count(x->x==true,split_needed)
        # if union([Set(l) for l in n]...) != Set(key)
        #     no_split_needed +=1 
        # end
        # i = 1
        # while no_split_needed > 1
        #     if split_needed[i] == true
        #         node = TreeTools.lca(tree2, n[i]...)
        #         ##there should be a split here but there is not one -> need to split
        #         mcc = node.data.dat["mcc"]
        #         for l in POT(node)
        #             if l.data.dat["mcc"] == mcc
        #                 l.data.dat["mcc"] = next_i
        #                 if isleaf(l)
        #                     mcc_map[l.label] = next_i
        #                 end
        #             end
        #         end
        #         next_i +=1 
        #         no_split_needed -=1
        #         ## else there is a split here, nothing needed
        #     else
        #         i +=1
        #     end
        # end
    end
    for (key, n) in needed_splits_13
        for s in n
            if !isroot(TreeTools.lca(tree3, s...)) && TreeTools.lca(tree3, s...).data.dat["shared_branch_constraint"]==true
                node = TreeTools.lca(tree3, s...)
                ##there should be a split here but there is not one -> need to split
                mcc = node.data.dat["mcc"]
                for l in POT(node)
                    if l.data.dat["mcc"] == mcc
                        l.data.dat["mcc"] = next_i
                        if isleaf(l)
                            mcc_map[l.label] = next_i
                        end
                    end
                end
                next_i +=1 
            end
        end

        # split_needed = [!isroot(TreeTools.lca(tree3, s...)) && TreeTools.lca(tree3, s...).data.dat["shared_branch_constraint"]==true for s in n]
        # no_split_needed = count(x->x==true,split_needed)
        # if union([Set(l) for l in n]...) != Set(key)
        #     no_split_needed +=1 
        # end
        # i = 1
        # while no_split_needed > 1
        #     if split_needed[i] == true
        #         node = TreeTools.lca(tree3, n[1]...)
        #         ##there should be a split here but there is not one -> need to split
        #         mcc = node.data.dat["mcc"]
        #         for l in POT(node)
        #             if l.data.dat["mcc"] == mcc
        #                 l.data.dat["mcc"] = next_i
        #                 if isleaf(l)
        #                     mcc_map[l.label] = next_i
        #                 end
        #             end
        #         end
        #         next_i +=1 
        #         no_split_needed -=1
        #         ## else there is a split here, nothing needed
        #     else
        #         i +=1
        #     end
        # end
    end
    MCC3_new = MCC_vector_from_dict(get_MCC_as_dict(mcc_map))
    return MCC3_new
end

"""
find_split(constraint_MCC::Vector{Vector{String}}, MCC::Vector{Vector{String}})

Find the mccs of the `constraint_MCC` that are not a subset of mccs in `MCC` and calculate their intersection with mccs in `MCC`.
In order for consistency to hold these mccs must be split up. This function returns a dictionary whre each mcc that needs to 
be split is mapped to its (non empty) intersections with sets in `MCC`.
"""
function find_split(constraint_MCC::Vector{Vector{String}}, MCC2::Vector{Vector{String}})
    add_subsets = Dict()
    for mcc1 in constraint_MCC
        for mcc2 in MCC2
            if issubset(Set{String}(mcc1), Set{String}(mcc2))
                break
            elseif !isempty(intersect(Set{String}(mcc2), Set{String}(mcc1)))
                if haskey(add_subsets, mcc1)
                    append!(add_subsets[mcc1], [intersect(mcc2, mcc1)])
                else
                    add_subsets[mcc1] = [intersect(mcc2, mcc1)]
                end
            end
        end
    end
    return add_subsets
end

"""
find_set_to_split(MCC_to_split::Vector{Vector{String}}, mcc_constraint::Vector{String})

Find and return the mcc in `MCC_to_split` that contains the mcc `mcc_constraint`.
"""
function find_set_to_split(MCC_to_split::Vector{Vector{String}}, mcc_constraint::Vector{String})
    for mcc in MCC_to_split
        if issubset(Set{String}(mcc_constraint), Set{String}(mcc))
            return mcc
        end
    end
end

function find_set_to_split(MCC_to_split::Dict{Int, Vector{String}}, mcc_constraint::Vector{String})
    for (key, mcc) in MCC_to_split
        if issubset(Set{String}(mcc_constraint), Set{String}(mcc))
            return key
        end
    end
end

"""
split_mcc(mcc_to_split::Vector{String}, mcc_constraint::Vector{String}, splits::Vector{Vector{String}}, tree::Tree)

Split the `mcc_to_split` as prescribed in the `splits` vector. `mcc_constraint` is a subset of `mcc_to_split`, 
The union of all leaves in the `splits` vector should be equal to the `mcc_constraint`, `mcc_to_split` may contain more leaves.
If there are `k` splits in the `splits` vector, k-1 splits must be introduced into `mcc_to_split`. 
For each split calculate the lca, make sure the lca does not have any children that are in another split (using `tree`). 
If this not the case introduce a split "above the lca" by assigning all it's children that are in the `mcc_to_split`
to a new mcc.
"""
function split_mcc(mcc_to_split::Vector{String}, mcc_constraint::Vector{String}, splits::Vector{Vector{String}}, tree::Tree)
    new_mcc_list = Vector{String}[]
    assigned_leaves = Set{String}()
    count = length(splits) - 1
    for s in splits
        lca_ = TreeTools.lca(tree, s...)
        if Set{String}([x.label for x in POTleaves(lca_) if x.label ∈ mcc_constraint]) == Set{String}(s) && count >0
            new_mcc = [x.label for x in POTleaves(lca_) if x.label ∈ mcc_to_split]
            append!(new_mcc_list, [sort(new_mcc, lt=TreeKnit.clt)])
            union!(assigned_leaves, Set{String}(new_mcc))
            count -=1
        end
    end
    if count >0
        print("there is an issue in the split_mcc function")
    end
    append!(new_mcc_list, [sort(collect(setdiff(Set{String}(mcc_to_split), assigned_leaves)), lt=TreeKnit.clt)])
    return sort(new_mcc_list, lt=TreeKnit.clt)
end

"""
split_MCCs(first, second, third, tree1, tree2)

For each MCC triplet, e.g. MCC_ab, MCC_ac and MCC_bc the transitivity relation should hold that 
if a branch is in an mcc in MCCab and that branch in an mcc in MCCac that branch shoudl also be in an mcc in MCCbc (consistency). 
If this relation does not hold the function performs the following:

- calculate the `constraint` of MCCab and MCCac (which branches are shared in both trees defined by `join_sets(MCCab, MCCac)`
- find the mccs of the `constraint` that are not a subset of mccs in MCCbc and calculate their intersection with mccs in MCCbc `find_split`
- randomly choose to introduce splits in MCCab or MCCac (this will lead to the constraint having these splits). If MCCac is chosen look at tree c, else look at tree b
- find which mcc the mccs that need to be split are contained in in MCCac (or MCCab) `find_set_to_split`
- split these mccs `split_mcc`.
"""
function split_MCCs(first, second, third, tree1, tree2)
    constraint = TreeKnit.join_sets([first, second])
    r = rand((1, 2))
    MCC_to_split = [first, second][r]
    tree = [tree1, tree2][r]
    constraint_split_dict = find_split(constraint, third)
    if !isempty(constraint_split_dict)
        MCC_to_split_dict = TreeKnit.get_MCC_as_dict(TreeKnit.get_mcc_map(MCC_to_split))
        next_mcc = length(MCC_to_split) + 1
        for mcc_constraint in keys(constraint_split_dict)
            mcc_to_split_key = find_set_to_split(MCC_to_split_dict, mcc_constraint)
            new_mccs = split_mcc(MCC_to_split_dict[mcc_to_split_key], mcc_constraint, constraint_split_dict[mcc_constraint], tree)
            MCC_to_split_dict[mcc_to_split_key] = new_mccs[1]
            for new_m in new_mccs[2:end]
                MCC_to_split_dict[next_mcc] = new_m
                next_mcc +=1
            end
        end
        MCC_to_split = MCC_vector_from_dict(MCC_to_split_dict)
    end
    if r==1
        return MCC_to_split, second
    else
        return first, MCC_to_split
    end
end
"""
fix_consist_sets!(pair_MCCs::MCC_set, trees::Vector{Tree{MiscData}})

Fix inconsistencies in `pair_MCCs` by splitting mccs as in `split_MCCs`.
"""
function fix_consist_sets!(pair_MCCs::MCC_set, trees::Vector{Tree{MiscData}}; topo=false)
    l_t = pair_MCCs.no_trees
    rounds = 1+ (l_t-3)*3
    not_const = TreeKnit.is_degenerate(pair_MCCs)
    if topo == true
        not_topo_const = is_topologically_degenerate(pair_MCCs, trees)
    else
        not_topo_const = false
    end
    r = 0
    while (not_const ==true || not_topo_const == true) && r < rounds
        for i in 1:(l_t-1)
            for j in (i+1):l_t
                for x in 1:l_t
                    if x ∉ Set([i, j]) 
                        if topo == true
                            third = fix_topologically_degenerate( get(pair_MCCs, (i, x)), get(pair_MCCs, (j, x)), get(pair_MCCs, (j, i)), trees[i], trees[j])
                            TreeKnit.add!(pair_MCCs, third, (i, j))
                        end
                        first, second, third = get(pair_MCCs, (i, x)), get(pair_MCCs, (j, x)), get(pair_MCCs, (j, i))
                        first, second = split_MCCs(first, second, third, trees[i], trees[j])
                        TreeKnit.add!(pair_MCCs, first, (i, x))
                        TreeKnit.add!(pair_MCCs, second, (j, x))
                    end
                end
            end
        end
        not_const = TreeKnit.is_degenerate(pair_MCCs)
        r +=1
    end
    if not_const
        @info "Cannot find a consistent ARG"
    end
    return pair_MCCs
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