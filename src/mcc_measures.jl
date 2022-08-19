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
    MCC degeneracy measure, the 
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

function fix_consist_sets!(pair_MCCs::MCC_set, trees::Vector{Tree{MiscData}}; rounds = 1)
    l_t = pair_MCCs.no_trees
    rounds = 1+ (l_t-3)*3
    not_const = TreeKnit.is_degenerate(pair_MCCs)
    r = 0
    while not_const ==true && r < rounds
        for i in 1:(l_t-1)
            for j in (i+1):l_t
                for x in 1:l_t
                    if x ∉ Set([i, j]) 
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