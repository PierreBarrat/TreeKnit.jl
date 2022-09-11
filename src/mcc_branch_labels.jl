using Plots

"""
leaf_mcc_map(MCCs::Vector{Vector{String}}, get_cluster_no =true)

Returns a dictionary of which MCC each leaf is in, if `get_cluster_no = true` returns a list of 
which MCCs contain more than one leaf. 
"""
function leaf_mcc_map(MCCs::Vector{Vector{String}}; get_cluster_no =false)
    if get_cluster_no
        mcc_map = Dict{String, Int}()
        cluster_no = Int[]
        for (i,mcc) in enumerate(MCCs)
            if length(mcc)>1
                append!(cluster_no, i)
            end
            for node in mcc
                mcc_map[node] = i
            end
        end
        return mcc_map, cluster_no
    else
        mcc_map = Dict{String, Int}()
        for (i,mcc) in enumerate(MCCs)
            for node in mcc
                mcc_map[node] = i
            end
        end
        return mcc_map
    end
end

function leaf_mcc_map(MCCs::Vector{Vector{Vector{String}}})
    sets = [union([Set([m... ]) for m in mcc]...) for mcc in MCCs]
    @assert union(sets...) == sets[1] ## make sure labels are the same in all trees
    mcc_map = Dict{String, Vector{Int}}()
    for MCC in MCCs
        for (i,mcc) in enumerate(MCC)
            for node in mcc
                if haskey(mcc_map, node)
                    append!(mcc_map[node], i)
                else
                    mcc_map[node] = [i]
                end
            end
        end
    end
    return mcc_map
end


"""
assign_mccs_PR!(n::TreeNode)

Assign `mcc`s to branches (i.e. their child node) by a Pre-order traversal starting at the root node `n`.
"""
function assign_mccs_PR!(n::TreeNode, k::Int)
	
    if isroot(n)
        n.data["mcc"] = []
        for pos in 1:k
            if !isempty(n.data["child_mccs"][pos]) && length(n.data["child_mccs"][pos])==1
                append!(n.data["mcc"], pop!(n.data["child_mccs"][pos]))
            else
                append!(n.data["mcc"], [nothing])
            end
        end
    else
        n.data["mcc"] =[]
        for pos in 1:k
            if n.anc.data["mcc"][pos] in n.data["child_mccs"][pos] # parent MCC part of children -> that is the MCC
                append!(n.data["mcc"], n.anc.data["mcc"][pos])
            elseif length(n.data["child_mccs"][pos])==1  # child is an MCC
                append!(n.data["mcc"], pop!(n.data["child_mccs"][pos]))
            else # no unique child MCC and no match with parent -> not part of an MCCs
                append!(n.data["mcc"], [nothing])
            end
        end
    end

	#delete!(n.data.dat, "child_mccs")

	if !isempty(n.child)
		for c in n.child
			assign_mccs_PR!(c, k)
		end
	end


end

function assign_mccs_PR!(t::Tree, k::Int)
	assign_mccs_PR!(t.root, k)
end


"""
assign_mccs_PR!(n::TreeNode)

Assign `mcc`s to branches (i.e. their child node) by a Pre-order traversal starting at the root node `n`.
"""
function assign_mccs_PR!(n::TreeNode)
    if isroot(n)
        if length(n.data["child_mccs"])==1
            n.data["mcc"] = collect(n.data["child_mccs"])[1]
        else
            n.data["mcc"] = nothing
        end
    else
        if n.anc.data["mcc"] in n.data["child_mccs"] # parent MCC part of children -> that is the MCC
            n.data["mcc"] = n.anc.data["mcc"]
        elseif length(n.data["child_mccs"])==1  # child is an MCC
            n.data["mcc"] = collect(n.data["child_mccs"])[1]
        else # no unique child MCC and no match with parent -> not part of an MCCs
            n.data["mcc"] = nothing
        end
    end
    #delete!(n.data.dat, "child_mccs")

    if !isempty(n.child)
        for c in n.child
            assign_mccs_PR!(c)
        end
    end

end

function assign_mccs_PR!(t::Tree)
	assign_mccs_PR!(t.root)
end



"""
mark_shared_branches!(filter::Union{Nothing, Vector{Vector{String}}}, t::Vararg{Tree})

Add a `shared_branch_constraint` parameter to the tree, branches with a `shared_branch_constraint` cannot 
have a recombination event occuring on them as they connect clades that should be together according to the 
input constraints (`constraint`). The constarint should be in the form of an MCC, where if nodes are in the same 
clade this means they cannot have a recombination event happen between them.

The function proceeds by allocating each node to the MCC it should be in using the Fitch algorithm, 
then branches which are in a MCC with 2 or more nodes are marked with `shared_branch_constraint`.
"""
function mark_shared_branches!(constraint::Union{Nothing, Vector{Vector{String}}}, t::Vararg{Tree})
	
	if isnothing(constraint)
		return "here"
	end

	mcc_map, cluster_no = leaf_mcc_map(constraint, get_cluster_no =true)
	# assign MCCs to leaves
    for tree in t
	    assign_mccs!(tree, mcc_map)

		for n in POT(tree)
			if n.data["mcc"] in cluster_no 
                m = n.data["mcc"]
                if (isroot(n) || m==n.anc.data["mcc"])
				    n.data["shared_branch_constraint"] = true
                elseif any([(c!=n && (c.data["mcc"]==m ||  any([l.data["mcc"] ==m for l in POTleaves(c)]))) for c ∈ n.anc.child])
                    n.data["shared_branch_constraint"] = true
                else
                    n.data["shared_branch_constraint"] = false
                end
            else
				n.data["shared_branch_constraint"] = false
			end
		end
	end
end

function assign_mccs!(tree_list::Vector{Tree{TreeTools.MiscData}}, mcc_map::Dict{String, Int}) 
	
	# assign MCCs to leaves
	for tree in tree_list
		# assign MCCs to leaves
        for leaf in tree.lleaves
            leaf.second.data["child_mccs"] = Set([mcc_map[leaf.second.label]])
            leaf.second.data["mcc"] = mcc_map[leaf.second.label]
        end
    
        # reconstruct MCCs with Fitch algorithm
        for n in POT(tree)
            if !n.isleaf
                common_mccs = intersect([c.data["child_mccs"] for c in n.child]...)
                if !isempty(common_mccs)
                    n.data["child_mccs"] = common_mccs
                else
                    n.data["child_mccs"] = union([c.data["child_mccs"] for c in n.child]...)
                end
            end
        end
        assign_mccs_PR!(tree)

	end
end

function assign_mccs!(tree::Tree{TreeTools.MiscData}, mcc_map::Dict{String, Int}) 
    return  assign_mccs!([tree], mcc_map) 
end

function assign_all_mccs!(tree::Tree{TreeTools.MiscData}, tree_list::Vector{Tree{TreeTools.MiscData}}, 
    MCCs_dict::MCC_set)
    
    len_tree_list = length(tree_list)
    MCCs_list = Vector{Vector{String}}[]
    for t in tree_list
        append!(MCCs_list, [get(MCCs_dict,(tree.label, t.label))])
    end
    mcc_map = leaf_mcc_map(MCCs_list)
    
    assign_all_mccs!(tree, len_tree_list, mcc_map)
    
end

function assign_all_mccs!(tree::Tree{TreeTools.MiscData}, len_tree_list::Int, mcc_map::Dict{String, Vector{Int}})
    # assign MCCs to leaves
    for leaf in tree.lleaves
        leaf.second.data["child_mccs"] = [Set([mcc_map[leaf.second.label][pos]]) for pos in 1:len_tree_list]
        leaf.second.data["mcc"] = mcc_map[leaf.second.label]
    end

    # reconstruct MCCs with Fitch algorithm
    for n in POT(tree)
        if !n.isleaf
            common_mccs = [intersect([c.data["child_mccs"][pos] for c in n.child]...) for pos in 1:len_tree_list]
            n.data["child_mccs"] = Set{Int}[]
            for pos in 1:len_tree_list
                if !isempty(common_mccs[pos])
                    append!(n.data["child_mccs"], [common_mccs[pos]])
                else
                    append!(n.data["child_mccs"], [union([c.data["child_mccs"][pos] for c in n.child]...)])
                end
            end
        end
    end

	assign_mccs_PR!(tree, len_tree_list)
    
end

function get_MCC_as_dict(mcc_map::Dict{String, Int})
    MCC_dict = Dict{Int, Vector{String}}()
    for (key, item) in mcc_map
        if haskey(MCC_dict, item)
            append!(MCC_dict[item], [key])
        else
            MCC_dict[item] = [key]
        end
    end
    return MCC_dict
end

function get_MCC_as_dict(mcc_map::Dict{String, Vector{Int}}, pos::Int)
    MCC_dict = Dict{Int, Vector{String}}()
    for (key, item) in mcc_map
        if haskey(MCC_dict, item[pos])
            append!(MCC_dict[item[pos]], [key])
        else
            MCC_dict[item[pos]] = [key]
        end
    end
    return MCC_dict
end

function MCC_vector_from_dict(dict_::Dict{Int, Vector{String}})
    MCCs_new = Vector{String}[]
    for (num, mcc) in dict_
        append!(MCCs_new, [sort(mcc)])
    end
    MCCs = TreeKnit.sort(MCCs_new, lt=TreeKnit.clt)
    return MCCs
end


function check_merge(trees, to_be_merged)
    function get_masked_splitlist(tree, masked_leaves)
        l = sort(collect(keys(tree.lleaves)))
        leafmap = Dict(leaf=>i for (i,leaf) in enumerate(l))
        # Compute mask : leaves that are descendents or `r`
        mask = zeros(Bool, length(l))
        set_mask(n) = if n.label ∈ masked_leaves mask[leafmap[n.label]] = true end
        map!(set_mask, tree.root)
        return  SplitList(tree.root, l, mask, leafmap)
    end
    S3 = get_masked_splitlist(trees[2], to_be_merged)
    S2 = get_masked_splitlist(trees[1], to_be_merged)
    if all([TreeKnit.iscompatible(s, S2) for s in S3])
        return true
    else
        return false
    end
end

function fix_consist!(MCCs, trees; merge=true, split=true, i=nothing)

    @assert length(trees)==2 && length(MCCs)==3
    constraint = join_sets([MCCs[1], MCCs[2]])

    ##randomly choose which MCC to put the split in if mccs in MCC3 cannot be merged
    ##this determines which tree to iterate over when fixing inconsistencies
    if isnothing(i)
        i = rand((2, 1))
    end
    tree = copy(trees[i])

    mark_shared_branches!(constraint, tree)

    next_mcc_i = length(MCCs[i]) + 1
    next_mcc_joint = length(constraint) +1
    next_mcc_3 = length(MCCs[3]) + 1

    #get dictionary of which mcc each leaf is in
    mcc_map = leaf_mcc_map([constraint, MCCs[i], MCCs[3]])
    assign_all_mccs!(tree, 2, mcc_map)

    #get dictionary of which leaves are in a mcc
    MCCs3_dict = get_MCC_as_dict(mcc_map, 3)
    MCC_constraint_dict = get_MCC_as_dict(mcc_map, 1)

    for n in POT(tree)
        if isroot(n)
            continue
        end
        if n.data.dat["shared_branch_constraint"]==true ##should be in a MCC
            if !(TreeKnit.is_branch_in_mccs(n, MCCs3_dict))
                if merge == true
                    ##if can be merged modify MCCs[3] instead of MCCs[i]
                    mcc = n.data.dat["mcc"][1]
                    #get all nodes that should be together according to constraint
                    desired_mcc_clade = MCC_constraint_dict[mcc]

                    #get which mccs these children are in in MCCs3
                    to_be_merged_mcc_set = Set([mcc_map[c][3] for c in desired_mcc_clade])
                    #get all nodes in these mccs
                    to_be_merged = Set{String}()
                    for m in to_be_merged_mcc_set
                        union!(to_be_merged, Set(MCCs3_dict[m]))
                    end
                    #check if the mccs can be merged in MCCs3, i.e. check if splitlists of the 
                    #merged mccs are compatible in the 2 trees
                    if check_merge(trees, to_be_merged)
                        #if can be merged all the nodes in the to be merged mccs should now
                        #be in one MCC together
                        for m in to_be_merged_mcc_set
                            if haskey(MCCs3_dict, next_mcc_3)
                                append!(MCCs3_dict[next_mcc_3], MCCs3_dict[m])
                            else
                                MCCs3_dict[next_mcc_3] = MCCs3_dict[m]
                            end
                            for n in MCCs3_dict[m]
                                mcc_map[n][3] = next_mcc_3
                            end
                            delete!(MCCs3_dict, m)
                        end
                        next_mcc_3 +=1
                    end
                elseif split == true
                    ##split MCCi to make transitivity hold
                    mcc = n.data.dat["mcc"]
                    for x in POTleaves(n)
                        if mcc_map[x.label][2]==mcc[2]
                            mcc_map[x.label][2] = next_mcc_i
                        end
                    end
                    next_mcc_i +=1
                end
            end
        end
    end
    #if MCCs have changed update them from dictionary
    if next_mcc_i > length(MCCs[i]) + 1
        MCCsi_dict = get_MCC_as_dict(mcc_map, 2)
        MCCs_new = Vector{String}[]
        for (num, mcc) in MCCsi_dict
            append!(MCCs_new, [sort(mcc)])
        end
        MCCs[i] = TreeKnit.sort(MCCs_new, lt=TreeKnit.clt)
    end
    if next_mcc_3 > length(MCCs[3]) + 1
        MCCs_new = Vector{String}[]
        for (num, mcc) in MCCs3_dict
            append!(MCCs_new, [sort(mcc)])
        end
        MCCs[3] = TreeKnit.sort(MCCs_new, lt=TreeKnit.clt)
    end
    return MCCs
end

function fix_consist!(pair_MCCs::MCC_set, trees::Vector{Tree{MiscData}}; rounds=1, verbose=false, merge=false)
    l_t = pair_MCCs.no_trees
    not_const = is_degenerate(pair_MCCs)
    rep = 0
    while not_const ==true && rep <rounds
        for i in 1:(l_t-1)
            for j in (i+1):l_t
                for x in 1:l_t
                    if x ∉ Set([i, j]) 
                        first = get(pair_MCCs, (i, x))
                        second = get(pair_MCCs, (j, x))
                        third = get(pair_MCCs, (j, i))
                        new_MCCs = fix_consist!([first, second, third], [trees[i], trees[j]], merge=merge)
                        add!(pair_MCCs, new_MCCs[1], (i, x))
                        add!(pair_MCCs, new_MCCs[2], (j, x))
                        add!(pair_MCCs, new_MCCs[3], (i, j))
                    end
                end
            end
        end
        rep +=1
        not_const = is_degenerate(pair_MCCs)
    end
    if not_const
        verbose && @info "Cannot find a consistent ARG"
    end
    return pair_MCCs
end


function get_recombination_sites(first_tree::Tree{TreeTools.MiscData}, tree_list::Vector{Tree{TreeTools.MiscData}}, 
    MCCs::MCC_set)

    MCC_lists = [get(MCCs, first_tree.label, t.label) for t in tree_list]
    
    len_trees = MCCs.no_trees -1 
    recombination_pairs_list = Dict{Int, Tuple{TreeNode, TreeNode}}[]
    for pos in 1:len_trees
        tree = tree_list[pos]
        checked_mccs = Set{Int}()
        recombination_sites = Dict{Int, TreeNode}()
        for (key, n) in tree.lleaves
            if n.data["mcc"] ∉ checked_mccs
                while n!=tree.root
                    if n.data["mcc"] != n.anc.data["mcc"]
                        recombination_sites[n.data["mcc"]] = n
                        break
                    end
                    n = n.anc
                end
            end
            push!(checked_mccs, n.data["mcc"])
            if length(checked_mccs) == length(MCC_lists[pos])
                break
            end
        end
        checked_mccs_first = Set{Int}()
        recombination_sites_first = Dict{Int, TreeNode}()
        for (key, n) in first_tree.lleaves
            if n.data["mcc"][pos] ∉ checked_mccs_first
                while n!=first_tree.root
                    if n.data["mcc"][pos] != n.anc.data["mcc"][pos]
                        recombination_sites_first[n.data["mcc"][pos]] = n
                        break
                    end
                    n = n.anc
                end
            end
            push!(checked_mccs, n.data["mcc"][pos])
            if length(checked_mccs_first) == length(MCC_lists[pos])
                break
            end
        end
        recombination_pairs = Dict{Int, Tuple{TreeNode, TreeNode}}()
        for mcc in keys(recombination_sites_first)
            if mcc in keys(recombination_sites)
                recombination_pairs[mcc] = (recombination_sites_first[mcc], recombination_sites[mcc])
            end
        append!(recombination_pairs_list, [recombination_pairs])
        end
    end
    return recombination_pairs_list
end
