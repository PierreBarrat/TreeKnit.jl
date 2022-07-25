using Plots

"""
get_mcc_map(MCCs::Vector{Vector{String}}, get_cluster_no =true)

Returns a dictionary of which MCC each leaf is in, if `get_cluster_no = true` returns a list of 
which MCCs contain more than one leaf. 
"""
function get_mcc_map(MCCs::Vector{Vector{String}}; get_cluster_no =false)
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

function get_mcc_map(MCCs::Vector{Vector{Vector{String}}})
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
PRT!(n::TreeNode)

Assign `mcc`s to branches (i.e. their child node) by a Pre-order traversal starting at the root node `n`.
"""
function PRT!(n::TreeNode, k::Int)
	
    if isroot(n)
        n.data["mcc"] = []
        for pos in 1:k
            if !isempty(n.data["child_mccs"][pos])
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

	delete!(n.data.dat, "child_mccs")

	if !isempty(n.child)
		for c in n.child
			PRT!(c, k)
		end
	end


end

function PRT!(t::Tree, k::Int)
	PRT!(t.root, k)
end

function PRT!(n::TreeNode)
    if isroot(n)
        if !isempty(n.data["child_mccs"]) && length(n.data["child_mccs"])==1
            n.data["mcc"] = pop!(n.data["child_mccs"])
        else
            n.data["mcc"] = nothing
        end
    else
        if n.anc.data["mcc"] in n.data["child_mccs"] # parent MCC part of children -> that is the MCC
            n.data["mcc"] = n.anc.data["mcc"]
        elseif length(n.data["child_mccs"])==1  # child is an MCC
            n.data["mcc"] = pop!(n.data["child_mccs"])
        else # no unique child MCC and no match with parent -> not part of an MCCs
            n.data["mcc"] = nothing
        end
    end
    delete!(n.data.dat, "child_mccs")

    if !isempty(n.child)
        for c in n.child
            PRT!(c)
        end
    end

end

function PRT!(t::Tree)
	PRT!(t.root)
end

"""
add_mask!(filter::Union{Nothing, Vector{Vector{String}}}, t::Vararg{Tree})

Add a `mask` parameter to the tree, branches with a `mask` cannot have a recombination event occuring 
on them as they connect clades that should be together according to the input constraints (`filter`).
The filter should be in the form of an MCC, where if nodes are in the same clade this means they cannot
have a recombination event happen between them.

The function proceeds by allocating each node to the MCC it should be in using the Fitch algorithm, 
then branches which are in a MCC with 2 or more nodes are marked with `mask`.
"""
function add_mask!(filter::Union{Nothing, Vector{Vector{String}}}, t::Vararg{Tree})
	
	if isnothing(filter)
		return []
	end
	
	mcc_map, cluster_no = get_mcc_map(filter, get_cluster_no =true)
	# assign MCCs to leaves
    for tree in t
	    assign_mccs!(mcc_map, tree)

		for n in POT(tree)
			if n.data["mcc"] in cluster_no && (isroot(n) || n.data["mcc"]== n.anc.data["mcc"])
				n.data["mask"] = true
			else
				n.data["mask"] = false
			end
		end
	end
end

function assign_mccs!(mcc_map::Dict{String, Int}, t::Vector{Tree{TreeTools.MiscData}}) 
	
	# assign MCCs to leaves
	for tree in t
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

		PRT!(tree)

	end
end

function assign_mccs!(mcc_map::Dict{String, Int}, t::Tree{TreeTools.MiscData}) 
    return  assign_mccs!(mcc_map, [t]) 
end

function assign_all_mccs!(tree::Tree{TreeTools.MiscData}, tree_list::Vector{Tree{TreeTools.MiscData}}, 
    MCCs_dict::Dict{Set{String}, Vector{Vector{String}}})
    
    len_tree_list = length(tree_list)
    MCCs_list = Vector{Vector{String}}[]
    for t in tree_list
        append!(MCCs_list, [MCCs_dict[Set([tree.label, t.label])]])
    end
    mcc_map = get_mcc_map(MCCs_list)
    
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

	PRT!(tree, len_tree_list)
    
end

function get_recombination_sites(first_tree::Tree{TreeTools.MiscData}, tree_list::Vector{Tree{TreeTools.MiscData}}, 
    MCCs::Union{Vector{Vector{Vector{String}}}, Dict{Set{String}, Vector{Vector{String}}}})

    if isa(MCCs, Dict)
        MCC_lists = Vector{Vector{String}}[]
        for t in tree_list
            append!(MCC_lists,  [MCCs[Set([first_tree.label, t.label])]])
        end
    else
        MCC_lists = MCCs
    end
    
    len_trees = length(tree_list)
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
