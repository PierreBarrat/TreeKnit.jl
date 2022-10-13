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