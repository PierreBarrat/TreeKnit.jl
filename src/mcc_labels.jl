

"""
get_mcc_map(MCCs::Vector{Vector{String}}, get_cluster_no =true)

Returns a dictionary of which MCC each leaf is in, if `get_cluster_no = true` returns a list of 
which MCCs contain more than one leaf. 
"""
function get_mcc_map(MCCs::Vector{Vector{String}}; get_cluster_no =true)
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


"""
PRT!(n::TreeNode)

Assign `mcc`s to branches (i.e. their child node) by a Pre-order traversal starting at the root node `n`.
"""
function PRT!(n::TreeNode)
	
	if isroot(n)
		if !isempty(n.data["child_mccs"])
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

		for n in POT(tree)
			if n.data["mcc"] in cluster_no && (isroot(n) || n.data["mcc"]== n.anc.data["mcc"])
				n.data["mask"] = true
			else
				n.data["mask"] = false
			end
		end
	end
end