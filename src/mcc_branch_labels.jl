function map_mccs(tree, MCCs; internals = true)
	leaf_mcc_map = map_mccs_leaves(MCCs)
	if internals
		fitch_up = map_mccs_fitch_up(tree, leaf_mcc_map)
		@debug "Upward pass of Fitch" ficth_up=fitch_up
		fitch_down = map_mccs_fitch_down(tree, fitch_up)
		return fitch_down
	else
		return leaf_mcc_map
	end
end

map_mccs(MCCs) = map_mccs_leaves(MCCs)

function map_mccs!(tree::Tree{TreeTools.MiscData}, MCCs; internals = true)
	mcc_map = map_mccs(tree, MCCs; internals)
	for (label, mcc) in mcc_map
		tree[label].data["mcc"] = mcc
	end

	return mcc_map
end
function map_mccs!(tree, MCCs; internals=true)
	@error "`map_mccs!` is only for trees with data attached to nodes,\
	 *i.e.* `Tree{TreeTools.MiscData}`.
	 Try `map_mccs` to get a map from nodes to mcc.
	"
	error("Incorrect method")
end

"""
	map_mccs_leaves(MCCs::Vector{Vector{<:AbstractString}})

Returns a dictionary of which MCC each leaf is in. MCCs are identified by their index.
"""
function map_mccs_leaves(MCCs::Vector{<:Vector{<:AbstractString}})
    leaf_mcc_map = Dict{String, Union{Int, Nothing}}() # Union to match with `map_mccs`
    for (i,mcc) in enumerate(MCCs)
        for node in mcc
            leaf_mcc_map[node] = i
        end
    end
    return leaf_mcc_map
end

function map_mccs_fitch_up(tree, leaf_mcc_map)
	fitch_up = Dict{String, Any}()

	# Assign all leaves
	for leaf in leaves(tree)
		fitch_up[leaf.label] = Set([leaf_mcc_map[leaf.label]])
	end

	# Go up the tree
	for n in Iterators.filter(!isleaf, POT(tree))
		if isroot(n)
			fitch_up[n.label] = intersect([fitch_up[c.label] for c in n.child]...)
		else
			common_mccs = intersect([fitch_up[c.label] for c in n.child]...)
			fitch_up[n.label] = if !isempty(common_mccs)
				common_mccs
			else
				union([fitch_up[c.label] for c in n.child]...)
			end
		end
	end

	return fitch_up
end

function map_mccs_fitch_down(tree, fitch_up)
	fitch_down = Dict{String, Union{Int, Nothing}}()
	_map_mccs_fitch_down!(fitch_down, tree.root, fitch_up)
	return fitch_down
end
function _map_mccs_fitch_down!(fitch_down, n, fitch_up)
	if isroot(n)
		# Root is in an MCC if it gets coherent messages from all children
		fitch_down[n.label] = if length(fitch_up[n.label]) == 1
			first(fitch_up[n.label])
		else
			nothing
		end
	else
		if length(fitch_up[n.label]) == 1
			# coherent messages from all children - leaves are treated here
			fitch_down[n.label] = first(fitch_up[n.label])
		elseif fitch_down[n.anc.label] in fitch_up[n.label]
			# Non empty intersection between ancestor MCC and messages from children
			fitch_down[n.label] = fitch_down[n.anc.label]
		else
			# n is not in an MCC ...
			fitch_down[n.label] = nothing
		end
	end

	for c in n.child
		_map_mccs_fitch_down!(fitch_down, c, fitch_up)
	end

	return nothing
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
	assign_mccs!(t::Vector{Tree{TreeTools.MiscData}}, mcc_map::Dict{String, Int})

Assign each node, (leaf and internal) node to the MCC that they are a part, takes dictionary `mcc_map` 
(leaf => MCC) as input and a tree `t`, if there is a conflict or it is unknown which MCC a node is part of 
this is labeled as None, this is a distinction to `mcc_map` which will throw an error if there is a conflict
between nodes, meaning that it cannot be called on unresolved trees. 

"""
function assign_mccs!(tree::Tree{TreeTools.MiscData}, mcc_map::Dict{String, Int}) 
	
	# assign MCCs to leaves
    for leaf in tree.lleaves
        leaf.second.data["child_mccs"] = Set([mcc_map[leaf.second.label]])
        leaf.second.data["mcc"] = mcc_map[leaf.second.label]
    end

    # reconstruct MCCs with Fitch algorithm
    for n in POT(tree)
        if isroot(n)
            n.data["child_mccs"] = intersect([c.data["child_mccs"] for c in n.child]...)
        else
            if !n.isleaf
                common_mccs = intersect([c.data["child_mccs"] for c in n.child]...)
                if !isempty(common_mccs)
                    n.data["child_mccs"] = common_mccs
                else
                    n.data["child_mccs"] = union([c.data["child_mccs"] for c in n.child]...)
                end
            end
        end
    end

    assign_mccs_PR!(tree)
end
