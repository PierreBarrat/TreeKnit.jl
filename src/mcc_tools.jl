"""
	map_mccs(tree, MCCs; internals = true)
	map_mccs(MCCs)

Create a map from nodes of `tree` to the mcc they belong to.
MCCs are identified using their index in the array `MCCs`.
Nodes outside of any mcc are mapped to `nothing`.
The output is a `Dict`.
Use `internals=false` or the second form to map leaves only.
"""
map_mccs(MCCs) = map_mccs_leaves(MCCs)
function map_mccs(tree, MCCs; internals = true)
	leaf_mcc_map = map_mccs_leaves(MCCs)
	if internals
		fitch_up = map_mccs_fitch_up(tree, leaf_mcc_map)
		#@debug "Upward pass of Fitch" ficth_up=fitch_up
		fitch_down = map_mccs_fitch_down(tree, fitch_up)
		return fitch_down
	else
		return leaf_mcc_map
	end
end

"""
	map_mccs!(tree::Tree{TreeTools.MiscData}, MCCs; internals = true)

Same as `map_mccs`, but stores information in the `data` field of tree nodes.
A key `child_mccs` is also added to `data`, for use in some internal functions.
"""
function map_mccs!(tree::Tree{TreeTools.MiscData}, MCCs; internals = true)
	leaf_mcc_map = map_mccs_leaves(MCCs)

	if internals
		# Compute the map and add it to nodes
		fitch_up = map_mccs_fitch_up(tree, leaf_mcc_map)
		#@debug "Upward pass of Fitch" ficth_up=fitch_up
		fitch_down = map_mccs_fitch_down(tree, fitch_up)
		for (label, mcc) in fitch_down
			tree[label].data["mcc"] = mcc
			tree[label].data["child_mccs"] = fitch_up[label]
		end
	else
		# Just use the leaf_map
		for (label, mcc) in leaf_mcc_map
			tree[label].data["mcc"] = mcc
			tree[label].data["child_mccs"] = Set([mcc])
		end
	end

	return nothing
end
function map_mccs!(tree, MCCs; internals=true)
	@error "`map_mccs!` is only for trees with data attached to nodes, *i.e.* `Tree{TreeTools.MiscData}`.
	 Try `map_mccs` to get a map from nodes to mcc.
	"
	error("Incorrect method")
end

"""
	map_mccs_leaves(MCCs::Vector{Vector{<:AbstractString}})

Returns a dictionary of which MCC each leaf is in.
MCCs are identified by their index in `MCCs`.
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
			# coherent messages from all children - part of an MCC
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
	get_leaves_order(tree, MCCs)

Return an ordering of leaves in MCCs given a tree.

Output `order::Vector`: `order[i][label] == n` means that `label` appears at position `n` in MCC `i` for the input tree.

"""
function get_leaves_order(tree, MCCs)
	order = [Dict() for i in 1:length(MCCs)]
	for (i,mcc) in enumerate(MCCs)
		mcc_order = Dict{String, Int}()
		_get_leaves_order!(mcc_order, tree.root, mcc)
		order[i] = mcc_order
	end

	return order
end

function _get_leaves_order!(order::Dict, n::TreeNode, mcc::Vector{<:AbstractString})
	if n.isleaf && in(n.label, mcc)
		order[n.label] = length(order) + 1
	else
		for c in n.child
			_get_leaves_order!(order, c, mcc)
		end
	end

	return nothing

end

"""
	sort_polytomies!(t1::Tree, t2::Tree, MCCs)

Sort nodes of `t2` such that leaves of `t1` and `t2` in the same MCC face each other.
In a tanglegram with nodes colored according to MCC, lines of the same color will not cross.
The order of `t1` serves as a guide and is left unchanged.
"""
function sort_polytomies!(t1::Tree{T}, t2::Tree{T}, MCCs) where T
	mcc_order = get_leaves_order(t1, MCCs)
	_mcc_map = map_mccs(t2, MCCs)
	sort_polytomies!(t2.root, MCCs, _mcc_map, mcc_order)
	node2tree!(t2, t2.root)
	return nothing
end

function sort_polytomies!(n::TreeNode, MCCs, mcc_map, mcc_order)
	if n.isleaf
		i = mcc_map[n.label] # in MCCs[i]
		return mcc_order[i][n.label], 1 # rank in the polytomy, number of nodes in clade
	else
		i = mcc_map[n.label] # Id of the current MCC (can be `nothing`)

		# Construct ranking or children
		# positive if they are in same MCC as current node, negative otherwise
		# r = 1
		rank = Array{Any}(undef, length(n.child))
		for (k, c) in enumerate(n.child)
			if mcc_map[c.label] == i
				# child is in mcc. we know its rank in the MCC
				rank[k] = sort_polytomies!(c, MCCs, mcc_map, mcc_order)
			else
				# child is not in mcc: only count its ladder rank
				_trash, rladder = sort_polytomies!(c, MCCs, mcc_map, mcc_order)
				rank[k] = (nothing, rladder)
				# r += 1
			end
		end

		# # Sort children
		_sort_children!(n, rank)

		# # Return minimum value of the rank of nodes in the MCC.
		# # Will be used by parent to rank the current node
		return findmin(r -> isnothing(r[1]) ? Inf : r[1], rank)[1], sum(x[2] for x in rank)
	end
end

function _sort_children!(n, rank)
	# Custom sorting function
	## if the nodes are in the same MCC, rank them according to that
	## otherwise, rank according to ladder rank (biggest clade at the beginning)
	function _isless(x,y)
		if !isnothing(x[1]) && !isnothing(y[1])
			return x[1] < y[1] # those two numbers can't be equal (rank in mcc)
		else
			return x[2] < y[2]
		end
	end

	children_order = sortperm(rank, lt = _isless)
	n.child = n.child[children_order]
end



## BEFORE REMOVING `write_mccs!`:
## Check whether it's used in annotating auspice json files
# """
#     write_mccs!(trees::Dict, MCCs::Dict, key=:mcc_id)

# Write MCCs id to field `data.dat[key]` of tree nodes. Expect `trees` indexed by single segments, and `MCCs` indexed by pairs of segments.
# """
# function write_mccs!(trees::Dict, MCCs::Dict, key=:mcc_id; overwrite=false)
#     for ((i,j), mccs) in MCCs
#         k = Symbol(key,"_$(i)_$(j)")
#         write_mccs!(trees[i], mccs, k, overwrite=overwrite)
#         write_mccs!(trees[j], mccs, k, overwrite=overwrite)
#     end
# end
# """
#     write_mccs!(t::Tree, MCCs, key=:mcc_id)

# Write MCCs id to field `data.dat[key]` of tree nodes.
# """
# function write_mccs!(t::Tree{TreeTools.MiscData}, MCCs, key=:mcc_id; overwrite=false)
#     for (i,mcc) in enumerate(MCCs)
#         for label in mcc
#             t.lleaves[label].data.dat[key] = i
#         end
#         for n in Iterators.filter(n->!n.isleaf, values(t.lnodes))
#             if is_branch_in_mcc(n, mcc)
#                 if !overwrite && haskey(n.data.dat, key)
#                     error("Node $(n.label) already has an MCC attributed")
#                 end
#                 n.data.dat[key] = i
#             end
#         end
#     end
#     nothing
# end

"""
    is_branch_in_mccs(n::TreeNode, mccs::Array)


Is the branch from `n` to `n.anc` in an element of `mccs`?
"""
function is_branch_in_mccs(n::TreeNode, mccs::Dict)
    for mcc in values(mccs)
        if is_branch_in_mcc(n, mcc)
            return true
        end
    end
    return false
end
function is_branch_in_mccs(n::TreeNode, mccs)
    for mcc in mccs
        if is_branch_in_mcc(n, mcc)
            return true
        end
    end
    return false
end

"""
    is_branch_in_mcc(n::TreeNode, mcc::Array{<:AbstractString})

Is the branch from `n` to `n.anc` in `mcc`?
The clade defined by `n` has to intersect with `mcc`, and this intersection should be strictly smaller `mcc`.
"""
function is_branch_in_mcc(n::TreeNode, mcc::Array{<:AbstractString,1})
    # Simple check
    if n.isleaf
        return length(mcc) > 1 && in(n.label, mcc)
    end

    i = count(c -> c.isleaf && in(c.label, mcc), n)

    return (i > 0 && i < length(mcc))
end
"""
    find_mcc_with_branch(n::TreeNode, mccs::Dict)

Find the mcc to which the branch from `n` to `n.anc` belongs. If `mccs` is an array, return the pair `(index, value)`. If it is a dictionary, return the pair `(key, value)`. If no such mcc exists, return `nothing`.
Based on the same idea that `is_branch_in_mcc`.
"""
function find_mcc_with_branch(n::TreeNode, mccs::Dict)
    cl = [x.label for x in POTleaves(n)]
    for (key,mcc) in mccs
        if !isempty(intersect(cl, mcc)) && !isempty(setdiff(mcc, intersect(cl, mcc)))
            return (key, mcc)
        end
    end
    return nothing
end
function find_mcc_with_branch(n::TreeNode, mccs::Array)
    cl = [x.label for x in POTleaves(n)]
    for (i,mcc) in enumerate(mccs)
        if !isempty(intersect(cl, mcc)) && !isempty(setdiff(mcc, intersect(cl, mcc)))
            return (i,mcc)
        end
    end
    return nothing
end

"""
    is_linked_pair(n1, n2, mccs)
    is_linked_pair(n1::T, n2::T, mccs::Dict{Any,Array{T,1}}) where T
    is_linked_pair(n1::T, n2::T, mccs::Array{Array{T,1},1}) where T

Can I join `n1` and `n2` through common branches only? Equivalent to: is there an `m` in `mccs` such that `in(n1,m) && in(n2,m)`?
"""
function is_linked_pair(n1::T, n2::T, mccs::Dict{Any,Array{T,1}}) where T
    for mcc in values(mccs)
        if in(n1, mcc)
            return in(n2, mcc)
        elseif in(n2, mcc)
            return false
        end
    end
    return false
end
function is_linked_pair(n1, n2, mccs)
    for mcc in values(mccs)
        if in(n1, mcc)
            return in(n2, mcc)
        elseif in(n2, mcc)
            return false
        end
    end
    return false
end
function is_linked_pair(n1::T, n2::T, mccs::Array{Array{T,1},1}) where T
    for mcc in mccs
        if in(n1, mcc)
            return in(n2, mcc)
        elseif in(n2, mcc)
            return false
        end
    end
    return false
end

"""
    find_mcc_with_node(n::String, mccs::Array{Array{<:AbstractString,1},1})

Find MCC to which `n` belongs.
"""
function find_mcc_with_node(n::String, mccs)
    for m in mccs
        if in(n, m)
            return m
        end
    end
    return nothing
end
find_mcc_with_node(n::TreeNode, mccs) = find_mcc_with_node(n.label, mccs)
find_mcc_with_node(n, mccs) = find_mcc_with_node(n.label, mccs)

