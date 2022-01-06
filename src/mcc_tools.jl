"""
	get_leave_order(tree, MCCs)

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
	mcc_map(tree::Tree, MCCS)

Compute a map between nodes of `tree` and `MCCs`.
Output `M`: `M[label] == i` means node `label` (internal or leaf) belongs to MCC number `i`
  if `typeof(i)::Int` or does not belong to an MCC if `isnothing(i)`.

## Important
`tree` has to be resolved using `MCCs` for this to work.
"""
function mcc_map(tree::Tree, MCCs)
	mcc_map = Dict()
	mcc_map!(mcc_map, tree.root, MCCs)
	return mcc_map
end
#=
Strategy: message passing going up the tree
Each node `n` returns to parent  `(i, cnt)` where:
- `i` is the MCC to which `n` belongs
- `cnt` is the number of nodes of `MCCs[i]` in the offsprings of `n`, including `n` itself
  *e.g.*: for a leaf, `cnt` is 1; for an internal, it can be higher.
This method allows us to know when we've found all nodes of a given MCC
=#
function mcc_map!(mcc_map::Dict, n::TreeNode, MCCs)
	if n.isleaf
		i = which_mcc(n.label, MCCs)
		@assert !isnothing(i) "Leaf $(n.label) does not belong to any MCC."
		mcc_map[n.label] = i
		if length(MCCs[i]) == 1
			return (nothing, 0) # we're done with MCC i
		else
			return (i, 1) # in MCC i, representing one node
		end
	else
		i = nothing # in which mcc is `n`
		M = 0 # how many nodes of this MCC does it has as offsprings
		for c in n.child
			ic, cnt = mcc_map!(mcc_map, c, MCCs)
			if !isnothing(ic)
				if !isnothing(i)
					# if !isnothing(ic), we're not done with MCC[ic]
					# Meaning n has to belong to it
					@assert i == ic "MCCs $i and $ic conflict \
					(resp. node $(n.label) and $(c.label)). \
					Are trees resolved using the MCCs?"
				end
				i = ic
				M += cnt
			end
		end
		mcc_map[n.label] = i

		if isnothing(i) || length(MCCs[i]) == M
			return (nothing, 0) # we're done with mcc i
		elseif length(MCCs[i]) > M
			return (i, M) # in MCC i, representing M nodes
		else
			error("Reached root $(n.label) of MCC $i and have not found all nodes \
				($(M) out of $(length(MCCs[i])))")
		end
	end
end

"""
	which_mcc(n, MCCs)

Return the index of the MCC to which `n` belongs.
If `i` is returned, `in(n, MCCs[i]`) is true.
"""
function which_mcc(n::AbstractString, MCCs)
	for (i, mcc) in enumerate(MCCs)
		if in(n, mcc)
			return i
		end
	end
	return nothing
end

function sort_polytomies!(t1::Tree, t2::Tree, MCCs)
	mcc_order = get_leaves_order(t1, MCCs)
	_mcc_map = mcc_map(t2, MCCs)
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

