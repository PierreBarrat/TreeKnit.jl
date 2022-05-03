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
					@assert i == ic "MCCs $i and $ic conflict (resp. node $(n.label) and $(c.label)).
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
			error("Reached root $(n.label) of MCC $i and have not found all nodes ($(M) out of $(length(MCCs[i])))")
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

################################################################################
################################# Confidence ###################################
################################################################################

"""
	distance_likelihood_ratio(d1, d2[, L1=1, L2=1])

Given two distances between two leaves `d1` and `d2` in two segment trees, estimate whether\
it is more likely that a reassortment happened along the way or not.

## Method

Given two paths of length `d1` and `d2`, with segments of length `L1` and `L2`, return the log-ratio of:
- the probability that the two paths are of different real lengths equal to respectively `d1` and `d2`.
- the probability that the two paths are of the same real length equal to `(d1*L1 + d2*L2) / (L1 + L2)`.

The ratio quantifies the confidence in the fact that there is a reassortment along the path.

*Note*: this assumes the mutation rate to be the same in the two segments.

## Optional sequence length

Sequence length for the two segments `L1` and `L2` can be provided. Ultimately, what \
matters are the ratio `L1/(L1+L2)` and `L2/(L1+L2)`, which quantify confidence to the \
observed `d1` and `d2`.

"""
function distance_likelihood_ratio(d1::Real, d2::Real, L1=1, L2=1)
	# expected nb of mutations on each segment if the path is shared
	n1s = (d1*L1 + d2*L2) / (L1 + L2) * L1
	n2s = (d1*L1 + d2*L2) / (L1 + L2) * L2
	# expected nb of mutations if the path is not shared
	n1 = d1*L1
	n2 = d2*L2

	llk = 0.
	# log ratios of Poisson probabilities
	if n1 != 0
		llk += - n1 + n1s + n1 * log(n1/n1s)
	else
		llk += n1s
	end
	if n2 != 0
		llk += - n2 + n2s + n2 * log(n2/n2s)
	else
		llk += n2s
	end

	return llk / (L1+L2)
end

distance_likelihood_ratio(d1::Missing, d2::Real, L1=1, L2=1) = missing
distance_likelihood_ratio(d1::Real, d2::Missing, L1=1, L2=1) = missing
distance_likelihood_ratio(d1::Missing, d2::Missing, L1=1, L2=1) = 0.
distance_likelihood_ratio(d1, d2, L1=1, L2=1) = 0.


function confidence_likelihood_ratio(MCCs, t1, t2; neighbours = :leaves)
	if neighbours == :leaves
		return map(m->_confidence_likelihood_ratio_leaves(m, t1, t2), MCCs)
	elseif neighbours == :joint
		ct1, ct2 = (copy(t1), copy(t2))
		resolve!(ct1, ct2, MCCs)
		mcc_maps = [mcc_map(ct1, MCCs), mcc_map(ct2, MCCs)]
		return map(m->_confidence_likelihood_ratio_joint(m, t1, t2, mcc_maps), MCCs)
	end
end
function _confidence_likelihood_ratio_joint(MCC, t1, t2, mcc_maps, L1=1, L2=1; verbose=false)
	verbose && @info "Confidence in MCC: " MCC
	A1 = lca(t1, MCC...)
	A2 = lca(t2, MCC...)
	if A1.isroot || A2.isroot # Infinite confidence in the root ...
		return Inf
	end

	neighbours = neighbour_joint_nodes(MCC, (t1, mcc_maps[1]), (t2, mcc_maps[2]))

	L = 0.
	for n in neighbours # n is an array of labels defining a joint internal node
		n1 = lca(t1, n...)
		n2 = lca(t2, n...)
		d1 = divtime(A1, n1)
		d2 = divtime(A2, n2)
		lk = distance_likelihood_ratio(d1, d2, L1, L2)
		L += lk
		verbose && @info "Contribution to Likelihood: " n, n1.label, n2.label, lk
	end

	return L / length(neighbours)
end
function _confidence_likelihood_ratio_leaves(MCC, t1, t2, L1=1, L2=1; verbose=false)
	verbose && @info "Confidence in MCC: " MCC
	A1 = lca(t1, MCC...)
	A2 = lca(t2, MCC...)
	if A1.isroot || A2.isroot # Infinite confidence in the root ...
		return Inf
	end

	# I could also look at the first neighbour joint nodes?
	neighbours = neighbour_leaves(MCC, t1, t2)

	L = 0.
	for n in neighbours
		n1 = t1.lnodes[n]
		n2 = t2.lnodes[n]
		d1 = divtime(A1, n1)
		d2 = divtime(A2, n2)
		lk = distance_likelihood_ratio(d1, d2, L1, L2)
		L += lk
		verbose && @info "Contribution to Likelihood: " n, lk
		# verbose && println
	end

	return L / length(neighbours)
end

function neighbour_joint_nodes(MCC, treemap, treemaps...)
	# Input `treemap` is a tuple `(tree, mcc_map)`
	S = _neighbour_joint_nodes(MCC, treemap[1], treemap[2])
	for tm in treemaps
		union!(S, _neighbour_joint_nodes(MCC, tm[1], tm[2]))
	end
	return S
end
function _neighbour_joint_nodes(MCC, tree, mp)
	# !! tree must be resolved !!
	R = lca(tree, MCC...)
	A = R.anc # Ancestor of root of MCC
	list = []
	for c in Iterators.filter(!=(R), A.child)
		l = []
		get_joint_offsprings!(l, c, mp)
		append!(list, l)
	end
	return sort(unique!(list))
end
function get_joint_offsprings!(list, A, mp)
	if !isnothing(mp[A.label])
		push!(list, map(x->x.label, POTleaves(A)))
		return nothing
	else
		for c in A.child
			get_joint_offsprings!(list, c, mp)
		end
		return nothing
	end
end

function neighbour_leaves(MCC, tree)
	A = lca(tree, MCC...) # MRCA of MCC
	S = map(x->x.label, POTleaves(A.anc))
	setdiff!(S, MCC)
	return S
end
function neighbour_leaves(MCC, tree, trees...)
	S = neighbour_leaves(MCC, tree)
	for t in trees
		union!(S, neighbour_leaves(MCC, t))
	end
	return S
end
