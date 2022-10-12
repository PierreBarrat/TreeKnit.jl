"""
	arg_from_trees(t1, t2, MCCs)

Construct an `ARG` from two trees and a set of MCCs.
"""
function arg_from_trees(it1, it2, MCCs)
	if !ismissing(TreeTools.branch_length(it1.root)) || !ismissing(TreeTools.branch_length(it2.root))
		@warn "Root node of one of the input trees has a non missing branch length.
		This may cause issues when setting the branch length of the ARG"
	end
	#
	t1 = copy(it1)
	t2 = copy(it2)
	# Resolve trees
	resolve!(t1, t2, MCCs)
	# Find shared nodes in both trees
	X1, X2 = shared_nodes(t1, t2, MCCs)
	# Introduce shared singletons of one into the other
	fix_shared_singletons!(t1, t2, X1, X2, MCCs)
	# Make an initial ARG with the first tree
	arg, lm1 = arg_from_tree(t1, X1, 1)
	# Add the second tree
	lm2 = add_tree_to_arg!(arg, t2, X2, lm1, 2)
	# Point their root to the same node
	# set_global_root!(arg)
	# Reverse label map
	rlm = reverse_label_map(arg, lm1, lm2)
	# Deal with branch lengths of the ARG
	set_branch_length!(arg, t1, t2, rlm)

	return arg, rlm, lm1, lm2
end

"""
	arg_from_tree(t::Tree, SN::Dict, clr=1)

Build an initial `ARG` using one tree.
"""
function arg_from_tree(t::Tree, SN::Dict, clr=1)
	label_map = Dict{String, Any}()
	arg = ARG(Dict(), Dict(), Dict(), Dict())

	# Intialize root
	a_R = ARGNode()
	add_color!(a_R, clr)
	set_root!(a_R, clr)
	t.root.isleaf && set_label!(a_R, t.root.label)
	label_map[t.root.label] = a_R.label
	arg.nodes[a_R.label] = a_R
	arg.roots[clr] = a_R

	# Go down in tree
	for t_c in t.root.child
		grow_arg_from_treenode!(arg, a_R, t_c, SN, label_map, clr)
	end

	return arg, label_map
end

"""
	grow_arg_from_treenode!(
		arg::ARG,
		a_a::ARGNode,
		t_n::TreeTools.TreeNode,
		SN::Dict,
		label_map::Dict,
		clr,
	)

Build an initial `ARG` using one tree.
Create a new `a_n::ARGNode` corresponding to `t_n::TreeNode`, and graft it onto the current
tip of the ARG `a_a`.
If a reassortment occurs above `t_n`, a hybrid node `a_h` is created and the grafting is done
such that we have `a_a` --> `a_h` --> `a_n`.

## Inputs
- `SN::Dict` is a dictionary giving the shared/non-shared status of nodes in the trees.
  It is obtained with the function `shared_nodes`.
- `label_map::Dict`: Pointing from `TreeNode` to `ARGNode`.
"""
function grow_arg_from_treenode!(
	arg::ARG,
	a_a::AbstractARGNode,
	t_n::TreeTools.TreeNode,
	SN::Dict,
	label_map::Dict,
	clr,
)
	treenodetype = SN[t_n.label][1] # Is t_n shared / mcc_root / non shared
	root_in_other_tree = SN[t_n.label][3] # Is t_n the root in the other tree

	# Creating an ARGNode to represent t_n
	## If t_n is the root of an MCC (and not the root of the other tree)
	## Create a hybrid node a_h, a node a_n for t_n, and do a_a > a_h > a_n
	## Otherwise, create a node a_n for t_n, and do a_a > a_n
	a_n, a_h = if treenodetype == :mcc_root && !root_in_other_tree
		a_h = ARGNode(; hybrid=true)
		a_n = ARGNode()
		add_color!(a_h, clr)
		add_color!(a_n, clr)
		graft!(a_a, a_h, clr)
		graft!(a_h, a_n, clr)

		(a_n, a_h)
	else
		a_n = ARGNode()
		add_color!(a_n, clr)
		graft!(a_a, a_n, clr)

		(a_n, nothing)
	end
	t_n.isleaf && set_label!(a_n, t_n.label)

	# Adding new nodes to `arg`
	arg.nodes[label(a_n)] = a_n
	if !isnothing(a_h)
		arg.nodes[label(a_h)] = a_h
		arg.hybrids[label(a_h)] = a_h
	end
	t_n.isleaf && (arg.leaves[label(a_n)] = a_n)

	# Adding new node to label_map
	label_map[t_n.label] = label(a_n)

	# Recursing down the tree
	for t_c in t_n.child
		grow_arg_from_treenode!(
			arg, a_n, t_c, SN, label_map, clr
		)
	end

	return nothing
end

"""
	add_tree_to_arg!(arg, t, SN, lmref, clr = 2)

Add tree to an existing `ARG`, using color `clr`.
This is only designed for an ARG of two trees.
"""
function add_tree_to_arg!(arg, t, SN, lmref, clr = 2)
	lm = Dict()

	# Deal with root of t
	a_R = if SN[t.root.label][1] == :non_shared
		ARGNode()
	else
		n2 = SN[t.root.label][2] # Label of t.root in ref tree
		arg.nodes[lmref[n2]]
	end
	add_color!(a_R, clr)
	set_root!(a_R, clr)

	arg.roots[clr] = a_R
	arg.nodes[label(a_R)] = a_R
	lm[t.root.label] = a_R.label

	# Go down t
	for t_c in t.root.child
		add_treenode_to_arg!(arg, a_R, t_c, SN, lm, lmref, clr)
	end

	return lm
end

function add_treenode_to_arg!(
	arg::ARG,
	a_a::AbstractARGNode,
	t_n::TreeTools.TreeNode,
	SN::Dict,
	lm::Dict,
	lmref::Dict,
	clr
)
	treenodetype = SN[t_n.label][1]
	root_in_other_tree = SN[t_n.label][3]

	a_n = if treenodetype == :non_shared
		# New node : we create a new argnode a_n
		a_n = ARGNode()
		add_color!(a_n, clr)
		graft!(a_a, a_n, clr)
		arg.nodes[label(a_n)] = a_n
		a_n
	elseif treenodetype == :mcc_root
		# t_n is either (i) the root of the other tree, or (ii) the point of a reassortment
		# In either case, there is already a node a_n representing it
		n2 = SN[t_n.label][2]
		a_n = arg.nodes[lmref[n2]]
		add_color!(a_n, clr)
		if root_in_other_tree
			# (i)
			graft!(a_a, a_n, clr)
		else
			# (ii): Find the corresponding hybrid node
			a_h = ancestor(a_n, othercolor(clr))
			@assert ishybrid(a_h) "$label(a_h)"
			add_color!(a_h, clr)
			graft!(a_a, a_h, clr)
			graft!(a_h, a_n, clr)
		end
		a_n
	else
		# t_n is shared and not the root of an MCC
		# --> there is already an ARGNode corresponding to it
		n2 = SN[t_n.label][2]
		a_n = arg.nodes[lmref[n2]]
		add_color!(a_n, clr)
		@assert ancestor(a_n, othercolor(clr)) == a_a
		graft!(a_a, a_n, clr)
		a_n
	end

	# Adding to label map
	lm[t_n.label] = label(a_n)

	# Recursing down the tree
	for t_c in t_n.child
		add_treenode_to_arg!(
			arg, a_n, t_c, SN, lm, lmref, clr
		)
	end

	return nothing
end

"""
Map from `ARGNode` to `(TreeNode, TreeNode)`.
"""
function reverse_label_map(arg, lm1, lm2)
	rlm = Dict()
	for an in nodes(arg)
		n1 = findfirst(==(label(an)), lm1)
		n2 = findfirst(==(label(an)), lm2)
		rlm[label(an)] = (n1, n2)
		# Some checks
		if isnothing(n1) && isnothing(n2)
			@assert ishybrid(an)
		elseif !isnothing(n1) && isnothing(n2)
			@assert color(an) == [true, false]
		elseif isnothing(n1) && !isnothing(n2)
			@assert color(an) == [false, true]
		else
			@assert color(an) == [true, true]
		end
	end

	return rlm
end
"""
	set_branch_length!(arg, t1, t2, X1, X2)

Set branch lengths of `arg` following a simple heuristic.
For a given node `a` of the ARG:
- If `a` and the branch above it are fully shared, average of tree branch lengths
- If `a` belongs only to one tree, then use this tree for the branch length
- If `a` is the root of an MCC and `a_h` the reassortment above it, then let `τ` be the minimum branch length above `a` in either of the trees.
  Place `a_h` at `τ/2` above `a`, and use individual trees for the length above `a_h`.

### Note:
Due to missing values propagating, we have to deal with them carefuly.
"""
function set_branch_length!(arg::ARG, t1, t2, lm)
	v1 = Dict()
	v2 = Dict()
	for a in nodes(arg)
		if !ishybrid(a)
			n1, n2  = lm[label(a)]
			if isnothing(n1)
				# Use n2 only
				set_branch_length!(a, t2.lnodes[n2].tau, 2)
			elseif isnothing(n2)
				# Use n1 only
				set_branch_length!(a, t1.lnodes[n1].tau, 1)
			else
				# Shared by n1 and n2
				# Is the node above a hybrid?
				if ishybrid(ancestor(a, 1))
					# MCC root
					a_h = ancestor(a, 1)
					τ = if ismissing(t1.lnodes[n1].tau)
						t2.lnodes[n2].tau
					elseif ismissing(t2.lnodes[n2].tau)
						t1.lnodes[n1].tau
					else
						min(t1.lnodes[n1].tau, t2.lnodes[n2].tau)
					end
					aτ = τ/2
					hτ1 = t1.lnodes[n1].tau - aτ
					hτ2 = t2.lnodes[n2].tau - aτ
					set_branch_length!(a, aτ, 1, 2)
					set_branch_length!(a_h, hτ1, 1)
					set_branch_length!(a_h, hτ2, 2)
				else
					# Fully shared node - if root of one of the trees, deal with missing tau
					τ1, τ2 = if ismissing(t1.lnodes[n1].tau)
						if t1.lnodes[n1].isroot
							(missing, t2.lnodes[n2].tau)
						else
							(t2.lnodes[n2].tau, t2.lnodes[n2].tau)
						end
					elseif ismissing(t2.lnodes[n2].tau)
						if t2.lnodes[n2].isroot
							(t1.lnodes[n1].tau, missing)
						else
							(t1.lnodes[n1].tau, t1.lnodes[n1].tau)
						end
					else
						τ = (t1.lnodes[n1].tau + t2.lnodes[n2].tau)/2
						(τ, τ)
					end
					set_branch_length!(a, τ1, 1)
					set_branch_length!(a, τ2, 2)
				end
			end
		end
	end
end

"""
	shared_nodes(t1::Tree, t2::Tree, MCCs)

Return two dictionaries containing the status of nodes for each tree.
*E.g.* `S::Dict` for `t1` such that for `n1` a node in `t1`, `S[n1] = (status, n2, isroot, id_mcc)` with:
- `status ∈ {:non_shared, :shared, :mcc_root, :shared_singleton}`
- `n2`: if `status==:shared` or `status==:mcc_root`, node in `t2`. Otherwise, `nothing`.
- `isroot` if `n2` is the root of `t2` (i.e. root of other tree)
- `id_mcc`: if `status!=:non_shared`, identifier for MCC
"""
function shared_nodes(t1::Tree, t2::Tree, MCCs)
	X1 = Dict()
	X2 = Dict()
	for (i,mcc) in enumerate(MCCs)
		shared_nodes!(X1, X2, t1, t2, mcc, i)
	end
	for n1 in internals(t1)
		if !haskey(X1, n1.label)
			X1[n1.label] = (:non_shared, nothing, false, nothing)
		end
	end
	for n2 in internals(t2)
		if !haskey(X2, n2.label)
			X2[n2.label] = (:non_shared, nothing, false, nothing)
		end
	end

	# Quick test of reflectivity
	for (k1, v1) in X1
		if v1[1] == :shared && (X2[v1[2]][1] != :shared || X2[v1[2]][2] != k1)
			error("Inconsistent shared nodes")
		end
	end

	return X1, X2
end
#=
Shared nodes `X` of tree `t`, using `t2` as a second tree, for a given mcc.
=#
function shared_nodes!(X1::Dict, X2::Dict, t1::Tree, t2::Tree, mcc, id_mcc)
	# To identify singletons, we'll need the splits restricted to the mcc
	Sm = TreeKnit.splits_in_mcc(mcc, t1, t2)

	# Deal with the root first
	A1 = lca(t1, mcc)
	A2 = lca(t2, mcc)
	X1[A1.label] = (:mcc_root, A2.label, A2.isroot, id_mcc)
	X2[A2.label] = (:mcc_root, A1.label, A1.isroot, id_mcc)

	# Now adding internal nodes
	for label in mcc
		# Going up the two trees at the same time starting from this leaf
		n1 = t1.lnodes[label]
		n2 = t2.lnodes[label]
		#  A leaf is always shared
		!haskey(X1, n1.label) && (X1[n1.label] = (:shared, n2.label, n2.isroot, id_mcc))
		!haskey(X2, n2.label) && (X2[n2.label] = (:shared, n1.label, n1.isroot, id_mcc))

		a1 = n1.anc
		a2 = n2.anc
		while (n1 != A1 || n2 != A2) && (!haskey(X1, a1.label) || !haskey(X2, a2.label))
			# Are a1/a2 shared / shared singleton?
			f1 = is_shared_singleton(a1, n1, Sm[1]) ? :shared_singleton : :shared
			f2 = is_shared_singleton(a2, n2, Sm[2]) ? :shared_singleton : :shared
			if f1 == :shared && f2 == :shared
				# Both are shared, they are equivalent. Go up in the trees
				!haskey(X1, a1.label) && (X1[a1.label] = (:shared, a2.label, a2.isroot, id_mcc))
				!haskey(X2, a2.label) && (X2[a2.label] = (:shared, a1.label, a1.isroot, id_mcc))
				n1 = a1
				n2 = a2
				a1 = n1.anc
				a2 = n2.anc
				continue
			end
			if f1 == :shared_singleton
				# a1 does not have an equivalent in t2, but is shared. Go up in t1
				X1[a1.label] = (:shared_singleton, nothing, false, id_mcc)
				n1 = a1
				a1 = n1.anc
			end
			if f2 == :shared_singleton
				# a2 does not have an equivalent in t1, but is shared. Go up in t2
				X2[a2.label] = (:shared_singleton, nothing, false, id_mcc)
				n2 = a2
				a2 = n2.anc
			end
		end
	end
end

"""
	is_shared_singleton(n::TreeNode, c::TreeNode, Sm::SplitList)

Check whether `n` is a shared node that would be a singleton in one of the trees.
This is the case if
	- the splits defined by `n` and `c` are the same when restricted to an mcc
	- the split defined by `n` is a leaf split when restricted to an mcc

It's necessary to check this second condition since leaves are not stored in `SplitList`
  objects.

## Input
  - node `n`
  - node `c`, child of `n`
  - `Sm::TreeTools.SplitList`: splits defined by the MCC `c` belongs to.
"""
function is_shared_singleton(n::TreeNode, c::TreeNode, Sm::SplitList)
	if c.anc != n
		error("$(c.label) is not a child of $(n.label)")
	end
	if c.isleaf
		return TreeTools.is_leaf_split(Sm.splitmap[n.label], Sm.mask)
	else
		return Base.isequal(Sm.splitmap[n.label], Sm.splitmap[c.label], Sm.mask)
	end
end

"""
	fix_shared_singletons!(t1, t2, X1, X2, MCCs)

Introduce shared singletons of one tree in the other, to make reconstructing the ARG easier.
"""
function fix_shared_singletons!(t1, t2, X1, X2, MCCs)
	for (i,mcc) in enumerate(MCCs)
		fix_shared_singletons!(t1, t2, X1, X2, mcc)
	end
end
"""
	fix_shared_singletons!(t, Xref, tref)

Introduce shared singletons of `tref` into `t`.
"""
function fix_shared_singletons!(t1, t2, X1, X2, mcc::Vector{String})
	A1 = lca(t1, mcc...)
	A2 = lca(t2, mcc...)
	id_mcc = X1[A1.label][4]
	@assert X1[A1.label][1] == :mcc_root "$(A1)"
	@assert X2[A2.label][1] == :mcc_root "$(A2)"
	@assert X1[A1.label][2] == A2.label
	@assert X2[A2.label][2] == A1.label

	v1 = Dict()
	v2 = Dict()
	# Go up from each leaf in the two trees
	for m in mcc
		n1 = t1.lnodes[m]
		n2 = t2.lnodes[m]
		while n1 != A1 && !haskey(v1, n1.label) && n2 != A2 && !haskey(v2, n2.label)
			# look up in the two trees
			a1 = n1.anc
			a2 = n2.anc
			f1 = (X1[a1.label][1] == :shared_singleton)
			f2 = (X2[a2.label][1] == :shared_singleton)
			if f1
				if f2
					# Both are shared.
					# Choose one and fix it
					if !ismissing(n1.tau) && !ismissing(n2.tau)
						if n1.tau < n2.tau
							# a1 in t2
							introduce_singleton!(n2, a2, a1, n1.tau, X2, X1)
						else
							# a2 in t1
							introduce_singleton!(n1, a1, a2, n2.tau, X1, X2)
						end
					else
						# a1 in t2 by default
						introduce_singleton!(n2, a2, a1, n1.tau, X2, X1)
					end

				else
					# Introduce a1 in t2
					introduce_singleton!(n2, a2, a1, n1.tau, X2, X1)
				end
			elseif f2
				# Introduce a2 in a1
				introduce_singleton!(n1, a1, a2, n2.tau, X1, X2)
			end
			# At this point, n1 and n2 should have the same ancestor
			v1[n1.label] = true
			v2[n2.label] = true
			n1 = n1.anc
			n2 = n2.anc
			@assert X1[n1.label][2] == n2.label "1 - $(n1.label) - $(X1[n1.label]) - $(n2.label)"
			@assert X2[n2.label][2] == n1.label "2 - $(n2.label) - $(X2[n2.label])"
		end
	end
	node2tree!(t1, t1.root)
	node2tree!(t2, t2.root)
end

function introduce_singleton!(n, a, sref, τ, X, Xref)
	_τ = ismissing(n.tau) ? missing : τ # Happens if one tree has times and the other not
	s = TreeTools.add_internal_singleton!(n, a, _τ, make_random_label("Singleton"))
	X[s.label] = (:shared, sref.label, false, X[n.label][4])
	Xref[sref.label] = (:shared, s.label, false, Xref[sref.label][4])

	return nothing
end
