"""
	arg_from_trees(t1, t2, MCCs)

Construct an `ARG` from two trees and a set of MCCs.
"""

function arg_from_trees(input_trees::Vector{Tree{T}}, MCCs::MCC_set) where T
    ##assume trees have been resolved compatibly with each other and MCCs are consistent
	if any([!ismissing(TreeTools.branch_length(t.root)) for t in input_trees])
		@warn "Root node of one of the input trees has a non missing branch length.
		This may cause issues when setting the branch length of the ARG"
	end
	#
	trees = [copy(t) for t in input_trees]
    X = Dict{Tuple{Int, Int}, Any}()
    for i in 1:(MCCs.no_trees-1)
        for j in (i+1):MCCs.no_trees
            # Find shared nodes in both trees
	        X[(i,j)] = SRG.shared_nodes(trees[i], trees[j], get(MCCs, (i,j)))
            # Introduce shared singletons of one into the other, add these to previous dict
			fix_multitree_singletons!(X, trees, i, j)
	        fix_shared_singletons!(trees, X, i,j, get(MCCs, (i,j)))
        end
    end
	# Make an initial ARG with the first tree
	arg, lm1 = arg_from_tree(trees[1], X, MCCs.no_trees, 1)
	lm_dict = Dict()
	lm_dict[1] = lm1
	# Add the second tree
    for k in 2:MCCs.no_trees
        lm = add_tree_to_arg!(arg, trees[k], X, lm_dict, MCCs.no_trees, k)
		lm_dict[k] = lm
        # Point their root to the same node
        # set_global_root!(arg)
        # Reverse label map
        ## rlm = reverse_label_map(arg, lm1, lm2)
        # Deal with branch lengths of the ARG
        ## set_branch_length!(arg, t1, t2, rlm)
    end

	return arg, lm_dict
end

function fix_multitree_singletons!(X_dict, trees, i, j)
	if i>1
		for (key, node) in X_dict[(i,j)][1]
			if node[1] == :shared_singleton
				for k in 1:(i-1)
					shared_node_tree_k = X_dict[(k,i)][2][key][2]
					if haskey(X_dict[(k,j)][1], shared_node_tree_k)
						shared_node_tree_j = X_dict[(k,j)][1][shared_node_tree_k][2]
						if X_dict[(i,j)][2][shared_node_tree_j][1] == :shared_singleton
							#mark as shared
							mcc_i = X_dict[(i,j)][1][trees[i].lnodes[key].anc.label][4]
							mcc_j = X_dict[(i,j)][2][trees[j].lnodes[shared_node_tree_j].anc.label][4]
							X_dict[(i,j)][1][key] = (:shared, shared_node_tree_j, false, mcc_i)
							X_dict[(i,j)][2][shared_node_tree_j] = (:shared, key, false, mcc_j)
						end
					end
				end
			end
		end
	end
end

"""
	fix_shared_singletons!(t1, t2, X1, X2, MCCs)

Introduce shared singletons of one tree in the other, to make reconstructing the ARG easier.
"""
function fix_shared_singletons!(trees, X_dict, i, j, MCC)
	for (k,mcc) in enumerate(MCC)
		fix_shared_singletons!(trees, X_dict, i, j, mcc)
	end
end
"""
	fix_shared_singletons!(t, Xref, tref)

Introduce shared singletons of `tref` into `t`.
"""
function fix_shared_singletons!(trees, X_dict, i, j, mcc::Vector{String})
    t1 = trees[i]
    t2 = trees[j]
    X1, X2 = X_dict[(i,j)]
    x1_toadd = [X_dict[(i, k)] for k in 2:(j-1) if i != k]
    trees1_toadd = [trees[k] for k in 2:(j-1) if i != k]
    x2_toadd = [X_dict[(k, j)] for k in 1:(i-1)]
    trees2_toadd = [trees[k] for k in 1:(i-1)]
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
			# println("$(n1.label) - $(n1.tau) - $(a1.label)")
			# println("$(n2.label) - $(n2.tau) - $(a2.label)")
			if f1
				if f2
					# Both are shared.
					# Choose one and fix it
					if !ismissing(n1.tau) && !ismissing(n2.tau)
						if n1.tau < n2.tau
							# a1 in t2
							introduce_singleton!(n2, a2, a1, n1.tau, X2, X1, x2_toadd, 2, trees2_toadd)
						else
							# a2 in t1
							introduce_singleton!(n1, a1, a2, n2.tau, X1, X2, x1_toadd, 1, trees1_toadd)
						end
					else
						# a1 in t2 by default
						introduce_singleton!(n2, a2, a1, n1.tau, X2, X1, x2_toadd, 2, trees2_toadd)
					end
				else
					# Introduce a1 in t2
					introduce_singleton!(n2, a2, a1, n1.tau, X2, X1, x2_toadd, 2, trees2_toadd)
				end
			elseif f2
				# Introduce a2 in a1
				introduce_singleton!(n1, a1, a2, n2.tau, X1, X2, x1_toadd, 1, trees1_toadd)
			end
			# At this point, n1 and n2 should have the same ancestor
			v1[n1.label] = true
			v2[n2.label] = true
			n1 = n1.anc
			n2 = n2.anc
			#@assert X1[n1.label][2] == n2.label "1 - $(n1.label) - $(X1[n1.label]) - $(n2.label)"
			#@assert X2[n2.label][2] == n1.label "2 - $(n2.label) - $(X2[n2.label])"
		end
	end
	node2tree!(t1, t1.root)
	node2tree!(t2, t2.root)
end

function introduce_singleton!(n, a, sref, τ, X, Xref, X_toadd, pos, trees)
	# println("Trying $(n.label) < $(sref.label) < $(a.label)")
	s = TreeTools.add_internal_singleton!(n, a, τ, SRG.make_random_label("Singleton"))
    for (x,t) in zip(X_toadd, trees)
        ##add this to the dictionary at position pos
        n_x = x[pos][n.label]
        a_x = x[pos][a.label]
        if (n_x[1]== :shared || n_x[1]== :shared_singleton) && (a_x[1]== :shared || a_x[1]== :shared_singleton || a_x[1]==:mcc_root)
            status = :shared 
            id = n_x[4]
        else
            status = :not_shared
			id = nothing
        end
        ##now add this to neighboring tree and add to dictionary at other pos
        nref_name = n_x[2]
        aref_name = a_x[2]
        sref_x = TreeTools.add_internal_singleton!(t.lnodes[nref_name], t.lnodes[aref_name], τ, SRG.make_random_label("Singleton"))
        x[pos][s.label] = (status, sref_x.label, false, id)
	    if pos ==1
            other_pos = 2
        else
            other_pos =1
        end
        x[other_pos][sref_x.label] = (status, s.label, false, id)
    end
	trees = [node2tree!(t, t.root) for t in trees]
	X[s.label] = (:shared, sref.label, false, X[n.label][4])
	Xref[sref.label] = (:shared, s.label, false, Xref[sref.label][4])

	return nothing
end

"""
	arg_from_tree(t::Tree, SN::Dict, clr=1)

Build an initial `ARG` using one tree.
"""
function arg_from_tree(t::Tree, X_dict::Dict, no_trees::Int, clr=1)
	label_map = Dict{String, Any}()
	arg = ARG(Dict(), Dict(), Dict(), Dict())

	# Intialize root
	a_R = ARGNode(;no_trees=no_trees)
	add_color!(a_R, clr)
	set_root!(a_R, clr)
	t.root.isleaf && set_label!(a_R, t.root.label)
	label_map[t.root.label] = a_R.label
	arg.nodes[a_R.label] = a_R
	arg.roots[clr] = a_R

	# Go down in tree
	for t_c in t.root.child
		grow_arg_from_treenode!(arg, a_R, t_c, X_dict, label_map, clr, no_trees)
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
	X_dict::Dict,
	label_map::Dict,
	clr,
	no_trees,
)
	treenodetype = [X_dict[(clr,k)][1][t_n.label][1] for k in 1:no_trees if clr<k] # Is t_n shared / mcc_root / non shared
	append!(treenodetype, [X_dict[(k,clr)][1][t_n.label][1] for k in 1:no_trees if k<clr] )
	root_in_other_tree = [X_dict[(clr,k)][1][t_n.label][3] for k in 1:no_trees if clr<k] # Is t_n the root in the other tree
	append!(root_in_other_tree, [X_dict[(k,clr)][1][t_n.label][3] for k in 1:no_trees if k<clr])

	# Creating an ARGNode to represent t_n
	## If t_n is the root of an MCC (and not the root of the other tree)
	## Create a hybrid node a_h, a node a_n for t_n, and do a_a > a_h > a_n
	## Otherwise, create a node a_n for t_n, and do a_a > a_n
	a_n, a_h = if any([t == :mcc_root && !r for (t,r) in zip(treenodetype, root_in_other_tree)])
		a_h = ARGNode(;hybrid=true, no_trees=no_trees)
		a_n = ARGNode(;no_trees=no_trees)
		add_color!(a_h, clr)
		add_color!(a_n, clr)
		graft!(a_a, a_h, clr)
		graft!(a_h, a_n, clr)

		(a_n, a_h)
	else
		a_n = ARGNode(;no_trees=no_trees)
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
			arg, a_n, t_c, X_dict, label_map, clr, no_trees
		)
	end

	return nothing
end

"""
	add_tree_to_arg!(arg, t, SN, lmref, clr = 2)

Add tree to an existing `ARG`, using color `clr`.
This is only designed for an ARG of two trees.
"""
function add_tree_to_arg!(arg, t, X_dict, lmref_dict, no_trees, clr = 2)
	lm = Dict()

	# Deal with root of t
	a_R = if all([X_dict[(k, clr)][2][t.root.label][1] == :non_shared for k in 1:(clr-1)])
		ARGNode(;no_trees=no_trees)
	else
		node_name = [(X_dict[(k, clr)][2][t.root.label][2], k) for k in 1:(clr-1)][1]
		n2 = node_name[1] # Label of t.root in ref tree
		pos = node_name[2]
		arg.nodes[lmref_dict[pos][n2]]
	end
	add_color!(a_R, clr)
	set_root!(a_R, clr)

	arg.roots[clr] = a_R
	arg.nodes[label(a_R)] = a_R
	lm[t.root.label] = a_R.label

	# Go down t
	for t_c in t.root.child
		add_treenode_to_arg!(arg, a_R, t_c, X_dict, lm, lmref_dict, no_trees, clr)
	end

	return lm
end


function add_treenode_to_arg!(
	arg::ARG,
	a_a::AbstractARGNode,
	t_n::TreeTools.TreeNode,
	X_dict::Dict,
	lm::Dict,
	lmref_dict::Dict,
	no_trees::Int, 
	clr
)
	treenodetype = [X_dict[(k, clr)][2][t_n.label][1] for k in 1:(clr-1)]
	root_in_other_tree = [X_dict[(k, clr)][2][t_n.label][3] for k in 1:(clr-1)]

	a_n = if all([t == :non_shared for t in treenodetype])
		# New node : we create a new argnode a_n
		a_n = ARGNode(;no_trees=no_trees)
		add_color!(a_n, clr)
		graft!(a_a, a_n, clr)
		arg.nodes[label(a_n)] = a_n
		a_n
	elseif any([t == :mcc_root for t in treenodetype])
		pos = findall( x -> x == :mcc_root, treenodetype)[1]
		# t_n is either (i) the root of the other tree, or (ii) the point of a reassortment
		# In either case, there is already a node a_n representing it
		n2 = X_dict[(pos, clr)][2][t_n.label][2]
		a_n = arg.nodes[lmref_dict[pos][n2]]
		add_color!(a_n, clr)
		if any([ r== true for r in root_in_other_tree])
			# (i)
			graft!(a_a, a_n, clr)
		else
			# (ii): Find the corresponding hybrid node
			pos = findall( x -> x == false, root_in_other_tree)[1]
			a_h = ancestor(a_n, pos)
			@assert ishybrid(a_h) "$label(a_h)"
			add_color!(a_h, clr)
			graft!(a_a, a_h, clr)
			graft!(a_h, a_n, clr)
		end
		a_n
	else
		# t_n is shared and not the root of an MCC
		# --> there is already an ARGNode corresponding to it
		n2 = X_dict[(1, clr)][2][t_n.label][2]
		a_n = arg.nodes[lmref_dict[1][n2]]
		add_color!(a_n, clr)
		@assert ancestor(a_n, 1) == a_a
		graft!(a_a, a_n, clr)
		a_n
	end

	# Adding to label map
	lm[t_n.label] = label(a_n)

	# Recursing down the tree
	for t_c in t_n.child
		add_treenode_to_arg!(
			arg, a_n, t_c, X_dict, lm, lmref_dict, no_trees, clr
		)
	end

	return nothing
end