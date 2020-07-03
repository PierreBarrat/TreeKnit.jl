"""
	trees_from_ARG(arg::ARG)
	trees_from_ARG!(ar::ARGNode, c::Int64)
	trees_from_ARG!(an::ARGNode, r::TreeNode, c::Int64)
"""
function trees_from_ARG(arg::ARG; prune_singletons = true)
	treelist = Array{Tree,1}(undef,0)
	for c in 1:arg.degree
		ar = arg.root[c]
		tr = trees_from_ARG!(ar, c)
		if prune_singletons
			push!(treelist, TreeTools.remove_internal_singletons!(node2tree(tr)))
		else
			push!(treelist, TreeTools.node2tree(tr))
		end
	end
	return treelist
end
function trees_from_ARG!(an::ARGNode, r::TreeNode, c::Int64)
	ic = findfirst(x->x[c], an.anccolor)
	n = TreeNode(an.data[ic],
		anc = r,
		isleaf = an.isleaf,
		isroot = false,
		label=an.label
		)
	for ac in ARGTools.get_children(an, c)
		tc = trees_from_ARG!(ac, n, c)
		push!(n.child, tc)
	end
	return n
end
function trees_from_ARG!(ar::ARGNode, c::Int64)
	if !ar.isroot[c]
		@error "Can't root tree on non-root `ARGNode` $(an.label)."
	end
	ic = findfirst(x->x[c], ar.anccolor)
		if !isnothing(ar.anc[ic])
			println("problem")
		end
	tr = TreeNode(ar.data[ic],
		anc = nothing,
		isleaf = ar.isleaf,
		isroot = true,
		label = ar.label
		)
	for ac in ARGTools.get_children(ar, c)
		tc = trees_from_ARG!(ac, tr, c)
		push!(tr.child, tc)
	end
	return tr
end

"""
	ARG_from_known_trees(trees::Vararg{Tree})

Construct an `ARG` from *known* trees. 
Labels of internal nodes in all input trees should correspond to an individual in the past population. Nodes with the same labels in different trees will be matched to the same `ARGNode` object. 

## Method 
Take `first(trees)` as a base `ARG`, and glue others to it. 
"""
function ARG_from_known_trees(trees::Vararg{Tree})
	reftree = trees[1]
	degree = length(trees)
	arg = ARG_from_tree(reftree, degree = degree)
	for i in 2:length(trees)
		glue_tree_to_ARG!(arg, trees[i], i)
	end
	return arg
end

function glue_tree_to_ARG!(arg::ARG, tree::Tree, color::Int64)
	glue_node_to_ARG!(arg, tree.root, nothing, color)
end

function glue_node_to_ARG!(arg::ARG, tn::TreeNode, ARGanc::Union{ARGNode, Nothing}, color::Int64)
	if haskey(arg.nodes, tn.label)
		# Update already existing ARG node
		an = arg.nodes[tn.label]
		an.color[color] = true
		# Special case of the root --> can't check label of ancestor in this case
		if isnothing(ARGanc)
			i = findfirst(x->x==nothing, an.anc)
			if !isnothing(i)
				an.anccolor[i][color] = true
			else
				push!(an.anc, ARGanc)
				push!(an.anccolor, _color(color, arg.degree))
			end
		else
			i = findfirst(x->x==ARGanc.label, [x.label for x in an.anc]) 
			if !isnothing(i)
				# if `ARGanc` is already the ancestor of `an`, set corresponding color to true
				an.anccolor[i][color] = true
			else
				# Else, add a new ancestor to `an`
				push!(an.anc, ARGanc)
				push!(an.anccolor, _color(color, arg.degree))
			end
		end
		push!(an.data, tn.data)
		push!(an.isroot, tn.isroot)
		an.isleaf = tn.isleaf
	else
		# Create a new ARG node
		an = ARGNode(degree=arg.degree,
			anc = [ARGanc],
			anccolor = [_color(color, arg.degree)],
			color = _color(color, arg.degree), 
			label = tn.label,
			data = [tn.data],
			isroot = _color(color, arg.degree)*tn.isroot,
			isleaf = tn.isleaf)
		arg.nodes[an.label] = an
	end
	# 
	for c in tn.child
		cn = glue_node_to_ARG!(arg, c, an, color)
		if !in(cn.label, [x.label for x in an.children])
			push!(an.children, cn)
		end
	end
	return an
end

"""
	ARG_from_tree(tree::Tree ; degree = 1)

Build an `ARG` from a single tree. If `degree > 1`, the output will not be a valid `ARG` object. 
"""
function ARG_from_tree(tree::Tree ; degree = 1)
	# Root
	arg = ARG(degree=degree)
	root = ARGNode(degree = degree,
		anc = [nothing], 
		anccolor = [_color(1, degree)],
		color = _color(1, degree),
		label = tree.root.label,
		data = [tree.root.data],
		isroot = _color(1, degree), 
		isleaf = tree.root.isleaf
		)
	arg.root[1] = root
	arg.nodes[root.label] = root
	# 
	for c in tree.root.child
		nc = ARGNode_from_TreeNode!(arg, c, root, degree=degree)
		push!(root.children, nc)
		arg.nodes[nc.label] = nc
	end
	return arg
end

"""
	ARGNode_from_TreeNode(tn::TreeNode, ARGanc::ARGNode; degree = 1)

Simple `ARGNode` from a `TreeNode` object, with ancestor `ARGanc` as first ancestor. If `degree > 1`, the first element of color is true and the others are false. 
"""
function ARGNode_from_TreeNode!(arg::ARG, tn::TreeNode, ARGanc::ARGNode; degree = 1)
	n = ARGNode(degree=degree,
		anc = [ARGanc],
		anccolor = [_color(1, degree)],
		color = _color(1, degree), 
		label = tn.label,
		data = [tn.data],
		isroot = _color(1, degree)*tn.isroot,
		isleaf = tn.isleaf)
	for c in tn.child
		nc = ARGNode_from_TreeNode!(arg, c, n, degree = degree)
		push!(n.children, nc)
		arg.nodes[nc.label] = nc
	end
	return n
end

# """
# """