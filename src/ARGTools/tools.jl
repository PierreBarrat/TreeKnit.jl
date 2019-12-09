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

"""
Count the number of nodes in `arg` that have more than one ancestor.
"""
function count_reassortments(arg::ARG)
	rea = [k for (k,v) in arg.nodes if length(v.anc) > 1]
	return length(rea), rea
end

"""
	regraft!(n::ARGNode, oldanc::ARGNode, newanc::ARGNode, color::Int64)

Regraft node `n` from `oldanc` to `newanc` for color `color`. Base routine. 
## Note
- `oldanc` has to be an ancestor of `n` for color `color`
- `newanc` and `n` have to be of color `color`
"""
function regraft!(n::ARGNode, oldanc::ARGNode, newanc::ARGNode, color::Int64)
	i_old = findfirst(x->x.label==oldanc.label, n.anc)
	i_child = findfirst(x->x.label==n.label, oldanc.children)
	if isnothing(i_old)
		error("Attempting to regraft node from uncorrect ancestor: $(oldanc.label) not ancestor of $(n.label)")
	elseif isnothing(i_child)
		error("Attempting to regraft node from uncorrect ancestor: $(n.label) not child of $(oldanc.label)")
	elseif length(oldanc.children) == 1
		error("Removing unique child from internal node ($(oldanc.label), $(n.label))")
	end
	if !n.anccolor[i_old][color]
		error("Attempting to regraft node from uncorrect ancestor: branch from $(n.label) to $(oldanc.label) is not of color $color.")
	end
	if !newanc.color[color]
		error("Attempting to regraft node to ancestor of incorrect color.")
	elseif !n.color[color]
		error("Attempting to regraft node for the wrong color.")
	end

	regraft!(n, oldanc, newanc, color, i_old)
end
"""
	regraft!(n::ARGNode, oldanc::ARGNode, newanc::ARGNode, color::Int64, i_old::Int64)

Core function for regrafting. Does not handle errors.
"""
function regraft!(n::ARGNode, oldanc::ARGNode, newanc::ARGNode, color::Int64, i_old::Int64)
	# Changes for n
	i_new = findfirst(x->x.label==newanc.label, n.anc) # Was `newanc` already an ancestor of `n` for another color? 
	spliceflag = sum(n.anccolor[i_old]) == 1 # Was `oldanc` the ancestor for `color` only? If yes, should be removed. If not, should be kept
	if isnothing(i_new) # We have to create a new ancestor for `n`
		push!(n.anc, newanc)
		push!(n.anccolor, _color(color))
	else # Just set the branch of the right color
		n.anccolor[i_new][color] = true
	end
	if !spliceflag 
		n.anccolor[i_old][color] = false
	else
		splice!(n.anccolor, i_old)
		splice!(n.anc, i_old)
	end

	# Changes for oldanc
	i_child = findfirst(x->x.label==n.label, oldanc.children)
	splice!(oldanc.children, i_child)

	# Changes for newanc
	i_child = findfirst(x->x.label==n.label, newanc.children)
	if isnothing(i_child) 
		push!(newanc.children, n)
	end
end

"""
Find the below situation
a1***
|	*
|	a2
|	*
n****
where "|" and "*" are two colors. Fix it. (See source code of function help for correct display)

## Method
Find triplets `(n, a1, a2)` such that `a1, a2` are ancestors of `n`, and `a1` is an ancestor of `a2` for `colors`. 
Then, the triplet is validated if the following conditions are met for two colors `clr1` and `clr2` of `n`. 
- `a1` is an ancestor of `n` for `clr1`. 
- `a2` is an ancestor of `n` for `clr2`. 
- `a1` is an ancestor of `a2` for `clr2`. 
- `a2` is *not* of color `clr1`. 
"""
function find_trivial_loop!(arg::ARG)
	nloops = 0
	for n in values(arg.nodes)
		# Finding loops that need fixing, then fixing. Fixing loops modifies indices of ancestors, so I have to proceed this way. 
		loops_to_fix = Array{Tuple{String, String, Int64, Int64},1}(undef, 0)
		for (i1, a1) in enumerate(n.anc)
			for (i2, a2) in enumerate(n.anc)
				flag, colors = is_ancestor(a1, a2)
				if flag # `a1` is an ancestor of `a2` for `colors`
					for clr1 in findall(n.color)
						for clr2 in findall(n.color)
							if n.anccolor[i1][clr1] && n.anccolor[i2][clr2] && colors[clr2] && !a2.color[clr1]
								nloops += 1
								push!(loops_to_fix, (a1.label, a2.label, clr1, clr2))
							end
						end
					end
				end
			end
		end
		for (l1, l2, clr1, clr2) in loops_to_fix
			println("Fixing loop ($l1, $(n.label), $l2, $clr1, $clr2) ")
			fix_trivial_loop!(arg.nodes[l1], n, arg.nodes[l2], clr1, clr2)
		end
	end
	println("`find_trivial_loop!`: Found $nloops loops to fix.")
end

"""
	fix_trivial_loop!(A::ARGNode, B::ARGNode, C::ARGNode, clr_1::Int64, clr_2::Int64)

Place `C` on the branch going up from `B` to `A` for `color`. Obviously, `C` should not already be of color `color`.  
A***			A     
|  *			|*  
|  C   -->  	C  
|  *			|*  
B***			B  
With "|" corresponding to `clr_1` and "*" to `clr_2`. (See source code of function help for correct display)
## Steps
- Set `C.color[clr_1]` to `true`
- Set branch from `C` to `A` to `clr_1` on top of `clr_2`. 
- Regraft `B` onto `C` for color `clr_1`. 
## Checks to perform
- `C` must be of color `clr_2` and not of `clr_1`. 
- `A` and `B` must be of color `clr_1` and `clr_2`.
- `B` must be a child of `A` and `C` for resp. `clr_1` and `clr_2`
- `C` must be a child of `A` for `clr_2`. 
"""
function fix_trivial_loop!(A::ARGNode, B::ARGNode, C::ARGNode, clr_1::Int64, clr_2::Int64)
	# Checks
	flag = true
	if !((!C.color[clr_1]) && C.color[clr_2])
		println("`C` is not of the correct color")
		flag = false
	end
	flag *= A.color[clr_1] && A.color[clr_2] && B.color[clr_1] && B.color[clr_2]
	flag *= is_ancestor(A, B, clr_1)
	flag *= is_ancestor(C, B, clr_2)
	flag *= is_ancestor(A, C, clr_2)
	if !flag
		error("Conditions for fixing trivial loop are not met.")
	end

	# Fixing
	C.color[clr_1] = true
	i = findfirst(x->x.label==A.label, C.anc)
	C.anccolor[i][clr_1] = true
	regraft!(B, A, C, clr_1)
end

"""

Check if `a` is an ancestor of `c` for colors in `color`. Also check whether `c` is in `a.children`. 		
"""
function is_ancestor(a::ARGNode, c::ARGNode, color::Vararg{Int64})
	i_anc = findfirst(x->!isnothing(x) && x.label==a.label, c.anc)
	flag = !isnothing(i_anc) && !isnothing(findfirst(x->x.label==c.label, a.children)) 
	if flag
		for clr in color
			flag *= c.anccolor[i_anc][clr]
		end
	end
	return flag
end

"""

Check if `a` is an ancestor of `c`. Return a `Bool` as well as an array of colors for which `a` is an ancestor of `c`. 		
"""
function is_ancestor(a::ARGNode, c::ARGNode)
	out = zeros(Bool, c.degree)
	flag = false
	for clr in findall(c.color)
		if is_ancestor(a,c,clr)
			out[clr] = true
			flag = true
		end
	end
	return flag, out
end

"""
"""
function is_ancestor(a::Nothing, c::ARGNode)
	flag = false
	out = zeros(Bool, c.degree)
	for i in findall(c.isroot)
		flag = true
		out[findall(c.anccolor[i])] .= true
	end
	return flag, out
end
function is_ancestor(a::Nothing, c::Nothing)
	return false, zeros(Bool, 0)
end
function is_ancestor(a::ARGNode, c::Nothing)
	return false, zeros(Bool, 0)
end