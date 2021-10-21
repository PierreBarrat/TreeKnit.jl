function trees_from_arg(arg::ARG)
	trees = Vector{Any}(undef, 2)
	for clr in 1:2
		trees[clr] = tree_from_arg(arg, clr)
	end
	return (trees[1], trees[2])
end

function tree_from_arg(arg::ARG, clr::Int)
	r = tree_from_arg!(nothing, arg.roots[clr], clr)
	return TreeTools.node2tree(r)
end

function tree_from_arg!(t_a, a_n::ARGNode, clr)
	@assert hascolor(a_n, clr)
	# Create new tree node
	n = TreeTools.TreeNode(; label = label(a_n), tau = branch_length(a_n, clr))
	if !isroot(a_n, clr)
		@assert !isnothing(t_a)
		TreeTools.graftnode!(t_a, n)
	end
	# Recurse
	for a_c in children(a_n, clr)
		tree_from_arg!(n, a_c, clr)
	end

	return n
end
