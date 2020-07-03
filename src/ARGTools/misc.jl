# Check that ancestors of n have the color stored in n.anccolor
function check_anc(n::ARGNode)
	for (a,ac) in zip(n.anc, n.anccolor)
		if !isnothing(a)
			for (i,c) in enumerate(ac)
				if c && !a.color[i]
					@error "Node $(n.label) : inconsistent ancestor colors"
				end
			end
		end
	end
end
# Check that the sum of colors of children of n is the degree of n
function check_children(n::ARGNode)
	if !n.isleaf
		x = zeros(Bool, length(n.color))
		for c in n.children
			i = findfirst(x->x==n, c.anc)
			x .= (|).(x, c.anccolor[i])
		end
		if x != n.color
			error("Node $(n.label) : inconsistent children colors")
		end
	end
end
function check_color_consistency(n::ARGNode)
	# Each color must be found once in n.anccolor
	for (ic, c) in enumerate(n.color)
		cnt = 0
		for ac in n.anccolor
			cnt += Int64(ac[ic])
		end
		if cnt != Int64(c)
			error("Node $(n.label): inconsistent self colors.")
		end
	end
end
#
function check_arg(arg::ARG)
	for n in values(arg.nodes)
		check_arg_node(n)
	end
end
function check_arg_node(n::ARGNode)
	# Check that ancestors of n have the color stored in n.anccolor
	check_anc(n)
	# Check that the sum of colors of children of n is the degree of n
	check_children(n)
	# Check that colors in n.anccolor are consistent with n.color
	check_color_consistency(n)
end


### 
function has_consistent_colors(arg::ARG)
	for n in values(arg.nodes)
		if !consistent_children_colors(n) || !consistent_ancestors_colors(n)
			return false
		end
	end
	return true
end
function consistent_children_colors(n::ARGNode)
	if !n.isleaf
		x = zeros(Bool, length(n.color))
		for c in n.children
			i = findfirst(x->x==n, c.anc)
			x .= (|).(x, c.anccolor[i])
		end
		if x != n.color
			return false
		end
	end
	return true
end
function consistent_ancestors_colors(n::ARGNode)
	for (a,ac) in zip(n.anc, n.anccolor)
		if !isnothing(a)
			for (i,c) in enumerate(ac)
				if c && !a.color[i]
					return false
				end
			end
		end
	end	
	return true
end