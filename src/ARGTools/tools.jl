function regraft!(n::ARGNode, oldanc::ARGNode, newanc::ARGNode, color::Array{Bool,1})
	for c in findall(color)
		regraft!(n, oldanc, newanc, c)
	end
end
"""
	regraft!(n::ARGNode, oldanc::ARGNode, newanc::ARGNode, color::Array{Bool,1})
	regraft!(n::ARGNode, oldanc::ARGNode, newanc::ARGNode, color::Int64)

Regraft node `n` from `oldanc` to `newanc` for color `color`. Base routine. 
## Note
- `oldanc` has to be an ancestor of `n` for color `color`
- `newanc` and `n` have to be of color `color`
## Warning  
Does not handle `Nothing`. Regrafting a root node **or** to `Nothing` will fail. 
"""
function regraft!(n::ARGNode, oldanc::ARGNode, newanc::ARGNode, color::Int64; w=false)
	i_old = findfirst(x->x==oldanc, n.anc)
	i_child = findfirst(x->x==n, oldanc.children)
	# n.label=="internal_23" && println(n.anc)
	if isnothing(i_old)
		error("Attempting to regraft node from uncorrect ancestor: $(oldanc.label) not ancestor of $(n.label)")
	elseif isnothing(i_child)
		error("Attempting to regraft node from uncorrect ancestor: $(n.label) not child of $(oldanc.label)")
	elseif length(oldanc.children) == 1
		w && @warn "Removing unique child from internal node ($(oldanc.label), $(n.label))"
	end
	if !n.anccolor[i_old][color]
		error("Attempting to regraft node from uncorrect ancestor: branch from $(n.label) to $(oldanc.label) is not of color $color.")
	end
	if !newanc.color[color]
		error("Attempting to regraft node to ancestor of incorrect color.")
	elseif !n.color[color]
		error("Attempting to regraft node for the wrong color.")
	end
	#
	regraft!(n, oldanc, newanc, color, i_old)
end
"""
	regraft!(n::ARGNode, oldanc::ARGNode, newanc::ARGNode, color::Int64, i_old::Int64)

Core function for regrafting. Does not handle errors.
"""
function regraft!(n::ARGNode, oldanc::ARGNode, newanc::ARGNode, color::Int64, i_old::Int64)
	# n.label == "internal_32" && println(n)
	# Changes for n and newanc
	graft!(n, newanc, color) 
	# Changes for n and oldanc
	cut_branch!(n, i_old, color)
end
"""
	graft!(n::ARGNode, a::ARGNode, clr::Int64)

Graft `n` onto `a` for color `clr`. 
"""
function graft!(n::ARGNode, a::ARGNode, clr::Int64)
	# Changes for `a`
	if !isnothing(a) && !is_ancestor(a, n)[1]
		push!(a.children, n)
	end
	# Changes for `n`
	ia = findfirst(x->x==a, n.anc) # Was `a` already an ancestor of `n` for another color? 
	if isnothing(ia) # We have to create a new ancestor for `n`
		push!(n.anc, a)
		push!(n.anccolor, _color(clr, length(n.color)))
	else # Just set the branch of the right color
		n.anccolor[ia][clr] = true
	end
end

"""
	prune!(n::ARGNode; nochildren=true)

Prune `n`. Return `n`. If `nochildren`, `n` should not have children when pruned. 
"""
function prune!(arg::ARG, n::ARGNode)
	if length(n.children)>0
		@error "Trying to prune node `n.label` with children."
	end
	while !isempty(n.anc)
		for c in findall(n.anccolor[1])
			cut_branch!(n, 1, c)
		end
	end
	delete!(arg.nodes, n.label)
	return n
end

"""
Cut branch from `n` to `n.anc[i]` for color `clr`. Color of `n.anc[i]` is changed if necessary. 1
# Note
`n` is not regrafted and may not have any ancestor left after this.   
`n.anc` may have no child after this. 
"""
function cut_branch!(n::ARGNode, i::Int64, clr::Int64)
	# Checks
	length(n.anc) < i && @error "Cannot delete anc #$i from $(n.label): too few ancestors"
	# Delete branches
	a = n.anc[i]
	n.anccolor[i][clr] = false
	if !(|)(n.anccolor[i]...)
		deleteat!(n.anc, i)
		deleteat!(n.anccolor, i)
		if !isnothing(a)
			ic = findfirst(x->x.label==n.label, a.children)
			deleteat!(a.children, ic)
		end
	end
	# Correct colors if necessary
	correct_color!(a)
end
function cut_branch!(n::ARGNode, idx::Array{Int64}, clr::Int64)
	length(n.anc) < max(idx...) && @error "Cannot delete anc #$(max(idx...)) from $(n.label): too few ancestors"
	# 
	todel = Int64[]
	alist = []
	for i in idx
		a = n.anc[i]
		push!(alist, a)
		n.anccolor[i][clr] = false
		if !(|)(n.anccolor[i]...)
			push!(todel, i)
			if !isnothing(a)
				ic = findfirst(x->x.label==n.label, a.children)
				deleteat!(a.children, ic)
			end
		end
	end
	deleteat!(n.anc, todel)
	deleteat!(n.anccolor, todel)
	for a in alist
		correct_color!(a)
	end
end

###
###


function correct_color!(n::ARGNode)
	if !n.isleaf
		x = zeros(Bool, length(n.color))
		for c in n.children
			i = findfirst(x->x==n, c.anc)
			x .= (|).(x, c.anccolor[i])
		end
		n.color = x
	end
end
correct_color!(n::Nothing) = nothing

"""
	has_singletons(arg)

Singletons are nodes `n` such that `length(n.children) == 1`.
"""
function has_singletons(arg)
	for n in values(arg.nodes)
		if length(n.children) == 1
			return true
		end
	end
	return false
end
"""
	prune_singletons!(arg::ARG; v=false)

Singletons are nodes `n` such that `length(n.children) == 1`.  
"""
function prune_singletons!(arg::ARG; v=false, Nit=1e3)
	npruned = 1
	nit = 0
	while has_singletons(arg) && nit < Nit
		for n in values(arg.nodes)
			if length(n.children) == 1
				# Graft the child onto ancestors
				v && println("Pruned $(n.label).\n Ancestors $([x.label for x in n.anc]).\n Child $(n.children[1].label).")
				for (i,a) in enumerate(n.anc)
					clr = n.anccolor[i]
					regraft!(n.children[1], n, a, clr)
				end
				prune!(arg, n)
				npruned += 1
			end
		end
		if length(arg.nodes) == 1
			@warn "ARG with one node only."
			break
		end
		nit += 1
	end
	# check_arg(arg)
	v && println("Pruned $npruned nodes")
	return nothing
end

"""
	has_lone_nodes(arg::ARG)

A lone node `n` is such that `!n.isleaf && length(n.children)==0`. This can arise when simulating an ARG.
"""
function has_lone_nodes(arg::ARG)
	for n in values(arg.nodes)
		if !n.isleaf && length(n.children) == 0
			return true
		end
	end
	return false
end
"""
	prune_lone_nodes!(arg::ARG; v=false)

A lone node `n` is such that `!n.isleaf && length(n.children)==0`. This can arise when simulating an ARG.
"""
function prune_lone_nodes!(arg::ARG; v=false, Nit=1e3)
	npruned = 1
	nit = 0
	while has_lone_nodes(arg) && nit < Nit
		for n in values(arg.nodes)
			if !n.isleaf && length(n.children) == 0
				prune!(arg, n)
				npruned += 1
			end
		end
		if length(arg.nodes) == 1
			@warn "ARG with one node only."
			break
		end
		nit += 1
	end
	v && println("Pruned $npruned nodes")
	return nothing
end


"""
	is_ancestor(a::ARGNode, c::ARGNode, color::Vararg{Int64})

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
	is_ancestor(a::ARGNode, c::ARGNode)

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
is_ancestor(a::Nothing, c::ARGNode) = ((|)(c.isroot...), c.isroot)
is_ancestor(a::Nothing, c::Nothing) = (false, zeros(Bool, 0))
is_ancestor(a::ARGNode, c::Nothing) = (false, zeros(Bool, 0))

"""
	get_children_index(a::ARGNode, clr::Int64)

Get index of children of `a` for color `clr`. 
"""
function get_children_index(a::ARGNode, clr::Int64)
	idx = Int64[]
	if !a.color[clr]
		error("$(a.label) is not of color $clr")
	end
	for (ic,c) in enumerate(a.children)
		i = findfirst(x->x==a, c.anc)
		if c.anccolor[i][clr]
			push!(idx, ic)
		end
	end
	return idx
end
"""
	get_children(a::ARGNode, clr::Int64)

Get children of `a` for color `clr`. 
"""
get_children(a::ARGNode, clr::Int64) = a.children[get_children_index(a, clr)]

#########
# """
# Find the below situation
# a1***
# |	*
# |	a2
# |	*
# n****
# where "|" and "*" are two colors. This is a trivial and useless split and coalescence. Fix it. (See source code of function help for correct display of above situation)

# ## Method
# Find triplets `(n, a1, a2)` such that `a1, a2` are ancestors of `n`, and `a1` is an ancestor of `a2` for `colors`. 
# Then, the triplet is validated if the following conditions are met for two colors `clr1` and `clr2` of `n`. 
# - `a1` is an ancestor of `n` for `clr1`. 
# - `a2` is an ancestor of `n` for `clr2`. 
# - `a1` is an ancestor of `a2` for `clr2`. 
# - `a2` is *not* of color `clr1`. 
# """
# function find_trivial_loop!(arg::ARG)
# 	nloops = 0
# 	for n in values(arg.nodes)
# 		# Find loops that need fixing, then fix. Fixing loops modifies indices of ancestors, so I have to proceed this way. 
# 		loops_to_fix = Array{Tuple{String, String, Int64, Int64},1}(undef, 0)
# 		for (i1, a1) in enumerate(n.anc)
# 			for (i2, a2) in enumerate(n.anc)
# 				flag, colors = is_ancestor(a1, a2)
# 				if flag && !isnothing(a1) # `a1` is an ancestor of `a2` for `colors`
# 					for clr1 in findall(n.color)
# 						for clr2 in findall(n.color)
# 							if n.anccolor[i1][clr1] && n.anccolor[i2][clr2] && colors[clr2] && !a2.color[clr1]
# 								nloops += 1
# 								push!(loops_to_fix, (a1.label, a2.label, clr1, clr2))
# 							end
# 						end
# 					end
# 				end
# 			end
# 		end
# 		for (l1, l2, clr1, clr2) in loops_to_fix
# 			println("Fixing loop ($l1, $(n.label), $l2, $clr1, $clr2) ")
# 			fix_trivial_loop!(arg.nodes[l1], n, arg.nodes[l2], clr1, clr2)
# 		end
# 	end
# 	println("`find_trivial_loop!`: Found $nloops loops to fix.")
# end

# """
# 	fix_trivial_loop!(A::ARGNode, B::ARGNode, C::ARGNode, clr_1::Int64, clr_2::Int64)

# Place `C` on the branch going up from `B` to `A` for `color`. Obviously, `C` should not already be of color `color`.  
# A***			A     
# |  *			|*  
# |  C   -->  	C  
# |  *			|*  
# B***			B  
# With "|" corresponding to `clr_1` and "*" to `clr_2`. (See source code of function help for correct display)
# ## Steps
# - Set `C.color[clr_1]` to `true`
# - Set branch from `C` to `A` to `clr_1` on top of `clr_2`. 
# - Regraft `B` onto `C` for color `clr_1`. 
# ## Checks to perform
# - `C` must be of color `clr_2` and not of `clr_1`. 
# - `A` and `B` must be of color `clr_1` and `clr_2`.
# - `B` must be a child of `A` and `C` for resp. `clr_1` and `clr_2`
# - `C` must be a child of `A` for `clr_2`. 
# """
# function fix_trivial_loop!(A::ARGNode, B::ARGNode, C::ARGNode, clr_1::Int64, clr_2::Int64)
# 	# Checks
# 	flag = true
# 	if !((!C.color[clr_1]) && C.color[clr_2])
# 		println("`C` is not of the correct color")
# 		flag = false
# 	end
# 	flag *= A.color[clr_1] && A.color[clr_2] && B.color[clr_1] && B.color[clr_2]
# 	flag *= is_ancestor(A, B, clr_1)
# 	flag *= is_ancestor(C, B, clr_2)
# 	flag *= is_ancestor(A, C, clr_2)
# 	if !flag
# 		error("Conditions for fixing trivial loop are not met.")
# 	end

# 	# Fixing
# 	C.color[clr_1] = true
# 	i = findfirst(x->x.label==A.label, C.anc)
# 	C.anccolor[i][clr_1] = true
# 	regraft!(B, A, C, clr_1)
# end