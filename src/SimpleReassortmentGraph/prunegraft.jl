function graft!(a, n, clr::Int, colors::Vararg{Int})
	graft!(a, n, clr)
	for c in colors
		graft!(a, n, c)
	end

	return nothing
end
function graft!(a, n, colors::Vector{Bool})
	for (i,c) in enumerate(colors)
		if c
			graft!(a, n, i)
		end
	end

	return nothing
end
function graft!(a::ARGNode, n::ARGNode, c::Int)
	# Check that `n` does not have an ancestor
	@assert !hasancestor(n, c)
	# Proceed: set `a` as ancestor of `n`
	set_ancestor!(n, a, c)
	# If `n` is not already in children of `a`, put it there
	if !haschild(a, n)
		add_child!(a, n)
	end

	return nothing
end

function graft!(a::RootNode, n::ARGNode, c::Int)
	@assert !hasancestor(n, c)
	set_ancestor!(n, a, c)

	return nothing
end

function set_root!(n::ARGNode, c::Int)
	graft!(RootNode(), n, c)
end

function unset_child!(a::ARGNode, n, clr)
	@assert ancestor(n, clr) == a
	@assert haschild(a, n)
	# Is a the ancestor of n for another color?
	# If yes, we're not doing anything
	for (i,c) in color(n)
		if c && (i != clr)
			return nothing
		end
	end
	# Else, remove n from the children of a
	deleteat!(findfirst(==(n), a.children))
end
unset_child!(a, n, clr) = nothing

function prune!(n::ARGNode, c::Int)
	a = ancestor(n, c)
	unset_child!(a, n, c)
	set_ancestor!(n, nothing, c)
end

# function set_global_root!(arg::ARG)
# 	gr = RootNode([])
# 	for (clr, r) in arg.roots
# 		@assert isroot(r, clr)
# 		prune!(r, clr)
# 		graft!(gr, r, clr)
# 	end
# end
