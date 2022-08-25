abstract type AbstractARGNode end

"""
	mutable struct ARGNode

Node of an Ancestral Reassortment Graph representing two genes / segments.

## Fields:
```
	no_trees :: Int = 2
	anc :: Vector{Union{Nothing, AbstractARGNode}} = fill(nothing, no_trees)# Always of length 2
	color :: Vector{Bool} = fill(false, no_trees)
	tau :: Vector{Union{Missing, Float64}} = [missing, missing]
	children :: Vector{ARGNode} = ARGNode[]
	label :: String = make_random_label()
	hybrid :: Bool = false # Hybrids can have the same ancestor for some reassortments
```
"""
@with_kw_noshow mutable struct ARGNode <: AbstractARGNode
	no_trees :: Int = 2
	anc :: Vector{Union{Nothing, AbstractARGNode}} = fill(nothing, no_trees)# Always of length 2
	color :: Vector{Bool} = fill(false, no_trees)
	tau :: Vector{Union{Missing, Float64}} = [missing, missing]
	children :: Vector{ARGNode} = ARGNode[]
	label :: String = make_random_label()
	hybrid :: Bool = false # Hybrids can have the same ancestor for some reassortments
end

struct RootNode <: AbstractARGNode
	# children :: Vector{ARGNode}
end

# Equality
isequal(a::ARGNode, b::ARGNode) = (a.label==b.label)
isequal(a::ARGNode, b::RootNode) = false
isequal(a::RootNode, b::ARGNode) = false
isequal(a::ARGNode, b::Nothing) = false
isequal(a::Nothing, b::ARGNode) = false
Base.:(==)(x::ARGNode, y::ARGNode) = isequal(x,y)
Base.hash(x::ARGNode, h::UInt) = hash(x.label, h)

# Access
ancestors(n::ARGNode) = n.anc
function ancestor(n::ARGNode, c::Int)
	@assert hascolor(n, c)
	return n.anc[c]
end

children(n::ARGNode) = n.children
children(n::RootNode) = n.children
"""
	children(n::ARGNode, clr)

Children of `n` for color `clr`. Child `c` must meet the conditions:
- `hascolor(c, clr)`
- `ancestor(c, clr) == n`
"""
function children(n::ARGNode, clr)
	Iterators.filter(children(n)) do c
		hascolor(c, clr) && ancestor(c, clr) == n
	end
end

color(n::ARGNode) = n.color
color(::RootNode) = nothing

label(n::ARGNode) = n.label
label(::RootNode) = "_RootNode_"
label(::Nothing) = nothing

branch_length(n::ARGNode) = n.tau
branch_length(n::RootNode) = 0.
tau(n) = branch_length(n)
function branch_length(n::ARGNode, clr)
	@assert hascolor(n, clr)
	return n.tau[clr]
end
function branch_length(n::RootNode, clr)
	@assert hascolor(n, clr)
	return 0.
end


## isroot / isleaf
"""
	isroot(n)

Is `n` the root for one color?
"""
isroot(n::ARGNode) = in(RootNode(), n.anc)
"""
	isroot(n, c)

Is `n` the root for color `c`?
"""
isroot(n::ARGNode, c::Int) = (n.anc[c] == RootNode())
function is_global_root(n::ARGNode)
	if isroot(n)
		for c in 1:length(n.anc)
			if !isroot(n, c)
				return false
			end
		end
		return true
	end
	return false
end
function is_partial_root(n::ARGNode)
	if !isroot(n)
		return false
	else
		return !is_global_root(n)
	end
end


isleaf(n::ARGNode) = isempty(children(n))

ishybrid(n::ARGNode) = n.hybrid
ishybrid(::RootNode) = false

degree(n) = sum(color(n))
isshared(n) = (degree(n) == length(color(n)))

# Tests
hasancestor(n::ARGNode, c::Int) = !isnothing(ancestor(n,c)) && !isroot(n,c)
function hasancestor(n::ARGNode, a::ARGNode)
	for x in ancestors(n)
		if x == a
			return true
		end
	end
	return false
end

function haschild(a::ARGNode, n::ARGNode)
	for c in children(a)
		if c == n
			return true
		end
	end
	return false
end

hascolor(n::ARGNode, c::Int) = n.color[c]
hascolor(::RootNode, c::Int) = true
hascolor(::Nothing, c::Int) = true
function share_color(n1::ARGNode, n2::ARGNode)
	for (c1,c2) in zip(color(n1), color(n2))
		if c1 && c2
			return true
		end
	end
	return false
end
share_color(n::ARGNode, r::RootNode) = true
share_color(r::RootNode, n::ARGNode) = true


# Set / add children / ancestors
function add_child!(n::ARGNode, c::ARGNode)
	@assert !haschild(n, c) && share_color(n, c)
	push!(children(n), c)
end

function set_ancestor!(n, a, c::Int)
	@assert hascolor(n, c) && hascolor(a, c) "Nodes do not share color $c: \n $(label(n)): $(color(n)) \n $(label(a)): $(color(a))"
	n.anc[c] = a
end
function set_ancestor!(n, a)
	for (i,c) in color(n)
		if c
			set_ancestor!(n, a, i)
		end
	end
end

function add_color!(n, c::Int)
	@assert !hascolor(n, c) "Node $(n.label) already has color $c"
	n.color[c] = true
end

"""
	set_branch_length!(n, τ, c, clrs...)

Set branch above `n` to length `τ` for all colors in `(c, clrs...)`.
"""
function set_branch_length!(n, τ::Real, c, clrs...)
	@assert hascolor(n, c) "$(n.label) is not of color $c\n $n"
	@assert hasancestor(n, c) "No ancestor of color $c for $(n.label)\n $(n)"
	n.tau[c] = τ + eps()
	for clr in clrs
		set_branch_length!(n, τ, clr)
	end
end
function set_branch_length!(n, ::Missing, c, clrs...)
	@assert hascolor(n, c) "$(n.label) is not of color $c\n $n"
	n.tau[c] = missing
	for clr in clrs
		set_branch_length!(n, missing, clr)
	end
end

set_label!(n::ARGNode, label) = (n.label = label)

# Utils.
function othercolor(c)
	@assert (c==1 || c==2) "Invalid color $c"
	return (c==1) ? 2 : 1
end
"""
	edgecolor(a, n)

Color of the edge going from `a` to `n`
"""
function edgecolor(a, n)
	clr = [false, false]
	for (i,c) in enumerate(color(n))
		if c && a == ancestor(n, i)
			clr[i] = true
		end
	end
	return clr
end



mutable struct ARG
	roots :: Dict{Int, ARGNode}
	nodes :: Dict{String, ARGNode}
	hybrids :: Dict{String, ARGNode}
	leaves :: Dict{String, ARGNode}
end
roots(arg::ARG) = arg.roots
nodes(arg::ARG) = values(arg.nodes)
hybrids(arg::ARG) = values(arg.hybrids)
leaves(arg::ARG) = values(arg.leaves)
