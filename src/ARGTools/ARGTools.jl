module ARGTools

using TreeTools

"""
	mutable struct ARGNode

An `ARG` is considered as a collage of `K` independent trees. 
Fields of the `ARGNode` structure: 
- `anc`: Array of `ARGNode` or `Nothing`. Each node in the ARG has between `1` and `K` ancestors. 
- `anc_color`: `Array{Array{Bool}}` of the same length as `anc`. For each ancestor, indicates the colors of the corresponding branch. 
- `children`
- `color`: Array of `Bool` of length `K`. It indicates the trees in which `ARGNode` exists. 
By construction, `sum(ARGNode.color)` should be equal to `sum(sum(ARGNode.anc_color))`. This means that each node has exactly one ancestor in each tree. For the same reason, `reduce((x,y)->x .+ y, n.anccolor)` should be all ones. 
- `degree`: length of `color`, number of trees current node exists. 
"""
mutable struct ARGNode
	anc::Array{Union{ARGNode,Nothing}}
	anccolor::Array{Array{Bool,1},1} # Should be renamed - `upbranchcolor` for instance. It's *not* the color of the ancestor, but the color of the branch leading to it. 
	children::Array{ARGNode}
	color::Array{Bool}
	degree::Int64
	label::String
	data::Array{TreeTools.EvoData}
	isroot::Array{Bool}
	isleaf::Bool
end
function ARGNode(; degree=1, 
	anc = Array{Union{ARGNode,Nothing}}(nothing, degree),
	anccolor = [_color(i, degree) for i in 1:length(anc)],
	children = Array{ARGNode}(undef, 0), 
	color = ones(Bool, degree),
	label = "",
	data = Array{TreeTools.EvoData}(undef, 0),
	isroot = ones(Bool, degree),
	isleaf = true
	)
	return ARGNode(anc, anccolor, children, color, degree, label, data, isroot, isleaf)
end

# Convenience functions
function _color(i::Int64, degree)
	out = zeros(Bool, degree)
	out[i] = true
	return out
end
function _color(idx::Array{Int64}, degree)
	out = zeros(Bool, degree)
	for i in idx
		out[i] = true
	end
	return out
end

"""
	mutable struct ARG

Structure containing a dictionary `nodes::Dict{String,ARGNode}` and a root node `root::Array{ARGNode}`. 
"""
mutable struct ARG
	degree::Int64
	root::Array{ARGNode}
	nodes::Dict{String,ARGNode}
end
function ARG(; degree=1,
	root = Array{ARGNode}(undef, degree),
	nodes = Dict{String, ARGNode}()
	)
	return ARG(degree, root, nodes)
end

include("tools.jl")


end

#=
Checks 
For n::ARGNode
reduce((x,y)->x .+ y, n.anccolor) should be all ones
sum(sum(n.anccolor)) == sum(n.color)

For A::ARG
All leaves should be of full color, *i.e.* `sum(n.color) == A.degree`
=#
