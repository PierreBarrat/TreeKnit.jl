export GraphNode, SplitNode, LeafNode

abstract type GraphNode end


mutable struct SplitNode <: GraphNode
	anc::Union{SplitNode,Nothing}
	child::Array{GraphNode,1}
	color::Int64
	conf::Array{Bool,1}
	isroot::Bool
end
function SplitNode( ; anc=nothing, child=Array{GraphNode,1}(undef, 0), color=0, conf=Array{Bool,1}(undef,0), isroot=true)
	return  SplitNode(anc, child, color, conf, isroot)
end

mutable struct LeafNode <: GraphNode
	anc::Array{SplitNode,1}
	conf::Array{Bool,1}
	index::Int64
end
function LeafNode(;anc=Array{SplitNode,1}(undef, 0), conf=Array{Bool,1}(undef,0), index=0)
	return LeafNode(anc, conf, index)
end

mutable struct Graph
	nodes::Array{GraphNode,1}
	leaves::Array{LeafNode,1}
	lleaves::Dict{String, LeafNode}
	internals::Array{SplitNode,1}
	labels::Array{String,1}
	labels_to_int::Dict{Any,Int64}
	K::Int64 
end
function Graph(; nodes=Array{GraphNode,1}(undef,0), leaves=Array{LeafNode,1}(undef,0), lleaves=Dict{String, LeafNode}(), internals=Array{SplitNode,1}(undef,0), labels=Array{String,1}(undef,0), labels_to_int=Dict{Any,Int64}(), K=0)
	return Graph(nodes, leaves, lleaves, internals, labels, labels_to_int, K)
end


