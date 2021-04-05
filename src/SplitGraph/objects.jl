export GraphNode, SplitNode, LeafNode

abstract type GraphNode end

"""
	mutable struct SplitNode <: GraphNode

Represent a split in the tree of a given segment. It is an internal node in the graph of splits. 
- `color` is the color of the segment that `SplitNode` refers to
- `conf` defines the split. It is an array of `Bool` of length `N` (number of leaves). All leaves `i` such that `conf[i]` are on one side of the split, the others are on the other side. 
- `anc::Union{SplitNode,Nothing}`
- `child::Array{GraphNode,1}`
- `isroot::Bool`
"""
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

"""
	mutable struct LeafNode <: GraphNode

Represent a leaf in the graph of splits. 
- `anc::Array{SplitNode,1}`: One `SplitNode` ancestor for each segment. 
- `conf::Array{Bool,1}`: there is only one `i` for which `conf[i]`, by definition.
- `index::Int64`: `conf[index]==true`
"""
mutable struct LeafNode <: GraphNode
	anc::Array{SplitNode,1}
	conf::Array{Bool,1}
	index::Int64
end
function LeafNode(;anc=Array{SplitNode,1}(undef, 0), conf=Array{Bool,1}(undef,0), index=0)
	return LeafNode(anc, conf, index)
end

"""
	mutable struct Graph

Graph of genealogic splits for a set of sequences consisting of many segments. Above leaf level, this graph consists of `K` trees of splits, where `K` is the number of segments. 
"""
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


"""
	struct OptArgs

Storing parameters for `SplitGraph.runopt` function.

### General
- `γ::Real = 3`   
- `itmax::Int64 = 15`: Maximal number of iterations of MCC / SA cycles
- `likelihood_sort::Bool = true`: sort equivalent configurations using likelihood test (based on branch length for now). 
- `resolve::Bool = true`: try to resolve trees while finding MCCs. 
- `seq_lengths = ones(Int64, 2)`: lengths of sequences that trees were built from
- `crossmap::Bool = false`: Wether cross-mappped mutations should be used to prune MCCs preventively. If `runopt` is called with a dictionary `trees::Dict{<:Any,<:Tree}` as input, the code will look at the number of suspicious mutations at `n.data.dat[:suspicious_muts][s]` where `n = trees[s]` (`s` is assumed to be an influenza segment). 
### Simulated annealing
- `Md::Real = 10`:  Number of SA iterations (per temperature) for a tree of `n` leaves is `ceil(Int64, n/Md)`
- `Tmin::Float64 = 1e-3`: Minimal temperature of SA
- `Tmax::Float64 = 1`: Maximal temperature of SA
- `dT::Float64 = 1e-2`: Temperature step
### Verbosity
- `verbose::Bool=false`: first level of verbosity
- `vv::Bool = false`: second level of verbosity
### Output
- `output = :mccs`: possible values `[:mccs, :mccs_df, :all]`
"""
@with_kw struct OptArgs
	γ::Real  = 3
	itmax::Int64 = 15
	likelihood_sort::Bool = true
	resolve::Bool = true
	seq_lengths = ones(Int64, 2)
	crossmap::Bool = false
	# For the annealing  
	Md::Real = 10 
	Tmin::Float64 = 1e-3
	Tmax::Float64 = 1.; @assert Tmax > Tmin
	dT::Float64 = 1e-2
	Trange = reverse(Tmin:dT:Tmax)
	sa_rep::Int64 = 1
	# Verbosity
	verbose::Bool = false
	vv::Bool = false
	# Output
	output = :mccs
end
getM(n,Md) = ceil(Int64, n/Md)
