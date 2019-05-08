export CrossMutations

"""
	CrossMutations

Store information of mutations when mapping sequences of one segment onto the tree of another. 
- `crossmut`: dictionnary storing mutations. Keys are of the form `"seg_i","seg_j"` (*i.e.* `Tuple{String,String}`), with the __first__ segment indicating the tree, and the __second__ the sequences. Values are arrays of mutations. 
- `suspicious`: dictionnary indicating the number of suspicious mutations. Keys are of the same form as that of `crossmut`. Values indicate the number of  corresponding mutations in `crossmut` that are suspicious. 
"""
mutable struct CrossMutations
	crossmut::Dict{Tuple{String, String}, Array{Tuple{Int64, Int64, Int64}, 1}}
	suspicious::Dict{Tuple{String, String}, Int64}
end
function CrossMutations()
	return CrossMutations(Dict{Tuple{String, String}, Array{Tuple{Int64, Int64, Int64}, 1}}(), Dict{Tuple{String, String}, Int64}())
end