"""
	write_mccs(of, MCCs::AbstractArray)

Write MCCs to file.
"""
function write_mccs(of, MCCs::AbstractArray, mode="w")
	open(of, mode) do w
		for (i,m) in enumerate(MCCs)
			for x in m[1:end-1]
				write(w, x * ",")
			end
			if i < length(MCCs)
				write(w, m[end] * "\n")
			else
				write(w, m[end])
			end
		end
	end

	return nothing
end


function read_mccs(file)
	map(eachline(file)) do m
       String.(split(m, ','))
    end
end


"""
function write_auspice_json(filepath, tree1::Tree{T}, tree2::Tree{T}, MCCs::Vector{Vector{String}}) where T 
Returns a .json file that can be used by auspice to plot a dendrogram of two trees 
with their MCCs drawn on in color. The format used is:
{
'nodes': {
	'node1': {branch_length: 0.0343, mcc_tree2: 1},
	'node2': {...}
}}
"""
function write_auspice_json(filepath, tree1::Tree{T}, tree2::Tree{T}, MCCs::Vector{Vector{String}})where T 
	trees = [copy(t) for t in [tree1, tree2]]
	for tree in trees
		other_tree = [t for t in trees if t.label!=tree.label][1]
		other_tree_name = other_tree.label
		assign_mccs!(tree, leaf_mcc_map(MCCs))
		full_filepath = filepath*"branch_lengths_"*tree.label*".json"
		start = true
		open(full_filepath, "w") do w
			write(w, "{ \"nodes\" : {\n")
			for n in values(tree.lnodes)
				if ismissing(n.tau)
					n.tau = 0.0
				end
				if !start
					write(w, ", \n")
				else
					start = false
				end
				write(w, "\""*n.label*"\": {\"branch_length\": "*string(n.tau)*", ")
				if !isnothing(n.data["mcc"])
					write(w, "\"mcc_"*other_tree_name*"\": "*string(n.data["mcc"])*"")
				else
					write(w, "\"mcc_"*other_tree_name*"\": null")
				end
				write(w, "}")
			end
			write(w, "}}")
		end
	end
	return nothing
end