"""
	write_mccs(of, MCCs::AbstractArray)

Write MCCs to file.
"""
function write_mccs(filePath, MCCs::AbstractArray, mode="w")
	file_type = split(filePath, ".")[2]
	if file_type !="dat"
		print("File type not supported")
	end
	open(filePath, mode) do w
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


function write_mccs(filePath, MCCs::MCC_set, mode="w")
	file_type = split(filePath, ".")[2]
	if file_type !="json"
		print("File type not supported")
	end

	open(filePath, mode) do w
		write(w, "{ \"MCC_dict\" : {\n")
		tree_pairs, Ms = iter_pairs(MCCs)
		for i in 1:binomial(MCCs.no_trees, 2)
			write(w, "\""*string(i)*"\": {\n \"trees\":[")
			s = sort(tree_pairs[i], lt=clt)
			write(w, "\""*split(s[1], ".")[1]*"\", \""*split(s[2], ".")[1]*"\"],")
			write(w, "\n\"mccs\": [")
			m = Ms[i]
			for x in m[1:end-1]
				write(w, string(x) * ",\n")
			end
			if i < length(Ms)
				write(w, string(m[end]) *"]\n},\n")
			else
				write(w, string(m[end])*"]\n}\n}\n}")
			end
		end
	end

	return nothing
end

"""
function write_auspice_json(filepath, trees::Vector{Tree{T}}, MCCs::MCC_set) where T 

Returns a .json file that can be used by auspice to plot a dendrogram of two trees 
with their MCCs drawn on in color. The format used is:

{
'nodes': {
	'node1': {branch_length: 0.0343, mcc_tree2: 1, mcc_tree3: nothing},
	'node2': {...}
}}
"""
function write_auspice_json(filepath, trees::Vector{Tree{T}}, MCCs::MCC_set) where T 
	for tree in trees
		other_trees = [trees[i] for i in 1:MCCs.no_trees if trees[i].label!=tree.label]
		other_tree_names = [t.label for t in other_trees]
		assign_all_mccs!(tree, other_trees, MCCs)
		full_filepath = filepath*"auspice_branch_lengths_"*tree.label*".json"
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
				if !isnothing(n.data["mcc"][1])
					write(w, "\"mcc_"*other_tree_names[1]*"\": "*string(n.data["mcc"][1])*"")
				else
					write(w, "\"mcc_"*other_tree_names[1]*"\": null")
				end
				for i in 2:(MCCs.no_trees-1)
					if !isnothing(n.data["mcc"][i])
						write(w, " ,\"mcc_"*other_tree_names[i]*"\": "*string(n.data["mcc"][i]))
					else
						write(w, "\"mcc_"*other_tree_names[i]*"\": null")
					end
				end
				write(w, "}")
			end
			write(w, "}}")
		end
	end
	return nothing
end

