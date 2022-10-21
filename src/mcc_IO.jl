"""
	write_mccs(of, MCCs::AbstractArray)

Write MCCs to file. A dat file can only be produced for one tree pair MCC, each line in the dat file
will represent one mcc.

If multiple MCCs should be written to a file they need to be written to a JSON file. This will have the form:

{ 
    "MCC_dict" : {
        "1": { 
            "trees":["a", "b"],
            "mccs": [["A"],["B","C"]]
            },
        "2": { 
            "trees":["a", "c"],
            "mccs": [["A","B","C"]]
            },
        ...
    }
}
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
with their MCCs drawn on in color. The format used is specified in 
https://github.com/nextstrain/augur/blob/master/augur/data/schema-export-v2.json

for more information see: 
https://docs.nextstrain.org/projects/auspice/en/stable/releases/v2.html?highlight=augur%20export%20v2#new-dataset-json-format
"""
function write_auspice_json(filepath, trees::Vector{Tree{T}}, MCCs::MCC_set) where T 
	start_string_1 = """{ \"version\": \"v2\",
		\"meta\": {
		\"updated\": \"2022-10-18\",
		\"colorings\": [
		"""
	start_string_2 = """],
		\"filters\": [],
		\"panels\": [
			\"tree\"
		]
		},
		\"tree\":
		"""
	for tree in trees
		other_trees = [copy(trees[i]) for i in 1:MCCs.no_trees if trees[i].label!=tree.label]
		other_tree_names = [t.label for t in other_trees]
		mcc_maps = [map_mccs(tree, get(MCCs, tree.label, otree.label)) for otree in other_trees]
		full_filepath = filepath*"auspice__"*tree.label*".json"
		open(full_filepath, "w") do w
			write(w, start_string_1)
			l = length(other_tree_names)
			for ot in other_tree_names
				order_ = sort([tree.label, ot])
				write(w, " { \"key\": \"mcc_"*order_[1]*"_"*order_[2]*"\",")
				write(w, "\"title\": \"mcc_"*order_[1]*"_"*order_[2]*"\",")
				write(w, "\"type\": \"ordinal\" }")
				l>1 ? write(w, ",") : nothing
				l -= 1
			end
			write(w, start_string_2)
			write_node_date!(w, tree.root, mcc_maps, tree.label, other_tree_names)
			write(w, "}")
		end
	end
	return nothing
end


function write_node_date!(w, n, mcc_maps, tree_name, other_tree_names; div =0)
	write(w, "{\"name\": \""*string(n.label)*"\",")
	if ismissing(n.tau) n.tau = 0.0 end
	div = n.tau + div
	write(w, " \"node_attrs\": { \"div\": "*string(div)*", ")
	l = length(other_tree_names)
	for (i,ot) in enumerate(other_tree_names)
		order_ = sort([tree_name, ot])
		mcc_ = (!isnothing(mcc_maps[i][n.label])) ? string(mcc_maps[i][n.label]) : "null"
		write(w, " \"mcc_"*order_[1]*"_"*order_[2]*"\": { \"value\":"*mcc_*"}" )
		l>1 ? write(w, ",") : nothing
		l -= 1
	end
	write(w, "}, \"branch_attrs\": {}")
	if !isleaf(n)
		write(w, ", \"children\": [")
		l = length(n.child)
		for c in n.child
			write_node_date!(w, c, mcc_maps, tree_name, other_tree_names; div=div)
			l>1 ? write(w, ",") : nothing
			l -= 1
		end
		write(w, "]")
	end
	write(w, "}")

	return nothing
end
	