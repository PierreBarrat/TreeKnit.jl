"""
	write_mccs(filename, MCCs::AbstractArray[, mode="w"])
	write_mccs(filename, MCCs::MCC_set[, mode="w"])

Write MCCs to file.

In the first form, the writing format is csv with one MCC per line.
This only makes sense for two MCCs

The second form is also adapted to multiple MCCs.
The format is JSON:

```
{ 
    "MCC_dict" : {
        "1": { 
            "trees":["a", "b"],
            "mccs": [["A"],["B","C"]],
            },
        "2": { 
            "trees":["a", "c"],
            "mccs": [["A","B","C"]],
            },
        ...
    }
}
```
"""
function write_mccs(filename, MCCs::AbstractArray, mode="w")
	file_type = split(filename, ".")[end]
	open(filename, mode) do w
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

function write_mccs(filename, MCCs::MCC_set, mode="w")
	file_type = split(filename, ".")[end]
	if file_type !="json"
		@warn "Got filename $filename, but MCCs are written as a json format."
	end

	json_string = write_mccs(MCCs)

	open(filename, mode) do f
		JSON3.pretty(f, json_string)
	end
end

function write_mccs(MCCs::MCC_set)

	MCC_dict = Dict{Int, Any}()
	tree_pairs, Ms = iter_pairs(MCCs)
	for i in 1:binomial(MCCs.no_trees, 2)
		MCC_dict[i] = Dict{String, Any}("trees"=> tree_pairs[i], "mccs"=>Ms[i])
	end
	MCC_to_write = Dict{String, Any}("MCC_dict" => MCC_dict)
	json_string = JSON3.write(MCC_to_write)

	return json_string
end

"""
	write_auspice_json(filepath, trees::Vector{Tree}, MCCs::MCC_set)

Returns a .json file that can be used by auspice to plot a dendrogram of two trees 
with their MCCs drawn on in color. The format used is specified in 
https://github.com/nextstrain/augur/blob/master/augur/data/schema-export-v2.json

for more information see: 
https://docs.nextstrain.org/projects/auspice/en/stable/releases/v2.html?highlight=augur%20export%20v2#new-dataset-json-format
"""
function write_auspice_json(filepath, trees::Vector{Tree{T}}, MCCs::MCC_set) where T 

	for tree in trees
		full_filepath = filepath*"auspice_"*tree.label*".json"
		
		json_string = get_auspice_json(tree, trees, MCCs)
		
		open(full_filepath, "w") do f
			JSON3.pretty(f, json_string)
		end
	end
	return nothing
end

function get_auspice_json(tree, tree_list, MCCs)
	other_trees = [copy(tree_list[i]) for i in 1:MCCs.no_trees if tree_list[i].label!=tree.label]
	other_tree_names = [t.label for t in other_trees]
	mcc_maps = [map_mccs(tree, get(MCCs, tree.label, otree.label)) for otree in other_trees]
	
	coloring_dict_list =[]
	for ot in other_tree_names
		order_ = sort([tree.label, ot])
		coloring_dict = Dict{String, String}(
						"key" => "mcc_"*order_[1]*"_"*order_[2],
						"title" => "mcc_"*order_[1]*"_"*order_[2],
						"type" => "ordinal")
		push!(coloring_dict_list, coloring_dict)
	end
	meta_dict = Dict{String, Any}(
				"updated" => "", 
				"colorings" => coloring_dict_list, 
				"filters" => String[], 
				"panels" => ["tree"])
	tree_dict = get_auspice_tree(tree.root, mcc_maps, tree.label, other_tree_names; div =0)
	auspice_dict =  Dict{String, Any}("version" => "v2", "meta" => meta_dict, "tree" => tree_dict)
	json_string = JSON3.write(auspice_dict)

	return json_string
end

function get_auspice_tree(n, mcc_maps, tree_name, other_tree_names; div =0)
    if ismissing(n.tau) n.tau = 0.0 end
	div = n.tau + div
    node_attr_dict = Dict{String, Any}("div" => div)
    for (i,ot) in enumerate(other_tree_names)
		order_ = sort([tree_name, ot])
		mcc_ = (!isnothing(mcc_maps[i][n.label])) ? string(mcc_maps[i][n.label]) : "null"
		node_attr_dict["mcc_"*order_[1]*"_"*order_[2]] = Dict("value" => mcc_)
	end
    tree_dict = Dict{String, Any}(
        "name" => string(n.label), 
        "node_attrs" => node_attr_dict,
        "branch_attrs" => Dict())
    if !isleaf(n)
        tree_dict["children"] = []
        for c in n.child
            child_tree_dict = get_auspice_tree(c, mcc_maps, tree_name, other_tree_names; div =div)
			push!(tree_dict["children"], child_tree_dict)
		end
    end
    return tree_dict
end
	
