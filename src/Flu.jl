module Flu

using RecombTools
using TreeTools

#
export nt_seq_lengths
export read_flu_trees

# Constants

const flu_lineages = ["h3n2"]
const flu_segments = ["ha", "na", "pb1", ]

const nt_seq_lengths = Dict(
	"h3n2" => Dict(
		"ha" => 1698,
		"na" => 1410,
		"pb1" => 2311,
	),
)

const substitution_rates = Dict(
	"h3n2" => Dict(
		"ha" => 3.66e-3,
		"na" => 2.87e-3,
	)
)

#
"""
	read_flu_trees(files::Dict, lineage = "h3n2")

Read trees stored in `values(files)` and return them as a dictionary indexed with
  `keys(files)`.
  Short branches are removed from trees.
"""
function read_flu_trees(files::AbstractDict, lineage = "h3n2")
	if !in(lineage, flu_lineages)
		error("Unknown flu lineage $lineage. Options are $(flu_lineages)")
	end

	trees = Dict{String, Tree}()
	for (segment, f) in files
		if !in(segment, flu_segments)
			error("Unknown flu segment $segment. Options are $(flu_segments)")
		end
		trees[segment] = read_tree(f)
		L = nt_seq_lengths[lineage][segment]
		TreeTools.delete_null_branches!(trees[segment]; threshold = 1/L/2)
	end

	return trees
end


end # module

