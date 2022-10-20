using Test
using TreeTools
using TreeKnit

MCCs_ref = [
	["A/NewYork/105/2002"],
	["A/NewYork/177/1999"],
	["A/NewYork/137/1999", "A/NewYork/138/1999"],
	["A/NewYork/52/2004", "A/NewYork/59/2003"],
	["A/NewYork/198/2003", "A/NewYork/199/2003", "A/NewYork/32/2003"],
	["A/NewYork/10/2004", "A/NewYork/11/2003", "A/NewYork/12/2003", "A/NewYork/13/2003", "A/NewYork/14/2003", "A/NewYork/15/2003", "A/NewYork/16/2003", "A/NewYork/17/2003", "A/NewYork/18/2003", "A/NewYork/19/2003", "A/NewYork/2/2003", "A/NewYork/21/2003", "A/NewYork/22/2003", "A/NewYork/23/2003", "A/NewYork/24/2003", "A/NewYork/25/2003", "A/NewYork/26/2003", "A/NewYork/27/2003", "A/NewYork/28/2003", "A/NewYork/29/2003", "A/NewYork/30/2003", "A/NewYork/31/2004", "A/NewYork/33/2004", "A/NewYork/34/2003", "A/NewYork/35/2003", "A/NewYork/36/2003", "A/NewYork/38/2003", "A/NewYork/39/2003", "A/NewYork/4/2003", "A/NewYork/40/2003", "A/NewYork/41/2003", "A/NewYork/42/2003", "A/NewYork/43/2003", "A/NewYork/44/2003", "A/NewYork/45/2003", "A/NewYork/46/2003", "A/NewYork/47/2003", "A/NewYork/48/2003", "A/NewYork/49/2003", "A/NewYork/5/2004", "A/NewYork/50/2003", "A/NewYork/51/2003", "A/NewYork/53/2003", "A/NewYork/54/2003", "A/NewYork/55/2003", "A/NewYork/56/2003", "A/NewYork/6/2004", "A/NewYork/60A/2003", "A/NewYork/61A/2003", "A/NewYork/62A/2003", "A/NewYork/63/2003", "A/NewYork/64/2003", "A/NewYork/65/2003", "A/NewYork/67/2003", "A/NewYork/69/2004", "A/NewYork/7/2003", "A/NewYork/70/2004", "A/NewYork/8/2003"],
]

t1 = read_tree("$(dirname(pathof(TreeKnit)))/../test/NYdata/tree_ha.nwk")
t2 = read_tree("$(dirname(pathof(TreeKnit)))/../test/NYdata/tree_na.nwk")

MCCs = computeMCCs(t1, t2)

@testset "computeMCCs on NY data" begin
	@test MCCs[1:end-1] == MCCs_ref
end
if MCCs[1:end-1] != MCCs_ref
	@warn "Found different MCCs for the NewYork data. Could indicate a problem..."
end

"""
	function check_sort_polytomies(t1, t2, MCCs)
		
Check that leaves in the same MCC are in the same order in tree1 and tree2 after 
calling ladderize and sort_polytomies on the resolved trees. This will make sure that 
lines between nodes in an MCC do not cross when The two trees are visualized as a dendrogram.
"""
function check_sort_polytomies(t1, t2, MCCs) 
	t1, t2, rS_strict = TreeKnit.resolve_strict(t1, t2, MCCs; tau = 0.)
	TreeTools.ladderize!(t1)
	TreeKnit.sort_polytomies_strict!(t1, t2, MCCs)
	leaf_map = map_mccs(MCCs)
	pos_in_mcc_t1 = Dict()
	pos = 1
	for leaf in POTleaves(t1)
		if haskey(pos_in_mcc_t1, leaf_map[leaf.label])
			append!(pos_in_mcc_t1[leaf_map[leaf.label]], pos)
		else
			pos_in_mcc_t1[leaf_map[leaf.label]] = [pos]
		end
		pos +=1 
	end

	pos_in_mcc_t2 = Dict()
	pos = 1
	for leaf in POTleaves(t2)
		if haskey(pos_in_mcc_t2, leaf_map[leaf.label])
			append!(pos_in_mcc_t2[leaf_map[leaf.label]], pos)
		else
			pos_in_mcc_t2[leaf_map[leaf.label]] = [pos]
		end
		pos +=1 
	end

	sorted = true
	for mcc in 1:length(MCCs)
		if (sort(pos_in_mcc_t1[mcc]) !=  pos_in_mcc_t1[mcc]) || (sort(pos_in_mcc_t2[mcc]) !=  pos_in_mcc_t2[mcc])
			sorted = false
			break
		end
	end
	return sorted
end

@testset "sort_polytomies_strict! on resolve_strict! NY trees" begin
	@test check_sort_polytomies(t1, t2, MCCs)
end

