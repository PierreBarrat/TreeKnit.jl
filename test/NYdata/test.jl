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
standard_trees = [copy(t1), copy(t2)]
parallel_trees = [copy(t1), copy(t2)]
t1_original = copy(t1)
t2_original = copy(t2)

MCCs = run_treeknit!(t1, t2, TreeKnit.OptArgs(rounds=1, pre_resolve=false))
MCC_pair = MCCs.mccs[Set([t1.label, t2.label])]

@testset "run_treeknit! on NY data" begin
	@test MCC_pair[1:end-1] == MCCs_ref
end
if MCC_pair[1:end-1] != MCCs_ref
	@warn "Found different MCCs for the NewYork data. Could indicate a problem..."
end

"""
	function check_sort_polytomies(t1, t2, MCCs)
		
Check that leaves in the same MCC are in the same order in tree1 and tree2 after 
calling ladderize and sort_polytomies on the resolved trees. This will make sure that 
lines between nodes in an MCC do not cross when the two trees are visualized as a tanglegram.
"""
function check_sort_polytomies(t1, t2, MCCs) 
	leaf_map = map_mccs(MCCs) ##map from leaf to mcc
	pos_in_mcc_t1 = Dict()
	for leaf in POTleaves(t1)
		if haskey(pos_in_mcc_t1, leaf_map[leaf.label])
			push!(pos_in_mcc_t1[leaf_map[leaf.label]], leaf.label)
		else
			pos_in_mcc_t1[leaf_map[leaf.label]] = [leaf.label]
		end
	end

	pos_in_mcc_t2 = Dict()
	for leaf in POTleaves(t2)
		if haskey(pos_in_mcc_t2, leaf_map[leaf.label])
			push!(pos_in_mcc_t2[leaf_map[leaf.label]], leaf.label)
		else
			pos_in_mcc_t2[leaf_map[leaf.label]] = [leaf.label]
		end
	end

	sorted = true
	for mcc in 1:length(MCCs)
		if (pos_in_mcc_t1[mcc] != pos_in_mcc_t2[mcc])
			sorted = false
			break
		end
	end
	return sorted
end

rS_strict = TreeKnit.resolve!(t1, t2, MCC_pair; tau = 0., strict=true)
TreeTools.ladderize!(t1)
TreeKnit.sort_polytomies!(t1, t2, MCC_pair; strict=true)
@testset "sort_polytomies! on strict resolve! NY trees" begin
	@test check_sort_polytomies(t1, t2, MCC_pair)
end

MCCs_standard = TreeKnit.run_standard_treeknit!(standard_trees, TreeKnit.OptArgs(rounds=1))
@testset "run_standard_treeknit! produces the same MCCs as run_treeknit! without preresolve" begin
	@test MCC_pair == get(MCCs_standard, (1,2))
	@test SplitList(t1) == SplitList(standard_trees[1])
end
@testset "run_standard_treeknit! correctly sorts polytomies" begin
	@test check_sort_polytomies(standard_trees[1], standard_trees[2], get(MCCs_standard, (1,2)))
end

MCCs_parallel = TreeKnit.run_parallel_treeknit!(parallel_trees, TreeKnit.OptArgs(rounds=1))
@testset "run_parallel_treeknit! produces the same MCCs as run_treeknit! without preresolve" begin
	@test MCC_pair == get(MCCs_parallel, (1,2))
	@test SplitList(t1) == SplitList(parallel_trees[1])
end
@testset "run_parallel_treeknit! correctly sorts polytomies" begin
	@test check_sort_polytomies(parallel_trees[1], parallel_trees[2], get(MCCs_parallel, (1,2)))
end

