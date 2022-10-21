using Test
using TreeTools
using TreeKnit, TreeKnit.SplitGraph

nwk1 = "((1:1,2:1):1,((3:1,4:1):1,5:1,6:1):1)"
nwk2 = "((6:1,2:1):1,((3:1,4:1):1,5:1,1:1):1)"

t1 = node2tree(TreeTools.parse_newick(nwk1), label = "a")
t2 = node2tree(TreeTools.parse_newick(nwk2), label = "b")

solutions = [
	[["1"], ["2"], ["3", "4", "5", "6"]],
	[["2"], ["6"], ["1", "3", "4", "5"]],
	[["1"], ["6"], ["2", "3", "4", "5"]],
	[["1"], ["2", "6"], ["3", "4", "5"]],
	[["6"], ["1", "2"], ["3", "4", "5"]],
]
@testset "Dyn-resolve" begin
	MCCs = TreeKnit.computeMCCs(
		t1, t2, OptArgs(likelihood_sort=false, γ=3);
	)
	@test in(MCCs, solutions)
end

@testset "tree naming" begin
	nwk = ["~/Documents/TreeKnit/TreeKnit.jl/examples/tree_h3n2_ha.nwk", "~/Documents/TreeKnit/TreeKnit.jl/examples/tree_h3n2_na.nwk"]
	fn = TreeKnit.get_tree_names(nwk)
	@test fn == ["tree_h3n2_ha.nwk", "tree_h3n2_na.nwk"]
	@test TreeKnit.make_output_tree_names(fn) == ["tree_h3n2_ha_resolved.nwk", "tree_h3n2_na_resolved.nwk"]

	nwk = ["~/Documents/TreeKnit/TreeKnit.jl/examples1/tree.nwk", "~/Documents/TreeKnit/TreeKnit.jl/examples2/tree.nwk"]
	fn = TreeKnit.get_tree_names(nwk)
	@test fn == ["tree_examples1.nwk", "tree_examples2.nwk"]
	@test TreeKnit.make_output_tree_names(fn) == ["tree_examples1_resolved.nwk", "tree_examples2_resolved.nwk"]
end
