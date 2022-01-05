using Test
using TreeTools
using TreeKnit, TreeKnit.SplitGraph

println("##### main #####")

nwk1 = "((1:1,2:1):1,((3:1,4:1):1,5:1,6:1):1)"
nwk2 = "((6:1,2:1):1,((3:1,4:1):1,5:1,1:1):1)"

t1 = node2tree(TreeTools.parse_newick(nwk1))
t2 = node2tree(TreeTools.parse_newick(nwk2))

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
