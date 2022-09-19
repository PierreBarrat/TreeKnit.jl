using Test
using TreeTools
using TreeKnit, TreeKnit.SplitGraph

println("##### main #####")

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

println("### IO ###")

outdir = "test_auspice"
if !isdir(outdir)
	mkdir(outdir)
end
t1, t2, rS = TreeKnit.resolve_strict(t1, t2, solutions[1])
TreeTools.ladderize!(t1)
TreeKnit.sort_polytomies_strict!(t1, t2, solutions[1])
TreeKnit.write_auspice_json(outdir * "/", t1, t2, solutions[1])
TreeKnit.write_newick(outdir * "/" * t1.label *".nwk", t1)
TreeKnit.write_newick(outdir * "/" * t2.label*".nwk", t2)
