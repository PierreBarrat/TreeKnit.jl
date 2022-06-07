using Test
using TreeTools
using TreeKnit, TreeKnit.SplitGraph
using Combinatorics

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
treelist = [t1, t2]
# treelist, copyleaves = TreeKnit.prepare_copies!(treelist);
# mcc = Dict()
# for tree_pair in Combinatorics.combinations(1:length(treelist), 2)
# 	mcc[tree_pair] = naive_mccs([treelist[tree_pair[1]], treelist[tree_pair[2]]], copyleaves[tree_pair])
# end
# mcc_names = TreeKnit.name_mcc_clades!(treelist, copyleaves, mcc)
# for t in treelist
# 	TreeKnit.remove_zero_copies!(t)
# end
# oa = OptArgs(resolve = false, likelihood_sort = false, γ = 2)
# g = trees2graph(treelist)
# oconfs, F, nfound = SplitGraph.sa_opt(trees2graph(treelist), γ=2, Trange=reverse(0.001:0.01:1.), M=10, rep=1, resolve=true)
# print("passed sa.opt")
# oconf, L = SplitGraph.sortconf(oconfs, treelist, g, oa.seq_lengths, mcc_names[[1,2]], oa.likelihood_sort, false)

# print("passed on their own\n")
# TreeKnit.runopt(oa, t1, t2; output = :mccs)
# print("runopt passed\n")

@testset "Dyn-resolve" begin
	MCCs = TreeKnit.computeMCCs(
		t1, t2, OptArgs(likelihood_sort=false, γ=3, verbose=true);
	)
	@test in(MCCs, solutions)
end
