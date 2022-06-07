using TreeKnit
using TreeKnit.SplitGraph
using Test
using TreeTools
using Combinatorics
t1 = read_tree("$(dirname(pathof(TreeKnit)))/..//test/splitgraph/3solutions/t1.nwk")
t2 = read_tree("$(dirname(pathof(TreeKnit)))/..//test/splitgraph/3solutions/t2.nwk")

t3 = node2tree(parse_newick("(A,B,(C,D)internal_2)"))
t4 = node2tree(parse_newick("(A,B,C,D)"))
t5 = node2tree(parse_newick("((A,C)internal_2,B,D)"))


@testset "test naive mcc for multitrees" begin
	@test naive_mccs([t3, t3], Set(keys(t3.lleaves))) == [["A", "B", "C", "D"]]
	##check works with mask
	@test naive_mccs([t3, t3], Set(["A", "B", "internal_2"])) == [["A", "B", "internal_2"]]
	##check works after being resolved (in this case should not be resolved)
	resolve!(t3, t4, t5)
	@test naive_mccs([t3, t4], Set(keys(t3.lleaves))) == [["A"], ["B"], ["C"], ["D"]]
	@test naive_mccs([t3, t5], Set(keys(t3.lleaves))) == [["A"], ["B"], ["C"], ["D"]]
	@test naive_mccs([t4, t5], Set(keys(t3.lleaves))) == [["A"], ["B"], ["C"], ["D"]]
end

treelist = [t1, t2]
treelist, copy_leaves = TreeKnit.prepare_copies!(treelist);
mcc = Dict()
for tree_pair in Combinatorics.combinations(1:length(treelist), 2)
	mcc[tree_pair] = naive_mccs([treelist[tree_pair[1]], treelist[tree_pair[2]]], copy_leaves[tree_pair])
end
TreeKnit.name_mcc_clades!(treelist, copy_leaves, mcc)
for t in treelist
	TreeKnit.remove_zero_copies!(t)
end
g = trees2graph(treelist)
oa = OptArgs(resolve = false, likelihood_sort = false, γ = 2)

@testset "No likelihood" begin
	local out = SplitGraph.sa_opt(g; Trange=oa.Trange, γ=2, M = 100, resolve = false)[1]
	sort!(out)
	@test out[1] == Bool[0,1,1]
	@test out[2] == Bool[1,0,1]
	@test out[3] == Bool[1,1,0]
end

@testset "With likelihood" begin
	local out = []
	for rep in 1:20
		tmp = SplitGraph.opttrees(
			t1, t2;
			γ=2, M=10, Trange=oa.Trange, likelihood_sort = true
		)
		if !in(tmp[1], out)
			push!(out, tmp[1])
		end
	end
	@test out[1] == [["C1","C2"]]
end