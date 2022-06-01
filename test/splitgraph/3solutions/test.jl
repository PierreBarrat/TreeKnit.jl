t1 = read_tree("$(dirname(pathof(TreeKnit)))/..//test/splitgraph/3solutions/t1.nwk")
t2 = read_tree("$(dirname(pathof(TreeKnit)))/..//test/splitgraph/3solutions/t2.nwk")

t3 = node2tree(parse_newick("(A,B,(C,D)internal_2)"))
l = Dict{String, TreeNode{TreeTools.EmptyData}}()
for k in keys(t3.lnodes)
	if t3.lnodes[k].label in ["A", "B", "internal_2"]
		l[k] = t3.lnodes[k]
	end
end

@testset "test naive mcc for multitrees" begin
	@test naive_mccs([t3, t3], t3.lleaves) == [["A", "B", "C", "D"]]
	##check works with mask
	@test naive_mccs([t3, t3], l) == [["A", "B", "internal_2"]]
end

treelist = deepcopy(Tree{TreeTools.EmptyData}[t1, t2])

mcc = naive_mccs(treelist)
mcc_names = TreeKnit.name_mcc_clades!(treelist, mcc)
for (i,t) in enumerate(treelist)
	treelist[i] = TreeKnit.reduce_to_mcc(t, mcc)
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