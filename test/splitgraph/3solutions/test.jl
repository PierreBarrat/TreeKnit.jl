t1 = read_tree("$(dirname(pathof(TreeKnit)))/..//test/splitgraph/3solutions/t1.nwk")
t2 = read_tree("$(dirname(pathof(TreeKnit)))/..//test/splitgraph/3solutions/t2.nwk")

treelist = deepcopy(Any[t1, t2])

mcc = naive_mccs(treelist)
mcc_names = TreeKnit.name_mcc_clades!(treelist, mcc)
for (i,t) in enumerate(treelist)
	treelist[i] = TreeKnit.reduce_to_mcc(t, mcc)
end

g = trees2graph(treelist)

oa = OptArgs(resolve = false, likelihood_sort = false, γ = 2)
Trange = reverse(oa.Tmin:oa.dT:oa.Tmax)

@testset "No likelihood" begin
	local out = SplitGraph.sa_opt(g; Trange=Trange, γ=2, M = 100, resolve = false)[1]
	sort!(out)
	@test out[1] == Bool[0,1,1]
	@test out[2] == Bool[1,0,1]
	@test out[3] == Bool[1,1,0]
end


@testset "With likelihood" begin
	local out = []
	for rep in 1:20
		tmp = SplitGraph.opttrees!(
			deepcopy(t1),deepcopy(t2);
			γ=2, M=10, Trange=Trange, likelihood_sort = true
		)
		if !in(tmp[1], out)
			push!(out, tmp[1])
		end
	end
	@test out[1] == [["C1","C2"]]
end
