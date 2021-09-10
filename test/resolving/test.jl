### A Pluto.jl notebook ###
# v0.12.21

# ╔═╡ 2ea1ccbc-8e0d-11eb-29c3-05a91ea18a6f
using RecombTools
using RecombTools.SplitGraph
using Test
using TreeTools


println("##### Basic resolving #####")

t1 = node2tree(TreeTools.parse_newick("(((A,B),C),D)"))
t2 = node2tree(TreeTools.parse_newick("(B,C,(A,D))"))
@testset "1" begin
	newsplits = resolve!(t1,t2)
	@test length(newsplits[1]) == 0
	@test length(newsplits[2]) == 0
end

t1 = read_tree("resolving/tree1.nwk")
t2 = read_tree("resolving/tree2.nwk")
@testset "2" begin
	newsplits = resolve!(t1,t2)
	@test newsplits[1].splits == [Split([1,2])]
	@test newsplits[2] == [
		["A1", "A2", "B", "C"],
		["D1", "D2"],
		["E", "F", "G", "H"],
	]
end


t1 = node2tree(TreeTools.parse_newick("(A,B,C,D)"))
t2 = node2tree(TreeTools.parse_newick("(A,(B,C,D))"))
t3 = node2tree(TreeTools.parse_newick("(A,B,(C,D))"))

@testset "Node naming" begin
	ns12 = resolve!(t1,t2)
	ns13 = resolve!(t1,t3)
	@test isempty(ns12[2])
	@test isempty(ns12[2])
	@test haskey(t1.lnodes, "RESOLVED_1")
	@test haskey(t1.lnodes, "RESOLVED_2")
	@test !haskey(t1.lnodes, "RESOLVED_3")
end


println("##### Resolve using MCCs #####")

# ╔═╡ 93ccb338-8e0d-11eb-0f36-49bebaf5570b
function infer_mccs(trees;kwargs...)
	RecombTools.runopt(
		OptArgs(;seq_lengths = Dict(s=>1 for s in keys(trees)), kwargs...),
		trees
	)
end

# ╔═╡ 45e695a6-8e0d-11eb-178d-97bf0cb45379
nwk1 = "(A:1.,B:1.,(C:1.,D:1.):1.)"

# ╔═╡ 6165e6e0-8e0d-11eb-1f87-b94b13b3ec95
nwk2 = "(((A:1.,B:1.):1.,C:1.):1.,D:100.)"

# ╔═╡ 6f0f4310-8e0d-11eb-1110-9d6d12d0cbae
t1 = node2tree(TreeTools.parse_newick(nwk1))

# ╔═╡ 7bdceb1a-8e0d-11eb-0160-0d60f81fd810
t2= node2tree(TreeTools.parse_newick(nwk2))

# ╔═╡ 90c5c30c-8e0e-11eb-00a0-45ccc42a2500
trees = Dict(1=>deepcopy(t1), 2=>deepcopy(t2))

# ╔═╡ fa1be92c-8e0d-11eb-0f5a-21209d7ffc68
@testset "Inner function _resolve_from_mccs!" begin
	ct = deepcopy(trees)
	nS = RecombTools._resolve_from_mccs!(infer_mccs, ct)
	@test isempty(nS[2])
	@test nS[1] == [["A", "B"]]
end

# ╔═╡ ef1e0a34-8e10-11eb-198d-1d31d25d8167
nwk3 = "(((A1:1,A2:1):1.,A3:1,A4:1):1.,B1:1.,B2:1,(C:1.,D:1.):1)"

# ╔═╡ 117d3fe6-8e11-11eb-02bf-6b21357867d7
t3 = node2tree(TreeTools.parse_newick(nwk3))

# ╔═╡ 1be36e60-8e11-11eb-2ca0-e59c20fa35a4
nwk4 = "(((((A1:1,A2:1,A3:1):1,A4:1):1.,(B1:1.,B2:1):1):1,C:1.):1.,D:100.)"

# ╔═╡ 50b097f8-8e11-11eb-1a43-ed66b32e8ee4
t4 = node2tree(TreeTools.parse_newick(nwk4))

# ╔═╡ 5fb91e64-8e11-11eb-3ac7-4906f631bd80
trees2 = Dict(1=>deepcopy(t3), 2=>deepcopy(t4))

# ╔═╡ 8d4d48b4-8e11-11eb-2445-f5196c3666b3
@testset "Outer function resolve_from_mccs!" begin
	ct2 = deepcopy(trees2)
	nS2 = RecombTools.resolve_from_mccs!(infer_mccs, ct2)
	@test nS2[2] == [["A1", "A2"]]
	@test nS2[1] == [["A1", "A2", "A3"], ["B1", "B2"], ["A1", "A2", "A3", "A4", "B1", "B2"]]
end


# Testing max-clique splits
# trees 2 and 3 try to resolve 1 in a different way.
nwk1 = "((A,B,C),(D,E,F,G))"
nwk2 = "(((A,B),C),(D,(E,(F,G))))"
nwk3 = "(((A,B),C),(D,(E,F),G))"
trees = Dict(i=>node2tree(parse_newick(nwk)) for (i,nwk) in enumerate([nwk1, nwk2, nwk3]))
MCCs = RecombTools._computeMCCs(infer_mccs, trees)
resolvable_splits = RecombTools.new_splits(trees, MCCs)
@testset "max-clique splits" begin
	S1 = RecombTools.max_clique_splits(resolvable_splits[1])
	S2 = RecombTools.max_clique_splits(resolvable_splits[2])
	S3 = RecombTools.max_clique_splits(resolvable_splits[3])
	@test length(S1) == 3 # would be 0 with compat
	@test length(S2) == 0 # is already resolved
	@test length(S3) == 1 # one split from tree 2 (depending if E or G is an MCC)
end


# ╔═╡ Cell order:
# ╠═2ea1ccbc-8e0d-11eb-29c3-05a91ea18a6f
# ╠═93ccb338-8e0d-11eb-0f36-49bebaf5570b
# ╠═45e695a6-8e0d-11eb-178d-97bf0cb45379
# ╠═6165e6e0-8e0d-11eb-1f87-b94b13b3ec95
# ╠═6f0f4310-8e0d-11eb-1110-9d6d12d0cbae
# ╠═7bdceb1a-8e0d-11eb-0160-0d60f81fd810
# ╠═90c5c30c-8e0e-11eb-00a0-45ccc42a2500
# ╠═fa1be92c-8e0d-11eb-0f5a-21209d7ffc68
# ╠═ef1e0a34-8e10-11eb-198d-1d31d25d8167
# ╠═117d3fe6-8e11-11eb-02bf-6b21357867d7
# ╠═1be36e60-8e11-11eb-2ca0-e59c20fa35a4
# ╠═50b097f8-8e11-11eb-1a43-ed66b32e8ee4
# ╠═5fb91e64-8e11-11eb-3ac7-4906f631bd80
# ╠═8d4d48b4-8e11-11eb-2445-f5196c3666b3
