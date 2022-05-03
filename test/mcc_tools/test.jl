### A Pluto.jl notebook ###
# v0.19.3

using Markdown
using InteractiveUtils

# ╔═╡ 055fada0-c7b0-11ec-0df1-35c1df17620f
begin
	using Revise
	using Pkg; Pkg.activate()
	using Test
	using TreeKnit
	using TreeTools
end

# ╔═╡ 37d26f44-b429-4412-b8b3-01065876efac
begin
	nwk1 = "(((((A,B),(R1,R2)),(S1,S2)),(Q1,Q2)),(P1,P2))"
	nwk2 = "((((R1,R2),(S1,S2)),((Q1,Q2),(A,B))),(P1,P2))"
	t1 = parse_newick_string(nwk1)
	t2 = parse_newick_string(nwk2)
end

# ╔═╡ b2d075d4-f776-4c96-a576-1114bfc07e2b
MCCs = computeMCCs(t1, t2)

# ╔═╡ 7320c38e-32db-4bea-9876-db161edeff06
@testset "Neighbour leaves" begin
	@test TreeKnit.neighbour_leaves(MCCs[1], t1) == ["R1", "R2"]
	@test TreeKnit.neighbour_leaves(MCCs[1], t2) == ["Q1", "Q2"]
	@test TreeKnit.neighbour_leaves(MCCs[1], t1, t2) == union(
		TreeKnit.neighbour_leaves(MCCs[1], t1),
		TreeKnit.neighbour_leaves(MCCs[1], t2)
	)
end

# ╔═╡ d0b8c7ec-7167-4de6-be75-7e46727e3d26
@testset "Neighbour joint nodes" begin
	mcc_maps = [TreeKnit.mcc_map(t1, MCCs), TreeKnit.mcc_map(t2, MCCs)]
	@test TreeKnit._neighbour_joint_nodes(MCCs[1], t1, mcc_maps[1]) == sort.([["R1", "R2"]])
	@test TreeKnit._neighbour_joint_nodes(MCCs[1], t2, mcc_maps[2]) == sort.([["Q1", "Q2"]])
end

# ╔═╡ Cell order:
# ╠═055fada0-c7b0-11ec-0df1-35c1df17620f
# ╠═37d26f44-b429-4412-b8b3-01065876efac
# ╠═b2d075d4-f776-4c96-a576-1114bfc07e2b
# ╠═7320c38e-32db-4bea-9876-db161edeff06
# ╠═d0b8c7ec-7167-4de6-be75-7e46727e3d26
