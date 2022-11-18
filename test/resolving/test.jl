### A Pluto.jl notebook ###
# v0.12.21

# ╔═╡ 2ea1ccbc-8e0d-11eb-29c3-05a91ea18a6f
using TreeKnit
using TreeKnit.SplitGraph
using Test
using TreeTools


println("### Basic resolving")

t1 = node2tree(TreeTools.parse_newick("(((A,B),C),D)"))
t2 = node2tree(TreeTools.parse_newick("(B,C,(A,D))"))
@testset "1" begin
	newsplits = resolve!(t1,t2)
	@test length(newsplits[1]) == 0
	@test length(newsplits[2]) == 0
end

t1 = read_tree("$(dirname(pathof(TreeKnit)))/..//test/resolving/tree1.nwk")
t2 = read_tree("$(dirname(pathof(TreeKnit)))/..//test/resolving/tree2.nwk")
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

println("### Resolve using pre-inferred MCCs")
t3 = read_tree("$(dirname(pathof(TreeKnit)))/..//test/resolving/tree3.nwk")
t4 = read_tree("$(dirname(pathof(TreeKnit)))/..//test/resolving/tree4.nwk")
mccs = read_mccs("$(dirname(pathof(TreeKnit)))/..//test/resolving/mccs34.dat")
@testset "Resolve with pre-inferred MCCs" begin
	rS = resolve!(t3, t4, mccs)
	@test rS[1] == [["F", "G"]]
	@test rS[2] == [["E", "F", "G"]]
end

##### map_split_to_tree functions

nwk1 = "(((A1,A2),(B1,B2),(C1,C2)),D,E)"
t = node2tree(TreeTools.parse_newick(nwk1)) #
S = SplitList(t)

# Let's pretend we found (A1,A2,B1,B2) and (C1,C2,D) to be MCCs
Smcc = SplitList(S.leaves)
append!(Smcc.splits, [Split([1,2,3,4]), Split([5,6,7])])


@testset "`map_split_to_tree` functions" begin
	Smapped = TreeKnit.map_splits_to_tree!(Smcc, t)
	@test leaves(Smapped, 1) == ["A1", "A2", "B1", "B2"]
	@test leaves(Smapped, 2) == ["A1", "A2", "B1", "B2", "C1", "C2", "D"]
	@test TreeTools.iscompatible(Smapped[1], S)
	@test TreeTools.iscompatible(Smapped[2], S)
end

println("### Strict resolving")

nwk1 = "((A,B,C,D),E)"
nwk2 = "(((A,B),C),(D,E))"
input_t1 = node2tree(TreeTools.parse_newick(nwk1; node_data_type=TreeTools.MiscData))
input_t2 = node2tree(TreeTools.parse_newick(nwk2; node_data_type=TreeTools.MiscData))

@testset "strict `map_split_to_tree` functions, Example 1" begin
	t1 = copy(input_t1)
	t2 = copy(input_t2)
	MCCs = TreeKnit.sort([["A", "B", "C"], ["D", "E"]], lt=TreeKnit.clt)
	new_splits = TreeKnit.new_splits(t1, MCCs, t2; strict=false)
	new_splits_strict = TreeKnit.new_splits(t1, MCCs, t2; strict=true)
	@test new_splits_strict ==[["A", "B"], ["A", "B", "C"]]
	@test new_splits == new_splits_strict

	t1 = copy(input_t1)
	t2 = copy(input_t2)
	MCCs = TreeKnit.sort([["A", "B", "C"], ["D"], ["E"]], lt=TreeKnit.clt)
	new_splits = TreeKnit.new_splits(t1, MCCs, t2; strict=false)
	new_splits_strict = TreeKnit.new_splits(t1, MCCs, t2; strict=true)
	@test isempty(new_splits_strict)
	@test new_splits ==[["A", "B"], ["A", "B", "C"]]
end

##assume true trees
true_nwk1 = "((((A,B),(C,D)),E),F)"
true_nwk2 = "((((A,B),(C,D)),E),F)"
true_t2 = node2tree(TreeTools.parse_newick(true_nwk2; node_data_type=TreeTools.MiscData))
MCCs = TreeKnit.sort([["A", "B", "E", "F"], ["C", "D"]], lt=TreeKnit.clt)

@testset "strict `map_split_to_tree` functions, Example 2" begin
	unresolved_nwk1 = "((((A,B),C,D),E),F)"
	unresolved_t1 = node2tree(TreeTools.parse_newick(unresolved_nwk1; node_data_type=TreeTools.MiscData))
	new_splits = TreeKnit.new_splits(unresolved_t1, MCCs, true_t2; strict=false)
	new_splits_strict = TreeKnit.new_splits(unresolved_t1, MCCs, true_t2; strict=true)
	@test new_splits_strict ==[["C", "D"]]
	@test new_splits ==[["C", "D"]]

	unresolved_nwk1 = "(((A,B),C,D,E),F)"
	unresolved_t1 = node2tree(TreeTools.parse_newick(unresolved_nwk1; node_data_type=TreeTools.MiscData))
	new_splits = TreeKnit.new_splits(unresolved_t1, MCCs, true_t2; strict=false)
	new_splits_strict = TreeKnit.new_splits(unresolved_t1, MCCs, true_t2; strict=true)
	@test new_splits_strict ==[["C", "D"]]
	@test new_splits ==[["C", "D"], ["A", "B", "E"]]

	unresolved_nwk1 = "((A,B),C,D,E,F)"
	unresolved_t1 = node2tree(TreeTools.parse_newick(unresolved_nwk1; node_data_type=TreeTools.MiscData))
	new_splits = TreeKnit.new_splits(unresolved_t1, MCCs, true_t2; strict=false)
	new_splits_strict = TreeKnit.new_splits(unresolved_t1, MCCs, true_t2; strict=true)
	@test isempty(new_splits_strict)
	@test new_splits  ==[["C", "D"], ["A", "B", "E"], ["A", "B", "E", "F"]]
end

##assume true trees
true_nwk1 = "(((3,2),(1,7)),(5,(4,6)))"
true_nwk2 = "((((3,2),(4,6)),(7,1)),5)"
unresolved_nwk2 = "(((3,2),(4,6)),(7,1),5)"
true_t1 = node2tree(TreeTools.parse_newick(true_nwk1; node_data_type=TreeTools.MiscData))
MCCs = TreeKnit.sort([["5"],["4"],["6"],["1", "2", "3", "7"]], lt=TreeKnit.clt)

@testset "strict `map_split_to_tree` functions, sister.mcc==nothing and anc.mcc==nothing" begin
	unres_t2 = node2tree(TreeTools.parse_newick(unresolved_nwk2; node_data_type=TreeTools.MiscData))
	new_splits_strict = TreeKnit.new_splits(unres_t2, MCCs, true_t1; strict=true)
	@test isempty(new_splits_strict) 
end

nwk_1 = "((A,B,C,D,E),X)"
nwk_2 = "((((A,B),C),(D,E)),X)"
t1 = node2tree(TreeTools.parse_newick(nwk_1; node_data_type=TreeTools.MiscData); label="t1")
t2 = node2tree(TreeTools.parse_newick(nwk_2; node_data_type=TreeTools.MiscData); label="t2")
MCCs = TreeKnit.sort([["A","B","C"],["D","E","X"]], lt=TreeKnit.clt)

@testset "strict_resolve - do not add unneeded splits" begin
	new_splits = TreeKnit.new_splits(t1, MCCs, t2; strict=false)
	new_splits_strict = TreeKnit.new_splits(t1, MCCs, t2; strict=true) 
	@test new_splits == [["A", "B"], ["A", "B", "C"], ["D", "E"]]
	@test new_splits_strict == [["A", "B"], ["A", "B", "C"]]

	t1_copy = copy(t1)
	t2_copy = copy(t2)
	rS = resolve!(t1_copy, t2_copy, MCCs; tau = 0.)
	@test write_newick(t1_copy) == "((((A,B)RESOLVED_1:0.0,C)RESOLVED_2:0.0,(D,E)RESOLVED_3:0.0)NODE_2,X)NODE_1:0;"
	
	t1_copy = copy(t1)
	t2_copy = copy(t2)
	rS = TreeKnit.resolve!(t1_copy, t2_copy, MCCs; tau = 0., strict=true)
	@test write_newick(t1_copy) == "((D,E,((A,B)RESOLVED_1:0.0,C)RESOLVED_2:0.0)NODE_2,X)NODE_1:0;"
	
	##check ladderize works the same for strict resolve
	TreeTools.ladderize!(t1_copy)
	@test write_newick(t1_copy) =="(X,(D,E,(C,(A,B)RESOLVED_1:0.0)RESOLVED_2:0.0)NODE_2)NODE_1:0;"
end

nwk_1 = "((A,B,C,(D,E)),X)"
nwk_2 = "((((A,B),C),(D,E)),X)"
t1 = node2tree(TreeTools.parse_newick(nwk_1; node_data_type=TreeTools.MiscData); label="t1")
t2 = node2tree(TreeTools.parse_newick(nwk_2; node_data_type=TreeTools.MiscData); label="t2")
t1_empty = node2tree(TreeTools.parse_newick(nwk_1; node_data_type=TreeTools.EmptyData); label="t1")
t2_empty = node2tree(TreeTools.parse_newick(nwk_2; node_data_type=TreeTools.EmptyData); label="t2")
MCCs = TreeKnit.sort([["A","B","C"],["D","E","X"]], lt=TreeKnit.clt)

@testset "check sort_polytomies! acts the same on strictly resolved trees" begin
	t1_copy = copy(t1_empty)
	t2_copy = copy(t2_empty)
	rS = resolve!(t1_copy, t2_copy, MCCs; tau = 0.)
	TreeTools.ladderize!(t1_copy)
	TreeKnit.sort_polytomies!(t1_copy, t2_copy, MCCs; strict=false)
	
	t1_strict = copy(t1)
	t2_strict = copy(t2)
	rS_strict = TreeKnit.resolve!(t1_strict, t2_strict, MCCs; tau = 0., strict=true)
	TreeTools.ladderize!(t1_strict)
	TreeKnit.sort_polytomies!(t1_strict, t2_strict, MCCs; strict=true)

	@test rS == rS_strict
	@test write_newick(t1_copy) == write_newick(t1_strict)
	@test write_newick(t2_copy) == write_newick(t2_strict)
end

nwk_1 = "(A,(C,D,E,B))"
nwk_2 = "(A,(E,(C,D),B))"
t1 = node2tree(TreeTools.parse_newick(nwk_1; node_data_type=TreeTools.MiscData); label="t1")
t2 = node2tree(TreeTools.parse_newick(nwk_2; node_data_type=TreeTools.MiscData); label="t2")
t1_empty = node2tree(TreeTools.parse_newick(nwk_1; node_data_type=TreeTools.EmptyData); label="t1")
t2_empty = node2tree(TreeTools.parse_newick(nwk_2; node_data_type=TreeTools.EmptyData); label="t2")
MCCs = TreeKnit.sort([["B"], ["A","C","D","E"]], lt=TreeKnit.clt)

function get_leaf_order(t)
	t_leaves =String[]
	for leaf in POTleaves(t)
		push!(t_leaves, leaf.label)
	end
	return t_leaves
end

@testset "check sort_polytomies! works when internal node could belong to more than one mcc" begin
	t1_strict = copy(t1)
	t2_strict = copy(t2)
	rS_strict = TreeKnit.resolve!(t1_strict, t2_strict, MCCs; tau = 0., strict=true)
	TreeTools.ladderize!(t1_strict)
	TreeKnit.sort_polytomies!(t1_strict, t2_strict, MCCs; strict=true)
	@test filter(e -> e!="B", get_leaf_order(t1_strict)) == filter(e -> e!="B", get_leaf_order(t2_strict))

	t1_strict = copy(t1)
	t2_strict = copy(t2)
	rS_strict = TreeKnit.resolve!(t2_strict, t1_strict, MCCs; tau = 0., strict=true)
	TreeTools.ladderize!(t2_strict)
	TreeKnit.sort_polytomies!(t2_strict, t1_strict, MCCs; strict=true)
	@test filter(e -> e!="B", get_leaf_order(t1_strict)) == filter(e -> e!="B", get_leaf_order(t2_strict))
end

nwk_1 = "(E,(B,C,D,A),K)"
nwk_2 = "(A,E,B,C,D,K)"
t1 = node2tree(TreeTools.parse_newick(nwk_1; node_data_type=TreeTools.MiscData); label="t1")
t2 = node2tree(TreeTools.parse_newick(nwk_2; node_data_type=TreeTools.MiscData); label="t2")
MCCs = TreeKnit.sort([["A","B","C","D","E"], ["K"]], lt=TreeKnit.clt)

@testset "check sort_polytomies! strict works on strictly resolved trees, 1" begin
	t1_strict = copy(t1)
	t2_strict = copy(t2)
	rS_strict = TreeKnit.resolve!(t1_strict, t2_strict, MCCs; tau = 0., strict=true)
	TreeTools.ladderize!(t2_strict)
	TreeKnit.sort_polytomies!(t2_strict, t1_strict, MCCs; strict=true)
	@test filter(e -> e!="K", get_leaf_order(t1_strict)) == filter(e -> e!="K", get_leaf_order(t2_strict))
end

t1 = node2tree(TreeTools.parse_newick("((A,B1),B2,C,D)"))
t2 = node2tree(TreeTools.parse_newick("((A,B1,B2,D),C)"))
MCCs = [["D"], ["A", "B1", "B2", "C"]]
rS_strict = TreeKnit.resolve!(t1, t2, MCCs; tau = 0., strict=true)
TreeTools.ladderize!(t1)
TreeKnit.sort_polytomies!(t1, t2, MCCs; strict=true)
@testset "check sort_polytomies! strict works on strictly resolved trees, 2" begin
	@test filter(e -> e!="D", get_leaf_order(t1)) == filter(e -> e!="D", get_leaf_order(t2))
end

nwk_1 = "(A,B,C,D)"
nwk_2 = "(A,C,(D,B))"
t1 = node2tree(TreeTools.parse_newick(nwk_1; node_data_type=TreeTools.MiscData); label="t1")
t2 = node2tree(TreeTools.parse_newick(nwk_2; node_data_type=TreeTools.MiscData); label="t2")
MCCs = TreeKnit.sort([["B"], ["A","C","D"]], lt=TreeKnit.clt)

@testset "check sort_polytomies! works when internal node could belong to more than one mcc" begin
	t1_strict = copy(t1)
	t2_strict = copy(t2)
	rS_strict = TreeKnit.resolve!(t1_strict, t2_strict, MCCs; tau = 0., strict=true)
	TreeTools.ladderize!(t1_strict)
	TreeKnit.sort_polytomies!(t1_strict, t2_strict, MCCs; strict=true)
end