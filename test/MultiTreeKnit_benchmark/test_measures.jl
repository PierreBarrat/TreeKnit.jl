using Test
using TreeTools
using TreeKnit

println("##### measures #####")

nwk1 = "((A,B),C);"
nwk2 = "(A,(B,C));"
nwk3 = "((A,B,D),C);"
nwk4 = "(((A,B),D),C);"

t1 = node2tree(TreeTools.parse_newick(nwk1))
t2 = node2tree(TreeTools.parse_newick(nwk2))
t3 = node2tree(TreeTools.parse_newick(nwk3))
t4 = node2tree(TreeTools.parse_newick(nwk4))

@testset "RF distance" begin
	@test TreeKnit.RF_distance(t1, t2) ==2
    @test TreeKnit.RF_distance(t3, t4) ==1
end

MCC1 = [["1"], ["2"], ["3", "4", "5", "6"]]
MCC2 = [["1", "2"], ["3", "4", "5", "6"]]
MCC3 = [["1", "2", "3", "4", "5", "6"]]

@testset "is_degenerate" begin
	@test TreeKnit.is_degenerate(3, [MCC1, MCC1, MCC2]) ==false
    @test TreeKnit.is_degenerate(3, [MCC1, MCC1, MCC3]) ==false
    @test TreeKnit.is_degenerate(3, [MCC2, MCC2, MCC1]) ==true
end

@testset "consistency_rate" begin
    nwk1 = "((A,B),C);"
    nwk2 = "(A,(B,C));"
    nwk3 = "(A,B,C);"
    t1 = node2tree(TreeTools.parse_newick(nwk1))
    t2 = node2tree(TreeTools.parse_newick(nwk2))
    t3 = node2tree(TreeTools.parse_newick(nwk3))
    input_trees = [copy(t1), copy(t2), copy(t3)]
    tree_names = ["a", "b", "c"]
    MCC_dict, trees = TreeKnit.infer_benchmark_MCCs(input_trees, tree_names)
    c = TreeKnit.consistency_rate(MCC_dict[["a", "b"]], MCC_dict[["a", "c"]], MCC_dict[["b", "c"]], trees)
    cfull = TreeKnit.consistency_rate([MCC_dict[["a", "b"]], MCC_dict[["a", "c"]], MCC_dict[["b", "c"]]], trees)
	@test c == 0
    @test c == cfull
    MCC12 = [["A"], ["B", "C"]]
    MCC13 = [["A", "B", "C"]]
    MCC23 = [["B"], ["B", "C"]]
    MCC14 = [["A", "B", "C"]]
    MCC24 = [["B"], ["B", "C"]]
    MCC34 = [["A", "B", "C"]]
    c1 = TreeKnit.consistent_mcc_triplets([MCC12, MCC13, MCC23], [t1, t2, t3])
    @test c1 ==1/3
    c = TreeKnit.consistency_rate(MCC12, MCC13, MCC23, [t1, t2, t3])
    @test c ==1/9
    c4 = TreeKnit.consistency_rate([MCC12, MCC13, MCC14, MCC23, MCC24, MCC34], [t1, t2, t3, t3])
    @test c4 < c
end
