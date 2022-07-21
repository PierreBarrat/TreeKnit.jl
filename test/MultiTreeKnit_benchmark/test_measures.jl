using Test
using TreeTools
using TreeKnit

println("##### measures #####")

nwk1 = "((A,B),C);"
nwk2 = "(A,(B,C));"
nwk3 = "((A,B,D),C);"
nwk4 = "(((A,B),D),C);"

t1 = node2tree(TreeTools.parse_newick(nwk1), label = "a")
t2 = node2tree(TreeTools.parse_newick(nwk2), label = "b")
t3 = node2tree(TreeTools.parse_newick(nwk3), label = "c")
t4 = node2tree(TreeTools.parse_newick(nwk4), label = "d")

@testset "RF distance" begin
	@test TreeKnit.RF_distance(t1, t2) ==2
    @test TreeKnit.RF_distance(t3, t4) ==1
end

MCC1 = [["1"], ["2"], ["3", "4", "5", "6"]]
MCC2 = [["1", "2"], ["3", "4", "5", "6"]]
MCC3 = [["1", "2", "3", "4", "5", "6"]]
MCC_dict = Dict(Set(["a", "b"]) => MCC1, Set(["a", "c"]) => MCC2, Set(["b", "c"]) => MCC2)

@testset "is_degenerate" begin
	@test TreeKnit.is_degenerate(MCC1, MCC1, MCC2) ==false
    @test TreeKnit.is_degenerate(MCC1, MCC1, MCC3) ==false
    @test TreeKnit.is_degenerate(MCC2, MCC2, MCC1) ==true
    @test TreeKnit.is_degenerate(MCC_dict) == true
end

@testset "consistency_rate" begin
    nwk1 = "((A,B),C);"
    nwk2 = "(A,(B,C));"
    nwk3 = "(A,B,C);"
    t1 = node2tree(TreeTools.parse_newick(nwk1), label = "a")
    t2 = node2tree(TreeTools.parse_newick(nwk2), label= "b")
    t3 = node2tree(TreeTools.parse_newick(nwk3), label= "c")
    input_trees = [deepcopy(t1), deepcopy(t2), deepcopy(t3)]
    MCC_dict = TreeKnit.infer_benchmark_MCCs!(input_trees, consistant=true)
    c = TreeKnit.consistency_rate(MCC_dict, input_trees)
    cfull = TreeKnit.consistency_rate(MCC_dict[Set(["a", "b"])], MCC_dict[Set(["a", "c"])], MCC_dict[Set(["c", "b"])], input_trees)
    @test c == 0
    @test c == cfull
    MCC12 = [["A"], ["B", "C"]]
    MCC13 = [["A", "B", "C"]]
    MCC23 = [["B"], ["A", "C"]]
    MCC14 = [["A", "B", "C"]]
    MCC24 = [["B"], ["A", "C"]]
    MCC34 = [["A", "B", "C"]]
    MCC_test_dict = Dict(Set(["a", "b"]) => MCC12, Set(["a", "c"]) => MCC13, Set(["a", "d"]) => MCC14, 
                    Set(["b", "c"]) => MCC23, Set(["b", "d"]) => MCC24, Set(["c", "d"]) => MCC34)
    c1 = TreeKnit.consistent_mcc_triplets([MCC12, MCC13, MCC23], [t1, t2, t3])
    @test c1 ==2/3 #Of the three nodes (B, C and NODE_1) in t1 that are in an MCC in MCC12 and MCC13, 2 are not in an MCC in MCC23 (B and NODE_1)
    c = TreeKnit.consistency_rate(MCC12, MCC13, MCC23, [t1, t2, t3])
    @test c ==(2/3 + 1/2)/3
    t4 = copy(t3)
    label!(t4, "d")
    c4 = TreeKnit.consistency_rate(MCC_test_dict, [t1, t2, t3, t4])
    @test c4 < c
end
