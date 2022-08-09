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
MCC_dict = MCC_set(3, ["a", "b", "c"], Dict(Set(["a", "b"]) => MCC1, Set(["a", "c"]) => MCC2, Set(["b", "c"]) => MCC2))

@testset "is_degenerate" begin
	@test TreeKnit.is_degenerate(MCC1, MCC1, MCC2) ==false
    @test TreeKnit.is_degenerate(MCC1, MCC1, MCC3) ==false
    @test TreeKnit.is_degenerate(MCC2, MCC2, MCC1) ==true
    @test TreeKnit.is_degenerate(MCC_dict) == true
end

@testset "accuracy consistency_rate, output TreeKnit" begin
    nwk1 = "((A,B),C);"
    nwk2 = "(A,(B,C));"
    nwk3 = "(A,B,C);"
    t1 = node2tree(TreeTools.parse_newick(nwk1, node_data_type=TreeTools.MiscData), label = "a")
    t2 = node2tree(TreeTools.parse_newick(nwk2, node_data_type=TreeTools.MiscData), label= "b")
    t3 = node2tree(TreeTools.parse_newick(nwk3, node_data_type=TreeTools.MiscData), label= "c")
    input_trees = [copy(t1), copy(t2), copy(t3)]
    MCC_dict = TreeKnit.infer_benchmark_MCCs!(input_trees, consistant=true)
    c = TreeKnit.consistency_rate(MCC_dict, input_trees)
    cfull = TreeKnit.consistency_rate(TreeKnit.iter_pairs(MCC_dict)[2]..., input_trees)
    @test c == 0
    @test c == sum(cfull)/3
    MCC12 = [["A"], ["B", "C"]]
    MCC13 = [["A", "B", "C"]]
    MCC23 = [["B"], ["A", "C"]]
    MCC14 = [["A", "B", "C"]]
    MCC24 = [["B"], ["A", "C"]]
    MCC34 = [["A", "B", "C"]]
    input_trees = [copy(t1), copy(t2), copy(t3)]
    t4 = copy(t3)
    label!(t4, "d")
    MCC_test_dict = MCC_set(4, ["a", "b", "c", "d"], Dict(Set(["a", "b"]) => MCC12, Set(["a", "c"]) => MCC13, Set(["a", "d"]) => MCC14, 
                    Set(["b", "c"]) => MCC23, Set(["b", "d"]) => MCC24, Set(["c", "d"]) => MCC34))
    c1 = TreeKnit.consistent_mcc_triplets([MCC12, MCC13, MCC23], [input_trees[2]])
    @test c1 ==1/2 #Of the 2 branches (between B, C and NODE_1) in t2 that are in an MCC in MCC12 and MCC13, 1 is not in MCC23 (between B and NODE_1)
    c = TreeKnit.consistency_rate(MCC12, MCC13, MCC23, input_trees)
    @test sum(c)/3 ==(1/2 + 1/2)/3
    c4 = TreeKnit.consistency_rate(MCC_test_dict, [input_trees..., t4])
    @test sum(c4)/3 < sum(c)/3
end

@testset "accuracy consistency_rate, unresolved trees" begin
    nwk_a = "(11,15,(4,12),(5,1),13,3,(7,10,14,2,6,8),9);"
    nwk_b = "((7,12,4,1,5),9,(13,3,10,14,2,6,8),11,15);"
    nwk_c = "(9,15,3,10,14,2,6,8,12,((13,1,5),(4,7,11)));"

    t_a = node2tree(TreeTools.parse_newick(nwk_a, node_data_type=TreeTools.MiscData), label = "a")
    t_b = node2tree(TreeTools.parse_newick(nwk_b, node_data_type=TreeTools.MiscData), label= "b")
    t_c = node2tree(TreeTools.parse_newick(nwk_c, node_data_type=TreeTools.MiscData), label= "c")

    ##real MCCs
    MCC_ab = [["11"], ["4"], ["7"], ["9"], ["1", "5"], ["10", "12", "13", "14", "15", "2", "3", "6", "8"]]
    MCC_ac = [["11"], ["12"], ["13"], ["4"], ["7"], ["9"], ["1", "5"], ["10", "14", "15", "2", "3", "6", "8"]]
    MCC_bc = [["11"], ["12"], ["13"], ["4"], ["7"], ["1", "5"], ["10", "14", "15", "2", "3", "6", "8", "9"]]
    MCC_dict = TreeKnit.convert_MCC_list_to_set(3, ["a", "b", "c"], [MCC_ab, MCC_ac, MCC_bc])
    c = TreeKnit.consistency_rate(MCC_dict, [t_a, t_b, t_c])
    @test c == 0
    @test TreeKnit.is_degenerate(MCC_dict) == false

    ##infered MCCs
    MCC_ab = [["7"], ["10", "11", "12", "13", "14", "15", "1", "2", "3", "4", "5", "6", "8", "9"]]
    MCC_ac = [["11"], ["12"], ["13"], ["4"], ["7"], ["1", "5"], ["2", "3", "6", "8", "9", "10", "15", "14"]]
    MCC_cb = [["11"], ["12"], ["13"], ["15"], ["4"], ["7"], ["9"], ["1", "5"], ["10", "14", "2", "3", "6", "8"]]
    iMCC_dict = TreeKnit.convert_MCC_list_to_set(3, ["a", "b", "c"], [MCC_ab, MCC_ac, MCC_cb])
    c = TreeKnit.consistency_rate(iMCC_dict, [t_a, t_b, t_c])
    @test c != 0
    @test TreeKnit.is_degenerate(iMCC_dict) == true
end

@testset "check is_degenerate(MCCs) == (consistency_rate!=0.0) holds" begin
    repeat = 0
    while repeat < 10
        trees, arg = TreeKnit.get_trees(3, 50, remove=true, c=0.75);
        label!(trees[1], "a")
        label!(trees[2], "b")
        label!(trees[3], "c")
        iMCCs = TreeKnit.get_infered_MCC_pairs!(trees, consistant = true, order="input", rev=false, constraint_cost=4.)
        @test TreeKnit.is_degenerate(iMCCs) == (TreeKnit.consistency_rate(iMCCs, trees)!=0.0)
        repeat += 1
    end
end

@testset "check is_degenerate(MCCs) == (consistency_rate!=0.0) holds after fix_consist!" begin
    repeat = 0
    while repeat < 10
        trees, arg = TreeKnit.get_trees(3, 50, remove=true, c=0.75);
        label!(trees[1], "a")
        label!(trees[2], "b")
        label!(trees[3], "c")
        iMCCs = TreeKnit.get_infered_MCC_pairs!(trees, consistant = true, order="input", rev=false, constraint_cost=4., force=true, force_rounds=3)
        @test TreeKnit.is_degenerate(iMCCs) == (TreeKnit.consistency_rate(iMCCs, trees)!=0.0)
        repeat += 1
    end
end