using TreeKnit
using TreeKnit.SplitGraph
using Test
using TreeTools


MCC_ab = [["5_0"], ["6_0", "10_0"], ["1_0", "2_0", "3_0", "4_0", "7_0", "8_0", "9_0"]]
MCC_ac = [["1_0", "4_0", "9_0"], ["10_0", "2_0", "3_0", "5_0", "6_0", "7_0", "8_0"]]
MCC_bc = [["5_0"], ["10_0", "6_0"], ["2_0", "7_0"], ["3_0", "8_0"], ["1_0", "4_0", "9_0"]]

nwk_a = "((8_0:1.0,3_0:1.0)internal_1:1.0,((9_0:1.0,4_0:1.0,1_0:1.0)RESOLVED_1:1.0,((10_0:1.0,6_0:1.0)RESOLVED_3:1.0,(5_0:1.0,(7_0:1.0,2_0:1.0)internal_4:1.0)internal_8:1.0)RESOLVED_2:1.0)internal_23:1.0)internal_26:1.0;"
nwk_b = "((5_0:1.0,(10_0:1.0,6_0:1.0)RESOLVED_4:1.0)internal_5:1.0,((8_0:1.0,3_0:1.0)RESOLVED_1:1.0,((7_0:1.0,2_0:1.0)RESOLVED_2:1.0,(9_0:1.0,4_0:1.0,1_0:1.0)RESOLVED_5:1.0)RESOLVED_3:1.0)internal_24:1.0)internal_27:1.0;"
nwk_c = "((9_0:1.0,4_0:1.0,1_0:1.0)internal_23:1.0,((8_0:1.0,3_0:1.0)internal_1:1.0,((10_0:1.0,6_0:1.0)RESOLVED_3:1.0,(5_0:1.0,(7_0:1.0,2_0:1.0)RESOLVED_1:1.0)RESOLVED_2:1.0)internal_11:1.0)internal_20:1.0)internal_26:1.0;" 
t_a_output = node2tree(TreeTools.parse_newick(nwk_a, node_data_type=TreeTools.MiscData), label = "a")
t_b_output = node2tree(TreeTools.parse_newick(nwk_b, node_data_type=TreeTools.MiscData), label = "b")
t_c_output = node2tree(TreeTools.parse_newick(nwk_c, node_data_type=TreeTools.MiscData), label = "c")
iMCCs = TreeKnit.convert_MCC_list_to_set(3, ["a", "b", "c"], [MCC_ab, MCC_ac, MCC_bc])
trees = [t_a_output, t_b_output, t_c_output]
@testset "fix_consist_sets! split test" begin
    for i in range(1,10)
        new_sets = TreeKnit.fix_consist_sets!(iMCCs, trees)
        changed_ab, changed_ac = new_sets.mccs[Set(["b", "a"])], new_sets.mccs[Set(["c", "a"])] 
        option1 = [["5_0"], ["10_0", "6_0"], ["2_0", "7_0"], ["1_0", "3_0", "4_0", "8_0", "9_0"]], [["1_0", "4_0", "9_0"], ["10_0", "2_0", "3_0", "5_0", "6_0", "7_0", "8_0"]]
        option2 = [["5_0"], ["6_0", "10_0"], ["1_0", "2_0", "3_0", "4_0", "7_0", "8_0", "9_0"]], [["2_0", "7_0"], ["1_0", "4_0", "9_0"], ["10_0", "3_0", "5_0", "6_0", "8_0"]]
        @test changed_ab== option1[1] && changed_ac == option1[2] || changed_ab== option2[1] && changed_ac == option2[2]
        @test sum(TreeKnit.consistency_rate(TreeKnit.iter_pairs(new_sets)[2]..., trees)) ==0 
    end
end
@testset "fix_consist_sets! always makes MCCs consistent for 3 trees" begin
    for i in range(1,10)
        trees, arg = TreeKnit.get_trees(3, 100, remove=true, c=0.2 + 0.6*rand(), ρ = 0.1);
        iMCCs = TreeKnit.get_infered_MCC_pairs!(trees, consistant = true, force=true, constraint_cost=4.)
        @test sum(TreeKnit.consistency_rate(iMCCs, trees)) == 0
    end
end

@testset "fix_consist_sets! always makes MCCs consistent for 4 trees" begin
    for i in range(1,10)
        trees, arg = TreeKnit.get_trees(4, 100, remove=true, c=0.2 + 0.6*rand(), ρ = 0.1);
        iMCCs = TreeKnit.get_infered_MCC_pairs!(trees, consistant = true, force=true, constraint_cost=4.)
        @test sum(TreeKnit.consistency_rate(iMCCs, trees)) == 0
    end
end