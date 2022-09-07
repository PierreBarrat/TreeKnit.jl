using TreeKnit
using TreeKnit.SplitGraph
using Test
using TreeTools


MCC_ab = [["5_0"], ["10_0", "6_0"], ["1_0", "2_0", "3_0", "4_0", "7_0", "8_0", "9_0"]]
MCC_ac = [["1_0", "4_0", "9_0"], ["10_0", "2_0", "3_0", "5_0", "6_0", "7_0", "8_0"]]
MCC_bc = [["5_0"], ["10_0", "6_0"], ["2_0", "7_0"], ["3_0", "8_0"], ["1_0", "4_0", "9_0"]]

nwk_a = "((8_0:1.0,3_0:1.0)internal_1:1.0,((9_0:1.0,4_0:1.0,1_0:1.0)RESOLVED_1:1.0,((10_0:1.0,6_0:1.0)RESOLVED_3:1.0,(5_0:1.0,(7_0:1.0,2_0:1.0)internal_4:1.0)internal_8:1.0)RESOLVED_2:1.0)internal_23:1.0)internal_26:1.0;"
nwk_b = "((5_0:1.0,(10_0:1.0,6_0:1.0)RESOLVED_4:1.0)internal_5:1.0,((8_0:1.0,3_0:1.0)RESOLVED_1:1.0,((7_0:1.0,2_0:1.0)RESOLVED_2:1.0,(9_0:1.0,4_0:1.0,1_0:1.0)RESOLVED_5:1.0)RESOLVED_3:1.0)internal_24:1.0)internal_27:1.0;"
nwk_c = "((9_0:1.0,4_0:1.0,1_0:1.0)internal_23:1.0,((8_0:1.0,3_0:1.0)internal_1:1.0,((10_0:1.0,6_0:1.0)RESOLVED_3:1.0,(5_0:1.0,(7_0:1.0,2_0:1.0)RESOLVED_1:1.0)RESOLVED_2:1.0)internal_11:1.0)internal_20:1.0)internal_26:1.0;" 
t_a_output = node2tree(TreeTools.parse_newick(nwk_a, node_data_type=TreeTools.MiscData), label = "a")
t_b_output = node2tree(TreeTools.parse_newick(nwk_b, node_data_type=TreeTools.MiscData), label = "b")
t_c_output = node2tree(TreeTools.parse_newick(nwk_c, node_data_type=TreeTools.MiscData), label = "c")
input_iMCCs = TreeKnit.convert_MCC_list_to_set(3, ["a", "b", "c"], [MCC_ab, MCC_ac, MCC_bc])
input_trees = [t_a_output, t_b_output, t_c_output]

@testset "fix_consist_sets! split test, topo=true" begin
    for i in range(1,10)
        trees = [copy(t) for t in input_trees]
        iMCCs = copy(input_iMCCs)
        new_sets = TreeKnit.fix_consist_sets!(iMCCs, trees; topo=true)
        changed_ab, changed_ac = new_sets.mccs[Set(["b", "a"])], new_sets.mccs[Set(["c", "a"])] 
        option1 = [["5_0"], ["10_0", "6_0"], ["2_0", "7_0"], ["1_0", "3_0", "4_0", "8_0", "9_0"]], [["3_0", "8_0"], ["1_0", "4_0", "9_0"], ["10_0", "2_0", "5_0", "6_0", "7_0"]]
        option1b = [["5_0"], ["10_0", "6_0"], ["2_0", "7_0"], ["3_0", "8_0"], ["1_0", "4_0", "9_0"]], [["1_0", "4_0", "9_0"], ["10_0", "2_0", "3_0", "5_0", "6_0", "7_0", "8_0"]]
        option2 = [["5_0"], ["10_0", "6_0"], ["1_0", "2_0", "3_0", "4_0", "7_0", "8_0", "9_0"]], [["2_0", "7_0"], ["3_0", "8_0"], ["1_0", "4_0", "9_0"], ["10_0", "5_0", "6_0"]]
        option2b = [["5_0"], ["10_0", "6_0"], ["3_0", "8_0"], ["1_0", "2_0", "4_0", "7_0", "9_0"]], [["2_0", "7_0"], ["1_0", "4_0", "9_0"], ["10_0", "3_0", "5_0", "6_0", "8_0"]]
        @test (changed_ab, changed_ac) == option1 || (changed_ab, changed_ac) == option2 || (changed_ab, changed_ac) == option1b || (changed_ab, changed_ac) == option2b
        @test sum(TreeKnit.consistency_rate(TreeKnit.iter_pairs(new_sets)[2]..., trees)) ==0 
        @test TreeKnit.is_topologically_degenerate(iMCCs, trees) == false
    end
end
@testset "fix_consist_sets! split test, topo=false" begin
    for i in range(1,10)
        trees = [copy(t) for t in input_trees]
        iMCCs = copy(input_iMCCs)
        new_sets = TreeKnit.fix_consist_sets!(iMCCs, trees; topo=false)
        changed_ab, changed_ac = new_sets.mccs[Set(["b", "a"])], new_sets.mccs[Set(["c", "a"])] 
        option1 = [["5_0"], ["10_0", "6_0"], ["2_0", "7_0"], ["1_0", "3_0", "4_0", "8_0", "9_0"]], [["1_0", "4_0", "9_0"], ["10_0", "2_0", "3_0", "5_0", "6_0", "7_0", "8_0"]]
        option2 = [["5_0"], ["10_0", "6_0"], ["1_0", "2_0", "3_0", "4_0", "7_0", "8_0", "9_0"]], [["2_0", "7_0"], ["1_0", "4_0", "9_0"], ["10_0", "3_0", "5_0", "6_0", "8_0"]]
        @test (changed_ab== option1[1] && changed_ac == option1[2] || changed_ab== option2[1] && changed_ac == option2[2])
        @test (sum(TreeKnit.consistency_rate(TreeKnit.iter_pairs(new_sets)[2]..., trees)) ==0)
        @test (TreeKnit.is_topologically_degenerate(iMCCs, trees))
    end
end
@testset "fix_consist_sets! always makes MCCs consistent for 3 trees" begin
    for i in range(1,10)
        true_trees, arg = TreeKnit.get_trees(2, 100, remove=false, c=0.75, ρ = 0.1);
        rMCCs = TreeKnit.get_real_MCCs(2, arg)
        trees = TreeKnit.remove_branches(true_trees; c=0.75)
        #TreeTools.print_tree_ascii(" ", trees[1])
        realMCCs_trees = [copy(t) for t in trees]

        rS = TreeKnit.resolve!(realMCCs_trees[2], realMCCs_trees[1], rMCCs[1])
        print("Distance true trees from realMCCs trees:"*string(sum([TreeKnit.RF_distance(true_trees[i], realMCCs_trees[i]) for i in 1:2])/2))
        print("Distance true trees from unresolved trees:"*string(sum([TreeKnit.RF_distance(true_trees[i], trees[i]) for i in 1:2])/2))
        prior_splits = [length(SplitList(t)) for t in trees]
        iMCCs = TreeKnit.get_infered_MCC_pairs!(trees, consistent = false, force_consist=false, constraint_cost=4., rounds=1)
        print("Distance true trees from infered trees:"*string(sum([TreeKnit.RF_distance(true_trees[i], trees[i]) for i in 1:2])/2))
        post_splits = [length(SplitList(t)) for t in trees]
        post_no_MCCs = [length(mcc[2]) for mcc in iMCCs.mccs]
        iMCCs = TreeKnit.fix_consist_sets!(iMCCs, trees; topo = false)
        post_fix_consist_no_MCCs = [length(mcc[2]) for mcc in iMCCs.mccs]
        print("Splits added by sequential:"*string(sum(post_splits - prior_splits)/2)*"\n")
        print("Additional MCC splits added by fix consist:"*string(sum(post_fix_consist_no_MCCs - post_no_MCCs)/2)*"\n")
        print("Additional MCC splits to fix topological issues:"*string(TreeKnit.count_new_splits_for_topo_degeneracy(iMCCs, trees)/2)*"\n")
        #@test sum(TreeKnit.consistency_rate(iMCCs, trees)) == 0
    end
end

@testset "fix_consist_sets! always makes MCCs consistent for 4 trees" begin
    for i in range(1,10)
        true_trees, arg = TreeKnit.get_trees(4, 100, remove=false, c=0.75, ρ = 0.1);
        trees = TreeKnit.remove_branches(true_trees; c=0.75)
        print("Distance true trees from input trees:"*string(sum([TreeKnit.RF_distance(true_trees[i], trees[i]) for i in 1:4])/4))
        prior_splits = [length(SplitList(t)) for t in trees]
        iMCCs = TreeKnit.get_infered_MCC_pairs!(trees, consistent = false, force_consist=false, constraint_cost=4.)
        print("Distance true trees from infered trees:"*string(sum([TreeKnit.RF_distance(true_trees[i], trees[i]) for i in 1:4])/4))
        post_splits = [length(SplitList(t)) for t in trees]
        post_no_MCCs = [length(mcc[2]) for mcc in iMCCs.mccs]
        iMCCs = TreeKnit.fix_consist_sets!(iMCCs, trees; topo = false)
        post_fix_consist_no_MCCs = [length(mcc[2]) for mcc in iMCCs.mccs]
        print("Splits added by sequential:"*string(sum(post_splits - prior_splits))*"\n")
        print("Additional MCC splits added by fix consist:"*string(sum(post_fix_consist_no_MCCs - post_no_MCCs))*"\n")
        print("Additional MCC splits to fix topological issues:"*string(TreeKnit.count_new_splits_for_topo_degeneracy(iMCCs, trees))*"\n")
        #@test sum(TreeKnit.consistency_rate(iMCCs, trees)) == 0
    end
end