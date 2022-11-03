using Test
using Dagger
using TreeTools
using TreeKnit
using TreeKnit.MTK

println("#### test parallelized recursive MultiTreeKnit ###")

nwk1 = "((A,B),C)R"
nwk2 = "(A,(B,C))R"
nwk3 = "(A,B,C)R"
t1 = node2tree(TreeTools.parse_newick(nwk1, node_data_type=TreeTools.MiscData), label = "a")
t2 = node2tree(TreeTools.parse_newick(nwk2, node_data_type=TreeTools.MiscData), label = "b")
t3 = node2tree(TreeTools.parse_newick(nwk3, node_data_type=TreeTools.MiscData), label = "c")
pre_trees = [t1, t2, t3]

@testset "3 trees, 1 round" begin
    trees = [copy(t) for t in pre_trees]
    MCCs = MTK.parallelized_compute_mccs!(trees, OptArgs(rounds=1, consistent=true))
    @test SplitList(trees[3]) == SplitList(trees[1])
    @test get(MCCs, (1, 2)) == get(MCCs, (3, 2)) 
    
    trees = [copy(t) for t in pre_trees]
    MCCs_slow = MTK.compute_mcc_pairs!(trees, OptArgs(rounds=1, consistent=true))
    @test SplitList(trees[3]) == SplitList(trees[1])
    @test get(MCCs_slow, (1, 2)) == get(MCCs_slow, (3, 2)) 
end

@testset "3 trees, 2 rounds" begin
    trees = [copy(t) for t in pre_trees]
    MCCs = MTK.parallelized_compute_mccs!(trees, OptArgs(consistent=true))
    @test SplitList(trees[3]) == SplitList(trees[1])
    @test get(MCCs, (1, 2)) == get(MCCs, (3, 2)) 

    trees = [copy(t) for t in pre_trees]
    MCCs_slow = MTK.compute_mcc_pairs!(trees, OptArgs(consistent=true))
    @test SplitList(trees[3]) == SplitList(trees[1])
    @test get(MCCs_slow, (1, 2)) == get(MCCs_slow, (3, 2)) 
end

nwk1 = "(((A,B),C),(D,E,F))R"
nwk2 = "((A,(B,C)),(D,E,F))R"
nwk3 = "((A,B,C),((D,E),F))R"
nwk4 = "((A,B,C),(D,(E,F)))R"
nwk5 = "((A,(B,C)),(D,(E,F)))R"
t1 = node2tree(TreeTools.parse_newick(nwk1, node_data_type=TreeTools.MiscData), label = "a")
t2 = node2tree(TreeTools.parse_newick(nwk2, node_data_type=TreeTools.MiscData), label = "b")
t3 = node2tree(TreeTools.parse_newick(nwk3, node_data_type=TreeTools.MiscData), label = "c")
t4 = node2tree(TreeTools.parse_newick(nwk4, node_data_type=TreeTools.MiscData), label = "d")
t5 = node2tree(TreeTools.parse_newick(nwk5, node_data_type=TreeTools.MiscData), label = "e")
pre_trees = [t1, t2, t3, t4]

split_list_final_t2 = [["D", "E"], ["D", "E", "F"], ["B", "C"], ["A", "B", "C"], ["A", "B", "C", "D", "E", "F"]]
split_list_final_t4 = [["E", "F"], ["D", "E", "F"], ["A", "B"], ["A", "B", "C"], ["A", "B", "C", "D", "E", "F"]]

@testset "4 trees, 1 round" begin
    trees = [copy(t) for t in pre_trees]
    MCCs = MTK.parallelized_compute_mccs!(trees, OptArgs(rounds=1, consistent=true))
    @test get(MCCs, (1, 2)) == get(MCCs, (3, 2)) 
    joint_ = MTK.MCC_join_constraint([get(MCCs, (1, 2)),  get(MCCs, (1, 4))])
    @test get(MCCs, (1, 4)) == get(MCCs, (3, 4))
    @test get(MCCs, (2, 4)) == joint_
    @test SplitList(trees[3]) == SplitList(trees[1])
    @test SplitList(trees[2])== split_list_final_t2
    @test SplitList(trees[4])== split_list_final_t4

    trees = [copy(t) for t in pre_trees]
    MCCs_slow = MTK.compute_mcc_pairs!(trees, OptArgs(rounds=1, consistent=true))
    @test get(MCCs_slow , (1, 2)) == get(MCCs_slow , (3, 2)) 
    joint_ = MTK.MCC_join_constraint([get(MCCs_slow , (1, 2)),  get(MCCs_slow , (1, 4))])
    @test get(MCCs_slow , (1, 4)) == get(MCCs_slow , (3, 4))
    @test get(MCCs_slow , (2, 4)) == joint_
    @test SplitList(trees[3]) == SplitList(trees[1])
    @test SplitList(trees[2])== split_list_final_t2
    @test SplitList(trees[4])== split_list_final_t4
end

@testset "4 trees, 2 rounds" begin
    trees = [copy(t) for t in pre_trees]
    MCCs = MTK.parallelized_compute_mccs!(trees, OptArgs(consistent=true))
    @test get(MCCs, (1, 2)) == get(MCCs, (3, 2)) 
    joint_ = MTK.MCC_join_constraint([get(MCCs, (1, 2)),  get(MCCs, (1, 4))])
    @test get(MCCs, (1, 4)) == get(MCCs, (3, 4))
    @test get(MCCs, (2, 4)) == joint_
    @test SplitList(trees[3]) == SplitList(trees[1])
    @test SplitList(trees[2])== split_list_final_t2
    @test SplitList(trees[4])== split_list_final_t4

    trees = [copy(t) for t in pre_trees]
    MCCs_slow = MTK.compute_mcc_pairs!(trees, OptArgs(consistent=true))
    @test get(MCCs_slow , (1, 2)) == get(MCCs_slow , (3, 2)) 
    joint_ = MTK.MCC_join_constraint([get(MCCs_slow , (1, 2)),  get(MCCs_slow , (1, 4))])
    @test get(MCCs_slow , (1, 4)) == get(MCCs_slow , (3, 4))
    @test get(MCCs_slow , (2, 4)) == joint_
    @test SplitList(trees[3]) == SplitList(trees[1])
    @test SplitList(trees[2])== split_list_final_t2
    @test SplitList(trees[4])== split_list_final_t4
end

pre_trees = [t1, t3, t2, t4]
split_list_final_t3 = [["D", "E"], ["D", "E", "F"], ["B", "C"], ["A", "B", "C"], ["A", "B", "C", "D", "E", "F"]]
split_list_final_t4 = [["E", "F"], ["D", "E", "F"], ["A", "B"], ["A", "B", "C"], ["A", "B", "C", "D", "E", "F"]]

@testset "4 trees (reordered), 1 round" begin
    trees = [copy(t) for t in pre_trees]
    MCCs = MTK.parallelized_compute_mccs!(trees, OptArgs(rounds=1, consistent=true))
    @test get(MCCs, (1, 3)) == get(MCCs, (3, 2)) 
    joint_ = MTK.MCC_join_constraint([get(MCCs, (1, 3)),  get(MCCs, (1, 4))])
    @test get(MCCs, (1, 4)) == get(MCCs, (2, 4))
    @test get(MCCs, (3, 4)) == joint_
    @test SplitList(trees[2]) == SplitList(trees[1])
    @test SplitList(trees[3])== split_list_final_t3
    @test SplitList(trees[4])== split_list_final_t4

    trees = [copy(t) for t in pre_trees]
    MCCs_slow = MTK.compute_mcc_pairs!(trees, OptArgs(rounds=1, consistent=true))
    @test get(MCCs_slow , (1, 3)) == get(MCCs_slow , (3, 2)) 
    joint_ = MTK.MCC_join_constraint([get(MCCs_slow , (1, 3)),  get(MCCs_slow , (1, 4))])
    @test get(MCCs_slow , (1, 4)) == get(MCCs_slow , (2, 4))
    @test get(MCCs_slow , (3, 4)) == joint_
    @test SplitList(trees[2]) == SplitList(trees[1])
    @test SplitList(trees[3])== split_list_final_t3
    @test SplitList(trees[4])== split_list_final_t4
end

@testset "4 trees (reordered), 2 rounds" begin
    trees = [copy(t) for t in pre_trees]
    MCCs = MTK.parallelized_compute_mccs!(trees, OptArgs(consistent=true))
    @test get(MCCs, (1, 3)) == get(MCCs, (3, 2)) 
    joint_ = MTK.MCC_join_constraint([get(MCCs, (1, 3)),  get(MCCs, (1, 4))])
    @test get(MCCs, (1, 4)) == get(MCCs, (2, 4))
    @test get(MCCs, (3, 4)) == joint_
    @test SplitList(trees[2]) == SplitList(trees[1])
    @test SplitList(trees[3])== split_list_final_t3
    @test SplitList(trees[4])== split_list_final_t4

    trees = [copy(t) for t in pre_trees]
    MCCs_slow = MTK.compute_mcc_pairs!(trees, OptArgs(consistent=true))
    @test get(MCCs_slow , (1, 3)) == get(MCCs_slow , (3, 2)) 
    joint_ = MTK.MCC_join_constraint([get(MCCs_slow , (1, 3)),  get(MCCs_slow , (1, 4))])
    @test get(MCCs_slow , (1, 4)) == get(MCCs_slow , (2, 4))
    @test get(MCCs_slow , (3, 4)) == joint_
    @test SplitList(trees[2]) == SplitList(trees[1])
    @test SplitList(trees[3])== split_list_final_t3
    @test SplitList(trees[4])== split_list_final_t4
end

pre_trees = [t1, t3, t2, t4, t5]

@testset "5 trees, 1 round" begin
    trees = [copy(t) for t in pre_trees]
    MCCs = MTK.parallelized_compute_mccs!(trees, OptArgs(consistent=true, rounds=1))
    joint_ = MTK.MCC_join_constraint([MTK.MCC_join_constraint([get(MCCs, (1, 3)),  get(MCCs, (1, 4))]), MTK.MCC_join_constraint([get(MCCs, (2, 3)),  get(MCCs, (2, 4))])])
    @test TreeKnit.is_MCC_subset(joint_, get(MCCs, (3,4)))
    joint_ = MTK.MCC_join_constraint([MTK.MCC_join_constraint([get(MCCs, (1, 5)),  get(MCCs, (1, 4))]), MTK.MCC_join_constraint([get(MCCs, (2, 5)),  get(MCCs, (2, 4))])])
    joint_ = MTK.MCC_join_constraint([MTK.MCC_join_constraint([get(MCCs, (3, 5)),  get(MCCs, (3, 4))]), joint_])
    @test TreeKnit.is_MCC_subset(joint_, get(MCCs, (5,4)))
end

pre_trees = [t1, t3, t2, t4, t5, t1, t2, t3]

@testset "8 trees, 1 round" begin
    trees = [copy(t) for t in pre_trees]
    MCCs = MTK.parallelized_compute_mccs!(trees, OptArgs(consistent=true, rounds=1))
    joint_ = MTK.MCC_join_constraint([MTK.MCC_join_constraint([get(MCCs, (1, 3)),  get(MCCs, (1, 4))]), MTK.MCC_join_constraint([get(MCCs, (2, 3)),  get(MCCs, (2, 4))])])
    @test TreeKnit.is_MCC_subset(joint_, get(MCCs, (3,4)))
    joint_ = MTK.MCC_join_constraint([MTK.MCC_join_constraint([get(MCCs, (1, 5)),  get(MCCs, (1, 4))]), MTK.MCC_join_constraint([get(MCCs, (2, 5)),  get(MCCs, (2, 4))])])
    joint_ = MTK.MCC_join_constraint([MTK.MCC_join_constraint([get(MCCs, (3, 5)),  get(MCCs, (3, 4))]), joint_])
    @test TreeKnit.is_MCC_subset(joint_, get(MCCs, (5,4)))
end

@testset "8 trees, 2 rounds" begin
    trees = [copy(t) for t in pre_trees]
    MCCs = MTK.parallelized_compute_mccs!(trees, OptArgs(consistent=true))
    @test typeof(MCCs) == MCC_set
end
