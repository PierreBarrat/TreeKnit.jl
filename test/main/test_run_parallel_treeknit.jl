using Test
using Dagger
using TreeTools
using TreeKnit

println("#### test parallelized recursive TreeKnit ###")

nwk1 = "((A,B),C)R;"
nwk2 = "(A,(B,C))R;"
nwk3 = "(A,B,C)R;"
t1 = parse_newick_string(nwk1, node_data_type=TreeTools.MiscData, label = "a")
t2 = parse_newick_string(nwk2, node_data_type=TreeTools.MiscData, label = "b")
t3 = parse_newick_string(nwk3, node_data_type=TreeTools.MiscData, label = "c")
pre_trees = [t1, t2, t3]

@testset "3 trees" begin
    trees = [copy(t) for t in pre_trees]
    resolve!(trees...)
    MCCs = TreeKnit.run_parallel_treeknit!(trees, OptArgs(3; method=:better_MCCs))
    @test SplitList(trees[3]) == SplitList(trees[1])

    trees = [copy(t) for t in pre_trees]
    resolve!(trees...)
    MCCs = TreeKnit.run_parallel_treeknit!(trees, OptArgs(3; method=:better_trees))
    @test SplitList(trees[3]) != SplitList(trees[1])
    
    trees = [copy(t) for t in pre_trees]
    resolve!(trees...)
    MCCs_slow = TreeKnit.run_standard_treeknit!(trees, OptArgs(3; method=:better_MCCs))
    @test SplitList(trees[3]) == SplitList(trees[1])

    trees = [copy(t) for t in pre_trees]
    resolve!(trees...)
    MCCs_slow = TreeKnit.run_standard_treeknit!(trees, OptArgs(3; rounds=1, method=:better_trees))
    @test SplitList(trees[3]) != SplitList(trees[1])
end

nwk1 = "(((A,B),C),(D,E,F))R;"
nwk2 = "((A,(B,C)),(D,E,F))R;"
nwk3 = "((A,B,C),((D,E),F))R;"
nwk4 = "((A,B,C),(D,(E,F)))R;"
nwk5 = "((A,(B,C)),(D,(E,F)))R;"
t1 = parse_newick_string(nwk1, label = "a")
t2 = parse_newick_string(nwk2, label = "b")
t3 = parse_newick_string(nwk3, label = "c")
t4 = parse_newick_string(nwk4, label = "d")
t5 = parse_newick_string(nwk5, label = "e")
pre_trees = [t1, t2, t3, t4]

pre_trees_for_preresolve = [t2, t3, t4, t5]

split_list_final_t2 = [["D", "E"], ["D", "E", "F"], ["B", "C"], ["A", "B", "C"], ["A", "B", "C", "D", "E", "F"]]
split_list_final_t4 = [["E", "F"], ["D", "E", "F"], ["A", "B"], ["A", "B", "C"], ["A", "B", "C", "D", "E", "F"]]


@testset "4 trees" begin
    trees = [copy(t) for t in pre_trees]
    resolve!(trees...)
    MCCs = TreeKnit.run_parallel_treeknit!(trees, OptArgs(4; method=:better_MCCs))
    @test SplitList(trees[3]) == SplitList(trees[1])
    @test SplitList(trees[2])== split_list_final_t2
    @test SplitList(trees[4])== split_list_final_t4

    trees = [copy(t) for t in pre_trees]
    resolve!(trees...)
    MCCs_slow = TreeKnit.run_standard_treeknit!(trees, OptArgs(4; method=:better_MCCs))
    @test SplitList(trees[3]) == SplitList(trees[1])
    @test SplitList(trees[2])== split_list_final_t2
    @test SplitList(trees[4])== split_list_final_t4

    trees = [copy(t) for t in pre_trees]
    resolve!(trees...)
    MCCs = TreeKnit.run_parallel_treeknit!(trees, OptArgs(4; method=:better_trees))
    @test all([SplitList(trees[i]) == SplitList(pre_trees[i]) for i in 1:4])

    trees = [copy(t) for t in pre_trees_for_preresolve]
    resolve!(trees...)
    MCCs = TreeKnit.run_parallel_treeknit!(trees, OptArgs(4; method=:better_trees))
    @test all([Split([2,3]) in SplitList(trees[i]) for i in 1:4])
end

pre_trees = [t1, t3, t2, t4]
split_list_final_t3 = [["D", "E"], ["D", "E", "F"], ["B", "C"], ["A", "B", "C"], ["A", "B", "C", "D", "E", "F"]]
split_list_final_t4 = [["E", "F"], ["D", "E", "F"], ["A", "B"], ["A", "B", "C"], ["A", "B", "C", "D", "E", "F"]]

@testset "4 trees (reordered)" begin
    trees = [copy(t) for t in pre_trees]
    resolve!(trees...)
    MCCs = TreeKnit.run_parallel_treeknit!(trees, OptArgs(4; method=:better_MCCs))
    @test SplitList(trees[2]) == SplitList(trees[1])
    @test SplitList(trees[3])== split_list_final_t3
    @test SplitList(trees[4])== split_list_final_t4

    trees = [copy(t) for t in pre_trees]
    resolve!(trees...)
    MCCs_slow = TreeKnit.run_standard_treeknit!(trees, OptArgs(4; method=:better_MCCs))
    @test SplitList(trees[2]) == SplitList(trees[1])
    @test SplitList(trees[3])== split_list_final_t3
    @test SplitList(trees[4])== split_list_final_t4
end

pre_trees = [t1, t3, t2, t4, t5, t1, t2, t3]

@testset "8 trees" begin
    trees = [copy(t) for t in pre_trees]
    resolve!(trees...)
    MCCs = TreeKnit.run_parallel_treeknit!(trees, OptArgs(8; method=:better_trees))
    @test all([SplitList(trees[i]) == SplitList(pre_trees[i]) for i in 1:4])
end
