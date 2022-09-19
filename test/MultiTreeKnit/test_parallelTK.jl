using Dagger, BenchmarkTools, StatsBase
using Plots, Test
using TreeKnit, TreeTools
print_time = false

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
    MCCs = TreeKnit.parallelized_compute_mccs!(trees, OptArgs(rounds=1))
    @test SplitList(trees[3]) == SplitList(trees[1])
    @test get(MCCs, (1, 2)) == get(MCCs, (3, 2)) 
    
    trees = [copy(t) for t in pre_trees]
    MCCs_slow = TreeKnit.compute_mcc_pairs!(trees, OptArgs(rounds=1))
    @test SplitList(trees[3]) == SplitList(trees[1])
    @test get(MCCs_slow, (1, 2)) == get(MCCs_slow, (3, 2)) 
end

if print_time == true
    print("Time parallelized, 1 round")
    trees = [copy(t) for t in pre_trees]
    @btime MCCs = TreeKnit.parallelized_compute_mccs!(trees, OptArgs(rounds=1))
    print("Time non-parallelized, 1 round")
    trees = [copy(t) for t in pre_trees]
    @btime TreeKnit.compute_mcc_pairs!(trees, OptArgs(rounds=1))
    print("Time parallelized, 1 round")
    trees = [copy(t) for t in pre_trees]
    @btime MCCs = TreeKnit.parallelized_compute_mccs!(trees, OptArgs(rounds=1))
    print("Time non-parallelized, 1 round")
    trees = [copy(t) for t in pre_trees]
    @btime TreeKnit.compute_mcc_pairs!(trees, OptArgs(rounds=1))
end

@testset "3 trees, 2 rounds" begin
    trees = [copy(t) for t in pre_trees]
    MCCs = TreeKnit.parallelized_compute_mccs!(trees, OptArgs())
    @test SplitList(trees[3]) == SplitList(trees[1])
    @test get(MCCs, (1, 2)) == get(MCCs, (3, 2)) 

    trees = [copy(t) for t in pre_trees]
    MCCs_slow = TreeKnit.compute_mcc_pairs!(trees, OptArgs())
    @test SplitList(trees[3]) == SplitList(trees[1])
    @test get(MCCs_slow, (1, 2)) == get(MCCs_slow, (3, 2)) 
end

if print_time == true
    print("Time parallelized, 2 rounds")
    trees = [copy(t) for t in pre_trees]
    @btime MCCs = TreeKnit.parallelized_compute_mccs!(trees, OptArgs())
    print("Time non-parallelized, 2 rounds")
    trees = [copy(t) for t in pre_trees]
    @btime  MCCs_slow = TreeKnit.compute_mcc_pairs!(trees, OptArgs())
    print("Time parallelized, 2 rounds")
    trees = [copy(t) for t in pre_trees]
    @btime MCCs = TreeKnit.parallelized_compute_mccs!(trees, OptArgs())
    print("Time non-parallelized, 2 rounds")
    trees = [copy(t) for t in pre_trees]
    @btime  MCCs_slow = TreeKnit.compute_mcc_pairs!(trees, OptArgs())
end

nwk1 = "(((A,B),C),(D,E,F))R"
nwk2 = "((A,(B,C)),(D,E,F))R"
nwk3 = "((A,B,C),((D,E),F))R"
nwk4 = "((A,B,C),(D,(E,F)))R"
t1 = node2tree(TreeTools.parse_newick(nwk1, node_data_type=TreeTools.MiscData), label = "a")
t2 = node2tree(TreeTools.parse_newick(nwk2, node_data_type=TreeTools.MiscData), label = "b")
t3 = node2tree(TreeTools.parse_newick(nwk3, node_data_type=TreeTools.MiscData), label = "c")
t4 = node2tree(TreeTools.parse_newick(nwk4, node_data_type=TreeTools.MiscData), label = "d")
pre_trees = [t1, t2, t3, t4]

split_list_final_t2 = [["D", "E"], ["D", "E", "F"], ["B", "C"], ["A", "B", "C"], ["A", "B", "C", "D", "E", "F"]]
split_list_final_t4 = [["E", "F"], ["D", "E", "F"], ["A", "B"], ["A", "B", "C"], ["A", "B", "C", "D", "E", "F"]]

@testset "4 trees, 1 round" begin
    trees = [copy(t) for t in pre_trees]
    MCCs = TreeKnit.parallelized_compute_mccs!(trees, OptArgs(rounds=1))
    @test get(MCCs, (1, 2)) == get(MCCs, (3, 2)) 
    joint_ = TreeKnit.join_sets([get(MCCs, (1, 2)),  get(MCCs, (1, 4))])
    @test get(MCCs, (1, 4)) == get(MCCs, (3, 4))
    @test get(MCCs, (2, 4)) == joint_
    @test SplitList(trees[3]) == SplitList(trees[1])
    @test SplitList(trees[2])== split_list_final_t2
    @test SplitList(trees[4])== split_list_final_t4

    trees = [copy(t) for t in pre_trees]
    MCCs_slow = TreeKnit.compute_mcc_pairs!(trees, OptArgs(rounds=1))
    @test get(MCCs_slow , (1, 2)) == get(MCCs_slow , (3, 2)) 
    joint_ = TreeKnit.join_sets([get(MCCs_slow , (1, 2)),  get(MCCs_slow , (1, 4))])
    @test get(MCCs_slow , (1, 4)) == get(MCCs_slow , (3, 4))
    @test get(MCCs_slow , (2, 4)) == joint_
    @test SplitList(trees[3]) == SplitList(trees[1])
    @test SplitList(trees[2])== split_list_final_t2
    @test SplitList(trees[4])== split_list_final_t4
end

if print_time == true
    print("Time parallelized, 1 round")
    trees = [copy(t) for t in pre_trees]
    @btime MCCs = TreeKnit.parallelized_compute_mccs!(trees, OptArgs(rounds=1))
    print("Time non-parallelized, 1 round")
    trees = [copy(t) for t in pre_trees]
    @btime MCCs_slow = TreeKnit.compute_mcc_pairs!(trees, OptArgs(rounds=1))
    print("Time parallelized, 1 round")
    trees = [copy(t) for t in pre_trees]
    @btime MCCs = TreeKnit.parallelized_compute_mccs!(trees, OptArgs(rounds=1))
    print("Time non-parallelized, 1 round")
    trees = [copy(t) for t in pre_trees]
    @btime MCCs_slow = TreeKnit.compute_mcc_pairs!(trees, OptArgs(rounds=1))
end

@testset "4 trees, 2 rounds" begin
    trees = [copy(t) for t in pre_trees]
    MCCs = TreeKnit.parallelized_compute_mccs!(trees, OptArgs())
    @test get(MCCs, (1, 2)) == get(MCCs, (3, 2)) 
    joint_ = TreeKnit.join_sets([get(MCCs, (1, 2)),  get(MCCs, (1, 4))])
    @test get(MCCs, (1, 4)) == get(MCCs, (3, 4))
    @test get(MCCs, (2, 4)) == joint_
    @test SplitList(trees[3]) == SplitList(trees[1])
    @test SplitList(trees[2])== split_list_final_t2
    @test SplitList(trees[4])== split_list_final_t4

    trees = [copy(t) for t in pre_trees]
    MCCs_slow = TreeKnit.compute_mcc_pairs!(trees, OptArgs())
    @test get(MCCs_slow , (1, 2)) == get(MCCs_slow , (3, 2)) 
    joint_ = TreeKnit.join_sets([get(MCCs_slow , (1, 2)),  get(MCCs_slow , (1, 4))])
    @test get(MCCs_slow , (1, 4)) == get(MCCs_slow , (3, 4))
    @test get(MCCs_slow , (2, 4)) == joint_
    @test SplitList(trees[3]) == SplitList(trees[1])
    @test SplitList(trees[2])== split_list_final_t2
    @test SplitList(trees[4])== split_list_final_t4
end

if print_time == true
    print("Time parallelized, 1 round")
    trees = [copy(t) for t in pre_trees]
    @btime MCCs = TreeKnit.parallelized_compute_mccs!(trees, OptArgs())
    print("Time non-parallelized, 1 round")
    trees = [copy(t) for t in pre_trees]
    @btime MCCs_slow = TreeKnit.compute_mcc_pairs!(trees, OptArgs())
    print("Time parallelized, 1 round")
    trees = [copy(t) for t in pre_trees]
    @btime MCCs = TreeKnit.parallelized_compute_mccs!(trees, OptArgs())
    print("Time non-parallelized, 1 round")
    trees = [copy(t) for t in pre_trees]
    @btime MCCs_slow = TreeKnit.compute_mcc_pairs!(trees, OptArgs())
end

pre_trees = [t1, t3, t2, t4]
split_list_final_t3 = [["D", "E"], ["D", "E", "F"], ["B", "C"], ["A", "B", "C"], ["A", "B", "C", "D", "E", "F"]]
split_list_final_t4 = [["E", "F"], ["D", "E", "F"], ["A", "B"], ["A", "B", "C"], ["A", "B", "C", "D", "E", "F"]]

@testset "4 trees (reordered), 1 round" begin
    trees = [copy(t) for t in pre_trees]
    MCCs = TreeKnit.parallelized_compute_mccs!(trees, OptArgs(rounds=1))
    @test get(MCCs, (1, 3)) == get(MCCs, (3, 2)) 
    joint_ = TreeKnit.join_sets([get(MCCs, (1, 3)),  get(MCCs, (1, 4))])
    @test get(MCCs, (1, 4)) == get(MCCs, (2, 4))
    @test get(MCCs, (3, 4)) == joint_
    @test SplitList(trees[2]) == SplitList(trees[1])
    @test SplitList(trees[3])== split_list_final_t3
    @test SplitList(trees[4])== split_list_final_t4

    trees = [copy(t) for t in pre_trees]
    MCCs_slow = TreeKnit.compute_mcc_pairs!(trees, OptArgs(rounds=1))
    @test get(MCCs_slow , (1, 3)) == get(MCCs_slow , (3, 2)) 
    joint_ = TreeKnit.join_sets([get(MCCs_slow , (1, 3)),  get(MCCs_slow , (1, 4))])
    @test get(MCCs_slow , (1, 4)) == get(MCCs_slow , (2, 4))
    @test get(MCCs_slow , (3, 4)) == joint_
    @test SplitList(trees[2]) == SplitList(trees[1])
    @test SplitList(trees[3])== split_list_final_t3
    @test SplitList(trees[4])== split_list_final_t4
end

@testset "4 trees (reordered), 2 rounds" begin
    trees = [copy(t) for t in pre_trees]
    MCCs = TreeKnit.parallelized_compute_mccs!(trees, OptArgs())
    @test get(MCCs, (1, 3)) == get(MCCs, (3, 2)) 
    joint_ = TreeKnit.join_sets([get(MCCs, (1, 3)),  get(MCCs, (1, 4))])
    @test get(MCCs, (1, 4)) == get(MCCs, (2, 4))
    @test get(MCCs, (3, 4)) == joint_
    @test SplitList(trees[2]) == SplitList(trees[1])
    @test SplitList(trees[3])== split_list_final_t3
    @test SplitList(trees[4])== split_list_final_t4

    trees = [copy(t) for t in pre_trees]
    MCCs_slow = TreeKnit.compute_mcc_pairs!(trees, OptArgs())
    @test get(MCCs_slow , (1, 3)) == get(MCCs_slow , (3, 2)) 
    joint_ = TreeKnit.join_sets([get(MCCs_slow , (1, 3)),  get(MCCs_slow , (1, 4))])
    @test get(MCCs_slow , (1, 4)) == get(MCCs_slow , (2, 4))
    @test get(MCCs_slow , (3, 4)) == joint_
    @test SplitList(trees[2]) == SplitList(trees[1])
    @test SplitList(trees[3])== split_list_final_t3
    @test SplitList(trees[4])== split_list_final_t4
end

for i in 1:100
    true_trees, arg = TreeKnit.get_trees(8, 15; ρ=(10^-0.1))
    labels = [t.label for t in true_trees]
    rMCCs = TreeKnit.convert_MCC_list_to_set(8, labels, TreeKnit.get_real_MCCs(8, arg))
    for no_trees in 2:8
        rand_order = sample(1:8, no_trees, replace = false)
        unresolved_trees = [copy(t) for t in true_trees[rand_order]]
        unresolved_trees = TreeKnit.remove_branches(unresolved_trees; c=0.75)

        i_trees = [copy(t) for t in unresolved_trees]
        try
            fc_i_MCCs = TreeKnit.get_infered_MCC_pairs!(i_trees, TreeKnit.OptArgs(consistent = true, force_consist=true, force_topo_consist=false, parallel=true))
        catch e
            print("n: "*string(no_trees))
            for t in true_trees
                print_tree_ascii(" ", t)
            end
            throw(error())
        end
    end
end

