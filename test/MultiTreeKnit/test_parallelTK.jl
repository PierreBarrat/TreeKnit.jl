using Dagger, BenchmarkTools
using Plots, Test
using TreeKnit, TreeTools

function run_step!(tree1::Tree, tree2::Tree, constraints)
#  function run_step!(tree1::Tree, tree2::Tree, constraints::Union{Nothing, Vector{Vector{Vector{String}}}})
    # if !isnothing(constraints)
    #     constraint = constraints[1]
    #     for mcc_con in constraints[2:end]
    #         constraint = TreeKnit.join_sets([constraint, mcc_con])
    #     end
    # else
    if !isnothing(constraints)
        constraint = fetch(constraints[1])
        for mcc_con in constraints[2:end]
            constraint = TreeKnit.join_sets([constraint, fetch(mcc_con)])
        end
    else
        constraint = nothing
    end
    MCC = TreeKnit.runopt(TreeKnit.OptArgs(;constraint=constraint), tree1, tree2; output = :mccs)
    rS = TreeKnit.resolve!(tree1, tree2, MCC)
    TreeTools.ladderize!(tree1)
    TreeKnit.sort_polytomies!(tree1, tree2, MCC)
    return MCC
end

function parallelized_compute_mccs(trees; consistant=true, rounds=2)
    l_t = length(trees)
    #pair_MCCs = MCC_set(l_t, [t.label for t in trees])
    parallel_MCCs = Dict()
    for r in 1:rounds
        for i in 1:(l_t-1)
            for j in (i+1):l_t
                MCC_list = nothing
                if consistant && (i>1 || r>1)
                    MCC_list = []
                    #MCC_list = Vector{Vector{String}}[]
                    if r>1
                        range = filter(e->e∉Set([i,j]), 1:l_t)
                    else
                        range = 1:(i-1)
                    end
                    for x in range
                        append!(MCC_list, [parallel_MCCs[Set([i, x])], parallel_MCCs[Set([j,x])]])
                        #Dagger.@spawn append!(MCC_list, [pair_MCCs[Set([trees[i].label, trees[x].label])]])
                        #Dagger.@spawn append!(MCC_list, [pair_MCCs[Set([trees[j].label,trees[x].label])]])
                    end
                end
                parallel_MCCs[Set([i,j])] = Dagger.@spawn run_step!(trees[i], trees[j], MCC_list)
                #pair_MCCs.mccs[Set([trees[i].label, trees[j].label])] = fetch(MCC)
            end
        end
    end
    pair_MCCs = MCC_set(l_t, [t.label for t in trees])
    for (key, mcc) in parallel_MCCs
        TreeKnit.add!(pair_MCCs, fetch(mcc), Tuple(key))
    end
    return pair_MCCs
end

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
    MCCs = parallelized_compute_mccs(trees, rounds=1)
    @test SplitList(trees[3]) == SplitList(trees[1])
    @test get(MCCs, (1, 2)) == get(MCCs, (3, 2)) 
    
    trees = [copy(t) for t in pre_trees]
    MCCs_slow = TreeKnit.compute_mcc_pairs(trees, rounds=1)
    @test SplitList(trees[3]) == SplitList(trees[1])
    @test get(MCCs_slow, (1, 2)) == get(MCCs_slow, (3, 2)) 
end

print("Time parallelized, 1 round")
trees = [copy(t) for t in pre_trees]
@btime MCCs = parallelized_compute_mccs(trees, rounds=1)
print("Time non-parallelized, 1 round")
trees = [copy(t) for t in pre_trees]
@btime TreeKnit.compute_mcc_pairs(trees, rounds=1)
print("Time parallelized, 1 round")
trees = [copy(t) for t in pre_trees]
@btime MCCs = parallelized_compute_mccs(trees, rounds=1)
print("Time non-parallelized, 1 round")
trees = [copy(t) for t in pre_trees]
@btime TreeKnit.compute_mcc_pairs(trees, rounds=1)

@testset "3 trees, 2 rounds" begin
    trees = [copy(t) for t in pre_trees]
    MCCs = parallelized_compute_mccs(trees, rounds=2)
    @test SplitList(trees[3]) == SplitList(trees[1])
    @test get(MCCs, (1, 2)) == get(MCCs, (3, 2)) 

    trees = [copy(t) for t in pre_trees]
    MCCs_slow = TreeKnit.compute_mcc_pairs(trees, rounds=2)
    @test SplitList(trees[3]) == SplitList(trees[1])
    @test get(MCCs_slow, (1, 2)) == get(MCCs_slow, (3, 2)) 
end

print("Time parallelized, 2 rounds")
trees = [copy(t) for t in pre_trees]
@btime MCCs = parallelized_compute_mccs(trees, rounds=2)
print("Time non-parallelized, 2 rounds")
trees = [copy(t) for t in pre_trees]
@btime  MCCs_slow = TreeKnit.compute_mcc_pairs(trees, rounds=2)
print("Time parallelized, 2 rounds")
trees = [copy(t) for t in pre_trees]
@btime MCCs = parallelized_compute_mccs(trees, rounds=2)
print("Time non-parallelized, 2 rounds")
trees = [copy(t) for t in pre_trees]
@btime  MCCs_slow = TreeKnit.compute_mcc_pairs(trees, rounds=2)

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
    MCCs = parallelized_compute_mccs(trees, rounds=1)
    @test get(MCCs, (1, 2)) == get(MCCs, (3, 2)) 
    joint_ = TreeKnit.join_sets([get(MCCs, (1, 2)),  get(MCCs, (1, 4))])
    @test get(MCCs, (1, 4)) == get(MCCs, (3, 4))
    @test get(MCCs, (2, 4)) == joint_
    @test SplitList(trees[3]) == SplitList(trees[1])
    @test SplitList(trees[2])== split_list_final_t2
    @test SplitList(trees[4])== split_list_final_t4

    trees = [copy(t) for t in pre_trees]
    MCCs_slow = TreeKnit.compute_mcc_pairs(trees, rounds=1)
    @test get(MCCs_slow , (1, 2)) == get(MCCs_slow , (3, 2)) 
    joint_ = TreeKnit.join_sets([get(MCCs_slow , (1, 2)),  get(MCCs_slow , (1, 4))])
    @test get(MCCs_slow , (1, 4)) == get(MCCs_slow , (3, 4))
    @test get(MCCs_slow , (2, 4)) == joint_
    @test SplitList(trees[3]) == SplitList(trees[1])
    @test SplitList(trees[2])== split_list_final_t2
    @test SplitList(trees[4])== split_list_final_t4
end

print("Time parallelized, 1 round")
trees = [copy(t) for t in pre_trees]
@btime MCCs = parallelized_compute_mccs(trees, rounds=1)
print("Time non-parallelized, 1 round")
trees = [copy(t) for t in pre_trees]
@btime MCCs_slow = TreeKnit.compute_mcc_pairs(trees, rounds=1)
print("Time parallelized, 1 round")
trees = [copy(t) for t in pre_trees]
@btime MCCs = parallelized_compute_mccs(trees, rounds=1)
print("Time non-parallelized, 1 round")
trees = [copy(t) for t in pre_trees]
@btime MCCs_slow = TreeKnit.compute_mcc_pairs(trees, rounds=1)

@testset "4 trees, 2 rounds" begin
    trees = [copy(t) for t in pre_trees]
    MCCs = parallelized_compute_mccs(trees, rounds=2)
    @test get(MCCs, (1, 2)) == get(MCCs, (3, 2)) 
    joint_ = TreeKnit.join_sets([get(MCCs, (1, 2)),  get(MCCs, (1, 4))])
    @test get(MCCs, (1, 4)) == get(MCCs, (3, 4))
    @test get(MCCs, (2, 4)) == joint_
    @test SplitList(trees[3]) == SplitList(trees[1])
    @test SplitList(trees[2])== split_list_final_t2
    @test SplitList(trees[4])== split_list_final_t4

    trees = [copy(t) for t in pre_trees]
    MCCs_slow = TreeKnit.compute_mcc_pairs(trees, rounds=2)
    @test get(MCCs_slow , (1, 2)) == get(MCCs_slow , (3, 2)) 
    joint_ = TreeKnit.join_sets([get(MCCs_slow , (1, 2)),  get(MCCs_slow , (1, 4))])
    @test get(MCCs_slow , (1, 4)) == get(MCCs_slow , (3, 4))
    @test get(MCCs_slow , (2, 4)) == joint_
    @test SplitList(trees[3]) == SplitList(trees[1])
    @test SplitList(trees[2])== split_list_final_t2
    @test SplitList(trees[4])== split_list_final_t4
end
print("Time parallelized, 1 round")
trees = [copy(t) for t in pre_trees]
@btime MCCs = parallelized_compute_mccs(trees, rounds=2)
print("Time non-parallelized, 1 round")
trees = [copy(t) for t in pre_trees]
@btime MCCs_slow = TreeKnit.compute_mcc_pairs(trees, rounds=2)
print("Time parallelized, 1 round")
trees = [copy(t) for t in pre_trees]
@btime MCCs = parallelized_compute_mccs(trees, rounds=2)
print("Time non-parallelized, 1 round")
trees = [copy(t) for t in pre_trees]
@btime MCCs_slow = TreeKnit.compute_mcc_pairs(trees, rounds=2)

pre_trees = [t1, t3, t2, t4]
split_list_final_t3 = [["D", "E"], ["D", "E", "F"], ["B", "C"], ["A", "B", "C"], ["A", "B", "C", "D", "E", "F"]]
split_list_final_t4 = [["E", "F"], ["D", "E", "F"], ["A", "B"], ["A", "B", "C"], ["A", "B", "C", "D", "E", "F"]]

@testset "4 trees (reordered), 1 round" begin
    trees = [copy(t) for t in pre_trees]
    MCCs = parallelized_compute_mccs(trees, rounds=1)
    @test get(MCCs, (1, 3)) == get(MCCs, (3, 2)) 
    joint_ = TreeKnit.join_sets([get(MCCs, (1, 3)),  get(MCCs, (1, 4))])
    @test get(MCCs, (1, 4)) == get(MCCs, (2, 4))
    @test get(MCCs, (3, 4)) == joint_
    @test SplitList(trees[2]) == SplitList(trees[1])
    @test SplitList(trees[3])== split_list_final_t3
    @test SplitList(trees[4])== split_list_final_t4

    trees = [copy(t) for t in pre_trees]
    MCCs_slow = TreeKnit.compute_mcc_pairs(trees, rounds=1)
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
    MCCs = parallelized_compute_mccs(trees, rounds=2)
    @test get(MCCs, (1, 3)) == get(MCCs, (3, 2)) 
    joint_ = TreeKnit.join_sets([get(MCCs, (1, 3)),  get(MCCs, (1, 4))])
    @test get(MCCs, (1, 4)) == get(MCCs, (2, 4))
    @test get(MCCs, (3, 4)) == joint_
    @test SplitList(trees[2]) == SplitList(trees[1])
    @test SplitList(trees[3])== split_list_final_t3
    @test SplitList(trees[4])== split_list_final_t4

    trees = [copy(t) for t in pre_trees]
    MCCs_slow = TreeKnit.compute_mcc_pairs(trees, rounds=2)
    @test get(MCCs_slow , (1, 3)) == get(MCCs_slow , (3, 2)) 
    joint_ = TreeKnit.join_sets([get(MCCs_slow , (1, 3)),  get(MCCs_slow , (1, 4))])
    @test get(MCCs_slow , (1, 4)) == get(MCCs_slow , (2, 4))
    @test get(MCCs_slow , (3, 4)) == joint_
    @test SplitList(trees[2]) == SplitList(trees[1])
    @test SplitList(trees[3])== split_list_final_t3
    @test SplitList(trees[4])== split_list_final_t4
end

