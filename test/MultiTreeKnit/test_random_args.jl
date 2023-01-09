using ARGTools
using MTKTools
using Random
using StatsBase

println("######## MultiTreeKnit on random ARGs ##########")

@testset "check is_degenerate(MCCs) == (consistency_rate!=0.0) holds for random multitrees" begin
    repeat = 0
    while repeat < 10
        trees, arg = MTKTools.get_trees(3, 50, remove=true, c=0.75);
        label!(trees[1], "a")
        label!(trees[2], "b")
        label!(trees[3], "c")
        iMCCs = MTK.get_infered_MCC_pairs!(trees)
        @test MTKTools.is_degenerate(iMCCs) == (MTKTools.consistency_rate(iMCCs, trees)!=0.0)
        repeat += 1
    end
end

function check_parallelized_TK()
    for i in 1:100
        true_trees, arg = MTKTools.get_trees(8, 15; ρ=(10^-0.1))
        labels = [t.label for t in true_trees]
        rMCCs = MCC_set(8, labels, MTKTools.get_real_MCCs(8, arg))
        for no_trees in 2:8
            rand_order = sample(1:8, no_trees, replace = false)
            unresolved_trees = [copy(t) for t in true_trees[rand_order]]
            unresolved_trees = MTKTools.remove_branches(unresolved_trees; c=0.75)

            i_trees = [copy(t) for t in unresolved_trees]
            try
                fc_i_MCCs = MTK.get_infered_MCC_pairs!(i_trees, TreeKnit.OptArgs(parallel=true))
            catch e
                print("n: "*string(no_trees))
                for t in true_trees
                    print_tree_ascii(" ", t)
                end
                throw(error())
            end
        end
    end
    return true
end

@testset "check parallelized TK runs correctly for random multitrees" begin
    @test check_parallelized_TK()
end