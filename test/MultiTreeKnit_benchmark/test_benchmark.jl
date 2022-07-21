using Test
using TreeTools
using TreeKnit

println("##### testing benchmark #####")

nwk1 = "((A,B),C);"
nwk2 = "(A,(B,C));"
nwk3 = "(A,B,C);"
t1 = node2tree(TreeTools.parse_newick(nwk1), label = "a")
t2 = node2tree(TreeTools.parse_newick(nwk2), label = "b")
t3 = node2tree(TreeTools.parse_newick(nwk3), label = "c")

solutions = [   [["A"], ["B", "C"]],
                [["B"], ["A", "C"]],
                [["C"], ["A", "B"]] 
            ]

@testset "benchmark with input trees" begin
    input_trees = [copy(t1), copy(t2), copy(t3)]
    MCC_dict, trees = TreeKnit.infer_benchmark_MCCs(input_trees, consistant=false)
    @test MCC_dict[Set(["a", "c"])] == [["A", "B", "C"]]
    @test in(MCC_dict[Set(["a", "b"])], solutions)
	@test SplitList(t1) == SplitList(trees[3])

    input_trees = [copy(t2), copy(t1), copy(t3)]
    MCC_dict, trees = TreeKnit.infer_benchmark_MCCs(input_trees, tree_names, consistant=false)
    @test MCC_dict[Set(["b", "c"])] == [["A", "B", "C"]]
    @test in(MCC_dict[Set(["a", "b"])], solutions)
	@test SplitList(t2) == SplitList(trees[3])

    input_trees = [copy(t3), copy(t1), copy(t2)]
    MCC_dict, trees = TreeKnit.infer_benchmark_MCCs(input_trees, tree_names, order="input", consistant=false)
    @test MCC_dict[Set(["a", "c"])] == [["A", "B", "C"]]
    @test in(MCC_dict[Set(["a", "b"])], solutions)
	@test SplitList(t1) == SplitList(trees[1])

end

@testset "consistent benchmark" begin
    input_trees = [copy(t1), copy(t2), copy(t3)]
    MCC_dict, trees = TreeKnit.infer_benchmark_MCCs(input_trees, tree_names, consistant=true)
    @test MCC_dict[Set(["a", "c"])] == [["A", "B", "C"]]
    @test MCC_dict[Set(["a", "b"])] == MCC_dict[["b", "c"]]
    @test SplitList(t1) == SplitList(trees[3])
    @test 0.0 == TreeKnit.consistency_rate(MCC_dict, trees)

    input_trees = [copy(t2), copy(t1), copy(t3)]
    MCC_dict, trees = TreeKnit.infer_benchmark_MCCs(input_trees, tree_names, consistant=true)
    @test 0.0 == TreeKnit.consistency_rate(MCC_dict, trees)

    input_trees = [copy(t3), copy(t1), copy(t2)]
    MCC_dict, trees = TreeKnit.infer_benchmark_MCCs(input_trees, tree_names, order="input", consistant=true)
    @test 0.0 == TreeKnit.consistency_rate(MCC_dict, trees)
end

@testset "benchmark input w/o name" begin
    input_trees = [copy(t1), copy(t2), copy(t3)]
    MCC_dict, trees = TreeKnit.infer_benchmark_MCCs(input_trees, consistant=false)
    @test MCC_dict[Set(["a", "b"])] == [["A", "B", "C"]]
    @test in(MCC_dict[Set(["a", "c"])], solutions)
	@test SplitList(t1) == SplitList(trees[1])
end

@testset "benchmark simulated trees" begin
    MCC_dict, trees = TreeKnit.infer_benchmark_MCCs(3, 10, debug=true, order="resolution")
    @test length(keys(MCC_dict)) == 4
    MCC_dict, trees = TreeKnit.infer_benchmark_MCCs(3, 10, debug=true, remove=true, order="resolution")
    @test length(keys(MCC_dict)) == 4
end

