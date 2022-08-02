using Test
using TreeTools
using TreeKnit

println("##### testing benchmark #####")

nwk1 = "((A,B),C);"
nwk2 = "(A,(B,C));"
nwk3 = "(A,B,C);"
t1 = node2tree(TreeTools.parse_newick(nwk1, node_data_type=TreeTools.MiscData), label = "a")
t2 = node2tree(TreeTools.parse_newick(nwk2, node_data_type=TreeTools.MiscData), label = "b")
t3 = node2tree(TreeTools.parse_newick(nwk3, node_data_type=TreeTools.MiscData), label = "c")

solutions = [   [["A"], ["B", "C"]],
                [["B"], ["A", "C"]],
                [["C"], ["A", "B"]] 
            ]

@testset "benchmark with input trees" begin
    input_trees = [copy(t1), copy(t2), copy(t3)]
    MCC_dict = TreeKnit.infer_benchmark_MCCs!(input_trees, consistant=false)
    @test MCC_dict[Set(["a", "c"])] == [["A", "B", "C"]]
    @test in(MCC_dict[Set(["a", "b"])], solutions)
	@test SplitList(t1) == SplitList(input_trees[3])

    input_trees = [copy(t2), copy(t1), copy(t3)]
    MCC_dict = TreeKnit.infer_benchmark_MCCs!(input_trees, consistant=false)
    @test MCC_dict[Set(["b", "c"])] == [["A", "B", "C"]]
    @test in(MCC_dict[Set(["a", "b"])], solutions)
	@test SplitList(t2) == SplitList(input_trees[3])

    input_trees = [copy(t3), copy(t1), copy(t2)]
    MCC_dict = TreeKnit.infer_benchmark_MCCs!(input_trees, order="input", consistant=false)
    @test MCC_dict[Set(["a", "c"])] == [["A", "B", "C"]]
    @test in(MCC_dict[Set(["a", "b"])], solutions)
	@test SplitList(t1) == SplitList(input_trees[1])

end

@testset "consistency shared_branch_constraint" begin
    input_trees = [copy(t1), copy(t2), copy(t3)]
    ot3 = copy(convert(Tree{TreeTools.MiscData}, input_trees[1]))
    ot2 = copy(convert(Tree{TreeTools.MiscData}, input_trees[2]))
    TreeKnit.mark_shared_branches!([["A"], ["B", "C"]], ot2, ot3)
    @test ot2.lnodes["NODE_1"].data.dat["shared_branch_constraint"] == false
    @test ot3.lnodes["NODE_1"].data.dat["shared_branch_constraint"] == true
    @test ot2.lleaves["A"].data.dat["shared_branch_constraint"] == false
    @test ot3.lleaves["A"].data.dat["shared_branch_constraint"] == false
    @test ot2.lleaves["B"].data.dat["shared_branch_constraint"] == true
    @test ot3.lleaves["B"].data.dat["shared_branch_constraint"] == true

    input_trees = [copy(t1), copy(t2), copy(t3)]
    for constraint in solutions
        mccs = TreeKnit.runopt(TreeKnit.OptArgs(;constraint=constraint), input_trees[1], input_trees[2]; output = :mccs)
        @test mccs == constraint
    end

end

@testset "consistent benchmark" begin
    input_trees = [copy(t1), copy(t2), copy(t3)]
    MCC_dict = TreeKnit.infer_benchmark_MCCs!(input_trees, consistant=true)
    @test MCC_dict[Set(["a", "c"])] == [["A", "B", "C"]]
    @test MCC_dict[Set(["a", "b"])] == MCC_dict[Set(["b", "c"])]
    @test SplitList(t1) == SplitList(input_trees[3])
    @test 0.0 == TreeKnit.consistency_rate(MCC_dict, input_trees)

    input_trees = [copy(t2), copy(t1), copy(t3)]
    MCC_dict = TreeKnit.infer_benchmark_MCCs!(input_trees, consistant=true)
    @test 0.0 == TreeKnit.consistency_rate(MCC_dict, input_trees)

    input_trees = [copy(t3), copy(t1), copy(t2)]
    MCC_dict = TreeKnit.infer_benchmark_MCCs!(input_trees, order="input", consistant=true)
    @test 0.0 == TreeKnit.consistency_rate(MCC_dict, input_trees)
end

@testset "benchmark simulated trees" begin
    MCC_dict = TreeKnit.infer_benchmark_MCCs(3, 10, debug=true, order="resolution")
    @test length(keys(MCC_dict)) == 4
    MCC_dict = TreeKnit.infer_benchmark_MCCs(3, 10, debug=true, remove=true, order="resolution")
    @test length(keys(MCC_dict)) == 4
end

