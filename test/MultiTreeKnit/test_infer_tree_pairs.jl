using Test
using TreeTools
using TreeKnit
using TreeKnit.MTK

println("##### testing MultiTreeKnit: get_infered_MCC_pairs! #####")

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

@testset "get_infered_MCC_pairs!, consistent=false" begin
    input_trees = [copy(t1), copy(t2), copy(t3)]
    MCC_dict = MTK.get_infered_MCC_pairs!(input_trees, OptArgs(rounds=1, consistent=false))
    @test get(MCC_dict,("a", "c")) == [["A", "B", "C"]]
    @test in(get(MCC_dict, ("a", "b")), solutions)
	@test SplitList(t1) == SplitList(input_trees[3])

    input_trees = [copy(t2), copy(t1), copy(t3)]
    MCC_dict = MTK.get_infered_MCC_pairs!(input_trees, OptArgs(rounds=1, consistent=false))
    @test get(MCC_dict, ("b", "c")) == [["A", "B", "C"]]
    @test in(get(MCC_dict, ("a", "b")), solutions)
	@test SplitList(t2) == SplitList(input_trees[3])

    input_trees = [copy(t3), copy(t1), copy(t2)]
    MCC_dict = MTK.get_infered_MCC_pairs!(input_trees, OptArgs(rounds=1, consistent=false))
    @test get(MCC_dict,("a", "c")) == [["A", "B", "C"]]
    @test in(get(MCC_dict,("a", "b")), solutions)
	@test SplitList(t1) == SplitList(input_trees[1])

end


@testset "consistency shared_branch" begin
    input_trees = [copy(t1), copy(t2), copy(t3)]
    ot3 = copy(convert(Tree{TreeTools.MiscData}, input_trees[1]))
    ot2 = copy(convert(Tree{TreeTools.MiscData}, input_trees[2]))
    MTK.map_shared_branches!([["A"], ["B", "C"]], ot2)
    MTK.map_shared_branches!([["A"], ["B", "C"]], ot3)
    @test ot2.lnodes["NODE_1"].data.dat["shared_branch"] == false
    @test ot3.lnodes["NODE_1"].data.dat["shared_branch"] == true
    @test ot2.lleaves["A"].data.dat["shared_branch"] == false
    @test ot3.lleaves["A"].data.dat["shared_branch"] == false
    @test ot2.lleaves["B"].data.dat["shared_branch"] == true
    @test ot3.lleaves["B"].data.dat["shared_branch"] == true

    input_trees = [copy(t1), copy(t2), copy(t3)]
    for constraint in solutions
        mccs = TreeKnit.runopt(TreeKnit.OptArgs(consistent=true), input_trees[1], input_trees[2], constraint; output = :mccs)
        @test mccs == constraint
    end

end

@testset "get_infered_MCC_pairs!, consistent=true" begin
    input_trees = [copy(t1), copy(t2), copy(t3)]
    MCC_dict = MTK.get_infered_MCC_pairs!(input_trees, OptArgs(rounds=1, consistent=true))
    @test get(MCC_dict,("a", "c")) == [["A", "B", "C"]]
    @test get(MCC_dict, ("a", "b")) == get(MCC_dict, ("b", "c"))
    @test SplitList(t1) == SplitList(input_trees[3])
    constraint = MTK.MCC_join_constraint([get(MCC_dict,("a", "c")), get(MCC_dict,("a", "b"))])
    @test MTK.is_MCC_subset(constraint, get(MCC_dict,("b", "c")))

    input_trees = [copy(t2), copy(t1), copy(t3)]
    MCC_dict = MTK.get_infered_MCC_pairs!(input_trees, OptArgs(rounds=1, consistent=true))
    constraint = MTK.MCC_join_constraint([get(MCC_dict,("a", "c")), get(MCC_dict,("a", "b"))])
    @test MTK.is_MCC_subset(constraint, get(MCC_dict,("b", "c")))

    input_trees = [copy(t3), copy(t1), copy(t2)]
    MCC_dict = MTK.get_infered_MCC_pairs!(input_trees, OptArgs(rounds=1, consistent=true))
    constraint = MTK.MCC_join_constraint([get(MCC_dict,("a", "c")), get(MCC_dict,("a", "b"))])
    @test MTK.is_MCC_subset(constraint, get(MCC_dict,("b", "c")))
end

