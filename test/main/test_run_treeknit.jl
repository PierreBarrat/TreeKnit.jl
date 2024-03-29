using Test
using TreeTools
using TreeKnit

println("##### testing MultiTreeKnit: get_infered_MCC_pairs! #####")

nwk1 = "((A,B),C);"
nwk2 = "(A,(B,C));"
nwk3 = "(A,B,C);"
t1 = parse_newick_string(nwk1, node_data_type=TreeTools.MiscData, label = "a")
t2 = parse_newick_string(nwk2, node_data_type=TreeTools.MiscData, label = "b")
t3 = parse_newick_string(nwk3, node_data_type=TreeTools.MiscData, label = "c")

solutions = [   [["A"], ["B", "C"]],
                [["B"], ["A", "C"]],
                [["C"], ["A", "B"]] 
            ]

@testset "get_infered_MCC_pairs!" begin
    input_trees = [copy(t1), copy(t2), copy(t3)]
    MCC_dict = run_treeknit!(input_trees, OptArgs(3; method= :better_MCCs))
    @test get(MCC_dict,("a", "c")) == [["A", "B", "C"]]
    @test in(get(MCC_dict, ("a", "b")), solutions)
	@test SplitList(t1) == SplitList(input_trees[3])

    input_trees = [copy(t1), copy(t2), copy(t3)]
    MCC_dict = run_treeknit!(input_trees, OptArgs(3; rounds=1, final_no_resolve=true))
    @test in(get(MCC_dict, ("a", "b")), solutions)
    @test in(get(MCC_dict, ("a", "c")), solutions)
	@test SplitList(t1) != SplitList(input_trees[3])

    input_trees = [copy(t2), copy(t1), copy(t3)]
    MCC_dict = run_treeknit!(input_trees, OptArgs(3; rounds=1, resolve=true, final_no_resolve=false))
    @test get(MCC_dict, ("b", "c")) == [["A", "B", "C"]]
    @test in(get(MCC_dict, ("a", "b")), solutions)
	@test SplitList(t2) == SplitList(input_trees[3])

    input_trees = [copy(t3), copy(t1), copy(t2)]
    MCC_dict = run_treeknit!(input_trees, OptArgs(3; rounds=1, resolve=true, final_no_resolve=false))
    @test get(MCC_dict,("a", "c")) == [["A", "B", "C"]]
    @test in(get(MCC_dict,("a", "b")), solutions)
	@test SplitList(t1) == SplitList(input_trees[1])

end


