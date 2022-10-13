println("##### testing constraints #####")

function join_sets_slow(input_sets::Vector{Vector{Vector{String}}})
    start_set = input_sets[1]
    for i in 2:length(input_sets)
        joint_sets = Vector{String}[]
        to_be_joint_set = input_sets[i]
        for s1 in start_set
            for s2 in to_be_joint_set
                if !isdisjoint(s1, s2)
                    push!(joint_sets, intersect(s1, s2))
                end
            end
        end
        start_set = joint_sets
    end
    return TreeKnit.sort(start_set; lt=TreeKnit.clt)
end


tree3 = convert(Tree{TreeTools.MiscData},node2tree(parse_newick("((A,B)i1,((C,D)i2,((E,(F1,F2)i6)i4,G)i5)i3)i7")))
tree1 = convert(Tree{TreeTools.MiscData},node2tree(parse_newick("((A,B)j1,(((C,D)j2,(E,(F1,F2)j3)j4)j5,G)j6)j7")))
tree2 = convert(Tree{TreeTools.MiscData},node2tree(parse_newick("((G,(A,B)k1)k2,((E,(C,D)k3)k4,(F1,F2)k5)k6)k7")))
MCC12 = TreeKnit.sort([["A", "B", "E", "F1", "F2"], ["G"], ["C", "D"]], lt=TreeKnit.clt)
MCC13 = TreeKnit.sort([["A", "B", "C", "D", "G"], ["E", "F1", "F2"]], lt=TreeKnit.clt)

constraint = TreeKnit.join_sets([MCC12, MCC13])

@testset "MCC join" begin
    @test join_sets_slow([MCC12, MCC13]) ==  TreeKnit.join_sets([MCC12, MCC13])
    @test constraint == [["G"], ["A", "B"], ["C", "D"], ["E", "F1", "F2"]]
end

@testset "mcc_map functions" begin
    mcc_map, cluster_no = TreeKnit.leaf_mcc_map(constraint, get_cluster_no =true)
    @test mcc_map == Dict("B" => 2, "A" => 2, "C" => 3, "D" => 3, "G" => 1, "E" => 4, "F1" => 4, "F2" => 4)
    mcc_double_map = TreeKnit.leaf_mcc_map([MCC12, MCC13])
    @test mcc_double_map == Dict("B" => [3, 2], "A" => [3, 2], "C" => [2, 2], "D" => [2, 2], "G" => [1, 2], "E" => [3, 1], "F1" => [3, 1], "F2" => [3, 1])
    TreeKnit.assign_mccs!(tree3, mcc_map)
    inner_node_dict = Dict("i1" => 2, "i2" => 3, "i3" => nothing, "i4" => 4, "i5" => nothing, "i6" => 4, "i7" => nothing)
    for node in nodes(tree3)
        if isleaf(node)
            @test node.data["mcc"] == mcc_map[node.label]
        else
            @test node.data["mcc"] == inner_node_dict[node.label]
        end
    end
end

@testset "mark_shared_branches functions" begin
    TreeKnit.mark_shared_branches!(constraint, tree3)
    true_labels = Set(["A", "B", "C", "D", "E", "F1", "F2", "i6"])
    for node in nodes(tree3)
        if node.label in true_labels
            @test node.data["shared_branch_constraint"] == true
        else
            @test node.data["shared_branch_constraint"] == false
        end
    end
end