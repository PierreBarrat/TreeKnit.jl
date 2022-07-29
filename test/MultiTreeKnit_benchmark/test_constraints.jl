using TreeKnit
using TreeTools
using Test

println("##### testing constraints #####")

function join_sets_slow(input_sets::Vector{Vector{Vector{String}}})
    start_set = input_sets[1]
    for i in 2:length(input_sets)
        joint_sets = Vector{String}[]
        to_be_joint_set = [Set{String}(s2) for s2 in input_sets[i]]
        for s1 in start_set
            nodes = length(s1)
            while nodes>0
                for s2 in to_be_joint_set
                    joint_set = intersect(Set{String}(s1), s2)
                    if !isempty(joint_set)
                        s2 = setdiff(s2,joint_set)
                        nodes -= length(joint_set)
                        append!(joint_sets, [sort(collect(joint_set))])
                    end
                end
            end
        end
        start_set = joint_sets
    end
    return TreeKnit.sort(start_set; lt=TreeKnit.clt)
end

trees, arg = TreeKnit.get_trees(3, 100)
rMCC_list = TreeKnit.get_real_MCCs(3, arg)
slow = join_sets_slow([rMCC_list[1], rMCC_list[2]])
fast = TreeKnit.join_sets([rMCC_list[1], rMCC_list[2]])

tree3 = convert(Tree{TreeTools.MiscData},node2tree(parse_newick("((A,B)i1,((C,D)i2,((E,(F1,F2)i6)i4,G)i5)i3)i7")))
tree1 = convert(Tree{TreeTools.MiscData},node2tree(parse_newick("((A,B)j1,(((C,D)j2,(E,(F1,F2)j3)j4)j5,G)j6)j7")))
tree2 = convert(Tree{TreeTools.MiscData},node2tree(parse_newick("((G,(A,B)k1)k2,((E,(C,D)k3)k4,(F1,F2)k5)k6)k7")))
MCC12 = TreeKnit.sort([["A", "B", "E", "F1", "F2"], ["G"], ["C", "D"]], lt=TreeKnit.clt)
MCC13 = TreeKnit.sort([["A", "B", "C", "D", "G"], ["E", "F1", "F2"]], lt=TreeKnit.clt)

constraint = TreeKnit.join_sets([MCC12, MCC13])

@testset "MCC join" begin
    @test slow == fast
    @test constraint == [["G"], ["A", "B"], ["C", "D"], ["E", "F1", "F2"]]
end

@testset "mcc_map functions" begin
    mcc_map, cluster_no = TreeKnit.get_mcc_map(constraint, get_cluster_no =true)
    @test mcc_map == Dict("B" => 2, "A" => 2, "C" => 3, "D" => 3, "G" => 1, "E" => 4, "F1" => 4, "F2" => 4)
    mcc_double_map = TreeKnit.get_mcc_map([MCC12, MCC13])
    @test mcc_double_map == Dict("B" => [3, 2], "A" => [3, 2], "C" => [2, 2], "D" => [2, 2], "G" => [1, 2], "E" => [3, 1], "F1" => [3, 1], "F2" => [3, 1])
    TreeKnit.assign_mccs!(mcc_map, tree3)
    inner_node_dict = Dict("i1" => 2, "i2" => 3, "i3" => nothing, "i4" => 4, "i5" => nothing, "i6" => 4, "i7" => nothing)
    for node in nodes(tree3)
        if isleaf(node)
            @test node.data["mcc"] == mcc_map[node.label]
        else
            @test node.data["mcc"] == inner_node_dict[node.label]
        end
    end
end

@testset "add_mask functions" begin
    TreeKnit.add_mask!(constraint, tree3)
    true_labels = Set(["A", "B", "C", "D", "E", "F1", "F2", "i6"])
    for node in nodes(tree3)
        if node.label in true_labels
            @test node.data["mask"] == true
        else
            @test node.data["mask"] == false
        end
    end
end

@testset "fix consistency functions" begin
    MCC23 = TreeKnit.sort([["A", "B", "C", "D", "G"], ["E"], ["F1", "F2"]], lt=TreeKnit.clt)
    tree3 = convert(Tree{TreeTools.MiscData},node2tree(parse_newick("((A,B)i1,((C,D)i2,((E,(F1,F2)i6)i4,G)i5)i3)i7")))
    MCCs = TreeKnit.fix_consist!([MCC12, MCC13, MCC23], [tree2, tree3])
    @test MCCs[3] == TreeKnit.sort([["A", "B", "C", "D", "G"], ["E", "F1", "F2"]], lt=TreeKnit.clt)


    MCC23 = TreeKnit.sort([["A", "B", "G"], ["C", "D", "E"], ["F1", "F2"]], lt=TreeKnit.clt)
    tree3 = convert(Tree{TreeTools.MiscData},node2tree(parse_newick("((A,B)i1,((C,D)i2,((E,(F1,F2)i6)i4,G)i5)i3)i7")))
    MCCs = TreeKnit.fix_consist!([MCC12, MCC13, MCC23], [tree2, tree3])
    @test MCCs[3] == MCC23
    @test MCCs[1] in [MCC12, [["G"], ["C", "D"], ["F1", "F2"], ["A", "B", "E"]], [["E"], ["G"], ["A", "B"], ["C", "D"], ["F1", "F2"]]]
    @test MCCs[2] in [MCC13, [["E"], ["F1", "F2"], ["A", "B", "C", "D", "G"]]]
end