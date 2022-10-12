using TreeKnit
using TreeKnit.SplitGraph
using Test
using TreeTools

println("##### Detect topological compatibility Issues #####")

t1 = node2tree(TreeTools.parse_newick("((((E1,E2),(B1,(B2,B3))),(C1,C2)),D1,D2)", node_data_type=TreeTools.MiscData))
t2 = node2tree(TreeTools.parse_newick("((B1,(B2,B3)),C1,C2,((E1,E2),D1,D2))", node_data_type=TreeTools.MiscData))


@testset "MCC compatible up to resolution" begin
    mcc = ["E1", "E2", "B1", "B2", "C1", "C2"]
    @test TreeKnit.is_topo_compatible(t1, t2, mcc) == true
    @test TreeKnit.is_full_topo_compatible(t1, t2, mcc) == false
    @test naive_mccs(t1, t2, mcc) == [["C1"], ["C2"], ["B1", "B2"], ["E1", "E2"]]
    @test naive_mccs(t1, t2, mcc; resolve=true) == [["B1", "B2", "C1", "C2", "E1", "E2"]]
end

@testset "MCC not compatible up to resolution" begin
    mcc = ["E1", "E2", "B1", "B2", "B3", "D1", "D2"]
    @test TreeKnit.is_topo_compatible(t1, t2, mcc) == false
    @test TreeKnit.is_full_topo_compatible(t1, t2, mcc) == false
    @test naive_mccs(t1, t2, mcc) ==  [["D1"], ["D2"], ["E1", "E2"], ["B1", "B2", "B3"]]
    @test naive_mccs(t1, t2, mcc; resolve=true) == [["D1"], ["D2"], ["E1", "E2"], ["B1", "B2", "B3"]]
end

@testset "MCC always compatible up to resolution" begin
    mcc = ["D1", "D2", "E1", "E2"]
    @test TreeKnit.is_topo_compatible(t1, t2, mcc) == true
    @test TreeKnit.is_full_topo_compatible(t1, t2, mcc) == true
    @test naive_mccs(t1, t2, mcc) ==  [["D1", "D2", "E1", "E2"]]
    @test naive_mccs(t1, t2, mcc; resolve=true) == [["D1", "D2", "E1", "E2"]]
end

t1_pre_TK = node2tree(TreeTools.parse_newick("((I,(A1,A2,B,C,D)),(E,F,G,H))", node_data_type=TreeTools.MiscData))
t1 = node2tree(TreeTools.parse_newick("((I,(A2,(A1,B),(C,D))),(E,F,(G,H)))", node_data_type=TreeTools.MiscData))
t2_pre_TK = node2tree(TreeTools.parse_newick("((A1,A2,B,C,D),((E,F,G,H),I))", node_data_type=TreeTools.MiscData))
t2 = node2tree(TreeTools.parse_newick("((A1,B,(A2,C),D),((F,(E,G),H),I))", node_data_type=TreeTools.MiscData))

@testset "MCC compatible up to resolution 2" begin
    mcc = ["A1", "A2", "B", "C"]
    @test TreeKnit.is_topo_compatible(t1, t2, mcc) == true
    @test TreeKnit.is_full_topo_compatible(t1, t2, mcc) == false
    @test naive_mccs(t1, t2, mcc) ==  [["A1"], ["A2"], ["B"], ["C"]]
    @test naive_mccs(t1, t2, mcc; resolve=true) == [["A1", "A2", "B", "C"]]
end

@testset "MCC not compatible up to resolution 2" begin
    mcc = ["E", "F", "G", "H"]
    @test TreeKnit.is_topo_compatible(t1, t2, mcc) == false
    @test TreeKnit.is_full_topo_compatible(t1, t2, mcc) == false
    @test naive_mccs(t1, t2, mcc) ==  [["E"], ["F"], ["G"], ["H"]]
    @test naive_mccs(t1, t2, mcc; resolve=true) == [["E"], ["F"], ["G"], ["H"]]
end

@testset "MCC not compatible up to resolution 3" begin
    mcc = ["A1", "A2", "B", "C", "D", "E", "F", "G", "H"]
    @test TreeKnit.is_topo_compatible(t1, t2, mcc) == false
    @test TreeKnit.is_full_topo_compatible(t1, t2, mcc) == false
    @test naive_mccs(t1, t2, mcc) ==  [["A1"], ["A2"], ["B"], ["C"], ["D"], ["E"], ["F"], ["G"], ["H"]]
    @test naive_mccs(t1, t2, mcc; resolve=true) == [["A1"], ["A2"], ["B"], ["C"], ["D"], ["E"], ["F"], ["G"], ["H"]]
end

@testset "topo_fix_mccs" begin
    MCCs_input = [["A1", "A2", "B", "C", "D", "E", "F", "G", "H"], ["I"]]
    MCCs = TreeKnit.topo_fix_mccs(t1, t2, MCCs_input; resolve=false)
    MCCs_true = TreeKnit.topo_fix_mccs(t1, t2, MCCs_input; resolve=true)
    @test MCCs == MCCs_true == [["A1"], ["A2"], ["B"], ["C"], ["D"], ["E"], ["F"], ["G"], ["H"], ["I"]]
    MCCs_input = [["A1", "A2", "B", "C"], ["D"], ["E", "F", "G", "H"], ["I"]]
    MCCs = TreeKnit.topo_fix_mccs(t1, t2, MCCs_input; resolve=false)
    MCCs_true = TreeKnit.topo_fix_mccs(t1, t2, MCCs_input; resolve=true)
    @test MCCs == [["A1"], ["A2"], ["B"], ["C"], ["D"], ["E"], ["F"], ["G"], ["H"], ["I"]]
    @test MCCs_true == [["D"], ["E"], ["F"], ["G"], ["H"], ["I"], ["A1", "A2", "B", "C"]]
end