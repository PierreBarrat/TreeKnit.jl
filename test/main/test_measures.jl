using Test
using TreeTools
using TreeKnit

println("##### measures #####")

nwk1 = "((A,B),C);"
nwk2 = "(A,(B,C));"
nwk3 = "((A,B,D),C);"
nwk4 = "(((A,B),D),C);"

t1 = node2tree(TreeTools.parse_newick(nwk1))
t2 = node2tree(TreeTools.parse_newick(nwk2))
t3 = node2tree(TreeTools.parse_newick(nwk3))
t4 = node2tree(TreeTools.parse_newick(nwk4))

@testset "RF distance" begin
	@test TreeKnit.RF_distance(t1, t2) ==2
    @test TreeKnit.RF_distance(t3, t4) ==1
end

MCC1 = [["1"], ["2"], ["3", "4", "5", "6"]]
MCC2 = [["1", "2"], ["3", "4", "5", "6"]]
MCC3 = [["1", "2", "3", "4", "5", "6"]]

@testset "is_degenerate" begin
	@test TreeKnit.is_degenerate(3, [MCC1, MCC1, MCC2]) ==false
    @test TreeKnit.is_degenerate(3, [MCC1, MCC1, MCC3]) ==false
    @test TreeKnit.is_degenerate(3, [MCC2, MCC2, MCC1]) ==true
end
