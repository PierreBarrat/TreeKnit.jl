using TreeKnit.SRG
using TreeTools
using Test

# Lots of intricated singletons

nwk1 = "(((((A:2,X3)AX3:1,X4)AX4:2,X5)AX5:2,B)AB:0.,X1,X2):0.R"
nwk2 = "((((A:1,X2)AX2:3,X1)AX1:3,B)AB:0.,X3,X4,X5):0.R"

t1 = parse_newick_string(nwk1);
t2 = parse_newick_string(nwk2);

MCCs = [["A","B"],["X1"],["X2"],["X3"],["X4"],["X5"]]

X1, X2 = SRG.shared_nodes(t1, t2, MCCs)

tc1, tc2, Xc1, Xc2 = begin
       tc1 = copy(t1)
       tc2 = copy(t2)
       Xc1 = copy(X1)
       Xc2 = copy(X2)
       SRG.fix_shared_singletons!(tc1, tc2, Xc1, Xc2, MCCs)
       tc1, tc2, Xc1, Xc2
end

@testset "fix_shared_singletons!" begin
	@test Xc1["AX3"][1] == :shared
	@test Xc1["AX3"][2][1:9] == "Singleton"
	@test Xc2[Xc1["AX3"][2]][2] == "AX3"

	@test Xc1["AX4"][1] == :shared
	@test Xc1["AX4"][2][1:9] == "Singleton"
	@test Xc2[Xc1["AX4"][2]][2] == "AX4"

	@test Xc1["AX5"][1] == :shared
	@test Xc1["AX5"][2][1:9] == "Singleton"
	@test Xc2[Xc1["AX5"][2]][2] == "AX5"
	#
	@test Xc2["AX1"][1] == :shared
	@test Xc2["AX1"][2][1:9] == "Singleton"
	@test Xc1[Xc2["AX1"][2]][2] == "AX1"

	@test Xc2["AX2"][1] == :shared
	@test Xc2["AX2"][2][1:9] == "Singleton"
	@test Xc1[Xc2["AX2"][2]][2] == "AX2"

	@test divtime(t1.lnodes["A"], t1.root) == divtime(tc1.lnodes["A"], tc1.root)
	@test divtime(t1.lnodes["AX3"], t1.root) == divtime(tc1.lnodes["AX3"], tc1.root)
	@test divtime(t1.lnodes["AX4"], t1.root) == divtime(tc1.lnodes["AX4"], tc1.root)
	@test divtime(t1.lnodes["AX5"], t1.root) == divtime(tc1.lnodes["AX5"], tc1.root)

	@test divtime(t2.lnodes["A"], t2.root) == divtime(tc2.lnodes["A"], tc2.root)
	@test divtime(t2.lnodes["AX2"], t2.root) == divtime(tc2.lnodes["AX2"], tc2.root)
	@test divtime(t2.lnodes["AX1"], t2.root) == divtime(tc2.lnodes["AX1"], tc2.root)
end

# Weird branch lenghts
nwk1 = "((A:1,B:1)AB:1,C:2)R:0"
nwk2 = "(A:4,(B:0,C:0):4)R:0" # from a resolution

t1 = parse_newick_string(nwk1);
t2 = parse_newick_string(nwk2);

MCCs = [["A"], ["B", "C"]]

X1, X2 = SRG.shared_nodes(t1, t2, MCCs)

tc1, tc2, Xc1, Xc2 = begin
       tc1 = copy(t1)
       tc2 = copy(t2)
       Xc1 = copy(X1)
       Xc2 = copy(X2)
       SRG.fix_shared_singletons!(tc1, tc2, Xc1, Xc2, MCCs)
       tc1, tc2, Xc1, Xc2
end
