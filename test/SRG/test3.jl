using TreeKnit
using TreeKnit.SRG
using Test
using TreeTools

nwk1 = "(A:3,(B:2,C:2)BC:1)R;"
nwk2 = "((A:3,B:3)AB:1,C:4)R;"

t1 = parse_newick_string(nwk1);
t2 = parse_newick_string(nwk2);

MCCs = [["A","B"], ["C"]]

X1, X2 = SRG.shared_nodes(t1, t2, MCCs)
SRG.fix_shared_singletons!(t1, t2, X1, X2, MCCs)
