using Test
using TreeKnit
using TreeKnit.SRG
using TreeTools


# nwkA= "((A,B),C);"
# nwkB = "(A,(B,C));"
# nwkC = nwkA

# tA = node2tree(TreeTools.parse_newick(nwkA), label = "a")
# tB = node2tree(TreeTools.parse_newick(nwkB), label = "b")
# tC = node2tree(TreeTools.parse_newick(nwkC), label = "c")
# trees = [tA, tB, tC]
# MCCs = TreeKnit.MCC_set(3, ["a", "b", "c"], Dict(Set(["a", "b"]) => [["A", "B"],["C"]], Set(["a", "c"]) => [["A", "B","C"]], Set(["b", "c"]) => [["A", "B"],["C"]]))

# arg, lm_dict = SRG.arg_from_trees(trees, MCCs)

# [length(nodes(t)) for t in trees] == [7,7,7]


nwkA= "((A,B),(C,D));"
nwkB = "(((A,B),C),D);"
nwkC = "(((A,B),D),C);"

tA = node2tree(TreeTools.parse_newick(nwkA), label = "a")
tB = node2tree(TreeTools.parse_newick(nwkB), label = "b")
tC = node2tree(TreeTools.parse_newick(nwkC), label = "c")
trees = [tA, tB, tC]
MCCs = TreeKnit.MCC_set(3, ["a", "b", "c"], Dict(Set(["a", "b"]) => [["A", "B", "D"],["C"]], Set(["a", "c"]) => [["A", "B","C"], ["D"]], Set(["b", "c"]) => [["A", "B"],["C", "D"]]))

arg, lm_dict = SRG.arg_from_trees(trees, MCCs)
