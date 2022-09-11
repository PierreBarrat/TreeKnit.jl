using TreeKnit, TreeTools

nwk_a = "(A,(B,C));"
nwk_b = "((A,B),C);"
nwk_c = "((A,C),B);"

t_a = node2tree(TreeTools.parse_newick(nwk_a, node_data_type=TreeTools.MiscData), label = "a")
t_b = node2tree(TreeTools.parse_newick(nwk_b, node_data_type=TreeTools.MiscData), label= "b")
t_c = node2tree(TreeTools.parse_newick(nwk_c, node_data_type=TreeTools.MiscData), label= "c")


##infered MCCs
MCC_ab = [["C"], ["A","B"]]
MCC_ac = [["B"], ["A", "C"]]
MCC_bc = [["A"], ["B", "C"]]

##as ["C"] is a unequal subset of ["A", "C"] "C" must have a reassortment event above it between trees b and c, the same follows for "B" and ["A", "B"]

pair_MCCs = TreeKnit.convert_MCC_list_to_set(3, ["a", "b", "c"], [MCC_ab, MCC_ac, MCC_bc])
trees = [t_a, t_b, t_c]
b = TreeKnit.is_topologically_degenerate(pair_MCCs, trees)
@assert b == true

##infered MCCs
MCC_ab = [["C"], ["A","B"]]
MCC_ac = [["B"], ["A", "C"]]
MCC_bc = [["A"], ["B"], ["C"]]

pair_MCCs = TreeKnit.convert_MCC_list_to_set(3, ["a", "b", "c"], [MCC_ab, MCC_ac, MCC_bc])
trees = [t_a, t_b, t_c]
b = TreeKnit.is_topologically_degenerate(pair_MCCs, trees)
@assert b == true

##infered MCCs
MCC_ab = [["C"], ["A","B"]]
MCC_ac = [["B"], ["A"], ["C"]]
MCC_bc = [["A"], ["B"], ["C"]]

pair_MCCs = TreeKnit.convert_MCC_list_to_set(3, ["a", "b", "c"], [MCC_ab, MCC_ac, MCC_bc])
trees = [t_a, t_b, t_c]
b = TreeKnit.is_topologically_degenerate(pair_MCCs, trees)
@assert b == false