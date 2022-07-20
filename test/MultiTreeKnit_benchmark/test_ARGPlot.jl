using Test
using TreeTools
using TreeKnit

println("##### testing ARGPlot #####")

nwk1 = "((A,B),C)R;"
nwk2 = "(A,(B,C))R;"
nwk3 = "(A,B,C)R;"
t1 = node2tree(TreeTools.parse_newick(nwk1))
t2 = node2tree(TreeTools.parse_newick(nwk2))
t3 = node2tree(TreeTools.parse_newick(nwk3))


TreeKnit.draw_ARG(
        t1, [],
        Dict()
    )
mcc_map = TreeKnit.get_mcc_map([[["A"], ["B", "C"]]])
ot1 = copy(convert(Tree{TreeTools.MiscData}, t1))
ot2 = copy(convert(Tree{TreeTools.MiscData}, t2))
TreeKnit.assign_mccs!(TreeKnit.get_mcc_map([["A"], ["B", "C"]]), [ot2])
TreeKnit.assign_all_mccs!(ot1, 1, mcc_map)
rec_dict = TreeKnit.get_recombination_sites(ot1,[ot2], [[["A"], ["B", "C"]]])

TreeKnit.draw_ARG(
        ot1, [ot2],
        rec_dict
    )