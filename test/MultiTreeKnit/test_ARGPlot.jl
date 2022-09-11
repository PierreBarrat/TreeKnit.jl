using Test
using TreeTools
using TreeKnit

println("##### testing ARGPlot #####")

nwk1 = "((A,B),C)R;"
nwk2 = "(A,(B,C))R;"
nwk3 = "(A,B,C)R;"
t1 = node2tree(TreeTools.parse_newick(nwk1, node_data_type=TreeTools.MiscData), label = "a")
t2 = node2tree(TreeTools.parse_newick(nwk2, node_data_type=TreeTools.MiscData), label= "b")
t3 = node2tree(TreeTools.parse_newick(nwk3, node_data_type=TreeTools.MiscData), label= "c")


TreeKnit.draw_ARG(
        copy(convert(Tree{TreeTools.MiscData}, t1)), nothing,
        nothing
)

mcc_map = TreeKnit.leaf_mcc_map([[["A"], ["B", "C"]]])
ot1 = copy(t1)
ot2 = copy(t2)
TreeKnit.assign_mccs!([ot2], TreeKnit.leaf_mcc_map([["A"], ["B", "C"]]))
TreeKnit.assign_all_mccs!(ot1, 1, mcc_map)
rec_dict = TreeKnit.get_recombination_sites(ot1,[ot2], MCC_set(2, ["a", "b"], Dict(Set(["a", "b"]) =>[["A"], ["B", "C"]])))

TreeKnit.draw_ARG(
        ot1, [ot2],
        rec_dict
    )

TreeKnit.draw_ARG(
        ot1, [ot2],
        rec_dict, draw_connections = true
    )
ot1 = copy(convert(Tree{TreeTools.MiscData}, t1))
label!(ot1, "a")
ot2 = copy(convert(Tree{TreeTools.MiscData}, t2))
label!(ot2, "b")
MCCs_dict = MCC_set(2, ["a", "b"], Dict(Set(["a", "b"]) => [["A"], ["B", "C"]]))

draw_ARG([ot1, ot2], MCCs_dict; draw_connections = true) 

t1 = read_tree("$(dirname(pathof(TreeKnit)))/../test/NYdata/tree_ha.nwk")
t2 = read_tree("$(dirname(pathof(TreeKnit)))/../test/NYdata/tree_na.nwk")
ot1 = copy(convert(Tree{TreeTools.MiscData}, t1))
ot2 = copy(convert(Tree{TreeTools.MiscData}, t2))

MCCs = computeMCCs(ot1, ot2)
MCCs_dict = MCC_set(2, ["ha", "na"], Dict(Set(["ha", "na"]) => MCCs))
label!(ot1, "ha")
label!(ot2, "na")
draw_ARG([ot1, ot2], MCCs_dict; label_nodes = false, draw_connections = true) 