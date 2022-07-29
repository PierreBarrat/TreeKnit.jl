using Test
using TreeTools
using TreeKnit


t1 = read_tree("$(dirname(pathof(TreeKnit)))/../test/FluForest/tree_ha.nwk")
t2 = read_tree("$(dirname(pathof(TreeKnit)))/../test/FluForest/tree_na.nwk")
t1 = copy(convert(Tree{TreeTools.MiscData}, t1), label="HA")
t2 = copy(convert(Tree{TreeTools.MiscData}, t2), label="NA")

MCCs = computeMCCs(t1, t2)
MCC_dict = Dict(Set([t1.label, t2.label]) => MCCs)
TreeKnit.draw_ARG([t1, t2], MCC_dict;
    label_nodes = true, draw_connections = false) 
print(MCCs)