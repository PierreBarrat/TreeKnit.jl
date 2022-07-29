using Test
using TreeTools
using TreeKnit


t1 = read_tree("$(dirname(pathof(TreeKnit)))/../test/FluForest/tree_ha.nwk")
t1.label = "HA"
print(t1.label)
t2 = read_tree("$(dirname(pathof(TreeKnit)))/../test/FluForest/tree_na.nwk")
t2.label="NA"
t3 = read_tree("$(dirname(pathof(TreeKnit)))/../test/FluForest/tree_pb1.nwk")
t3.label = "PB1"
t4 = read_tree("$(dirname(pathof(TreeKnit)))/../test/FluForest/tree_pb2.nwk")
t4.label="PB2"
t1 = copy(convert(Tree{TreeTools.MiscData}, t1))
t2 = copy(convert(Tree{TreeTools.MiscData}, t2))
t3 = copy(convert(Tree{TreeTools.MiscData}, t3))
t4 = copy(convert(Tree{TreeTools.MiscData}, t4))

# MCCs = computeMCCs(t1, t2)
# MCC_dict = Dict(Set([t1.label, t2.label]) => MCCs)
# TreeKnit.draw_ARG([t1, t2], MCC_dict;
#     label_nodes = false, draw_connections = false) 
# TreeKnit.write_mccs("mcc.dat", MCCs)

MCC_dict = TreeKnit.infer_benchmark_MCCs!([t1, t2, t3], consistant=false)
TreeKnit.draw_ARG([t1, t2, t3], MCC_dict;
    label_nodes = false, draw_connections = false) 
TreeKnit.write_mccs("mcc13.dat", MCCs[Set(["HA", "PB1"])])
TreeKnit.write_mccs("mcc23.dat", MCCs[Set(["NA", "PB1"])])