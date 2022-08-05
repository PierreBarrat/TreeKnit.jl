using Test
using TreeTools
using TreeKnit


t1 = read_tree("$(dirname(pathof(TreeKnit)))/../test/FluForest_small/tree_ha_6y.nwk")
t2 = read_tree("$(dirname(pathof(TreeKnit)))/../test/FluForest_small/tree_na_6y.nwk")
t3 = read_tree("$(dirname(pathof(TreeKnit)))/../test/FluForest_small/tree_pb1_6y.nwk")
t4 = read_tree("$(dirname(pathof(TreeKnit)))/../test/FluForest_small/tree_pb2_6y.nwk")
t1 = copy(convert(Tree{TreeTools.MiscData}, t1))
t2 = copy(convert(Tree{TreeTools.MiscData}, t2))
t3 = copy(convert(Tree{TreeTools.MiscData}, t3))
t4 = copy(convert(Tree{TreeTools.MiscData}, t4))
t1.label = "HA"
t2.label="NA"
t3.label = "PB1"
t4.label="PB2"

# MCCs = computeMCCs(t1, t2, TreeKnit.OptArgs(itmax=50, nT= 3000, γ = 2.25))
# #lowering the γ to 1.75 leads to a slightly higher number of MCCs that are found
# MCC_dict = Dict(Set([t1.label, t2.label]) => MCCs)
# #TreeKnit.draw_ARG([t1, t2], MCC_dict;
# #    label_nodes = false, draw_connections = false) 
# TreeKnit.write_mccs("mcc12.dat", MCCs)

MCC_dict = TreeKnit.infer_benchmark_MCCs!([t1, t2, t3], consistant=false, rounds=1)
#TreeKnit.draw_ARG([t1, t2, t3], MCC_dict;
#    label_nodes = false, draw_connections = false) 
TreeKnit.write_mccs("mcc12.dat", MCC_dict[Set(["HA", "NA"])])
TreeKnit.write_mccs("mcc13.dat", MCC_dict[Set(["HA", "PB1"])])
TreeKnit.write_mccs("mcc23.dat", MCC_dict[Set(["NA", "PB1"])])

MCC_dict = TreeKnit.infer_benchmark_MCCs!([t1, t2, t3], consistant=true, force=true, force_rounds=20)
#TreeKnit.draw_ARG([t1, t2, t3], MCC_dict;
#    label_nodes = false, draw_connections = false) 
TreeKnit.write_mccs("mcc12_c.dat", MCC_dict[Set(["HA", "NA"])])
TreeKnit.write_mccs("mcc13_c.dat", MCC_dict[Set(["HA", "PB1"])])
TreeKnit.write_mccs("mcc23_c.dat", MCC_dict[Set(["NA", "PB1"])])