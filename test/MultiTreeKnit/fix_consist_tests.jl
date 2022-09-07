using TreeTools
using TreeKnit


# #nwk1 = "((2_0:3194.443263974577,(5_0:3194.443263974577,(6_0:29.449927592527846,1_0:29.449927592527846)internal_1:3164.993336382049)RESOLVED_1:0.0)internal_11:39119.403233106874,(3_0:6819.317781996674,((7_0:1013.66826726748,9_0:1013.66826726748)RESOLVED_2:0.0,(10_0:1013.66826726748,(8_0:1013.66826726748,4_0:1013.66826726748)RESOLVED_3:0.0)RESOLVED_4:0.0)internal_5:5805.6495147291935)internal_16:35494.52871508478)internal_21:0;"
# #nwk1 = "((6_0:9311.428208782838,(5_0:5784.485018822639,(10_0:960.5199719022602,3_0:960.5199719022602)internal_4:4823.965046920379)internal_13:3526.9431899602005)internal_14:23043.676758974947,((8_0:137.58973093427062,9_0:137.58973093427062)internal_2:4663.72096306566,((7_0:4801.310693999931,1_0:4801.310693999931)RESOLVED_2:0.0,(2_0:4801.310693999931,4_0:4801.310693999931)RESOLVED_1:0.0)RESOLVED_3:0.0)internal_12:27553.794273757856)internal_22:0;"
# nwk1 = "(((10_0:144.16922758232607,5_0:144.16922758232607)internal_1:2542.039081117728,(6_0:2686.208308700054,(8_0:2686.208308700054,2_0:2686.208308700054)RESOLVED_2:0.0)RESOLVED_3:0.0)internal_10:21259.572478349197,((9_0:152.12177046535655,3_0:152.12177046535655)internal_2:3912.1950487985705,(7_0:4064.316819263927,(4_0:615.1812363058868,1_0:615.1812363058868)internal_4:3449.13558295804)RESOLVED_1:0.0)internal_12:19881.463967785323)internal_16:0;"
# ta = node2tree(TreeTools.parse_newick(nwk1, node_data_type=TreeTools.MiscData), label = "a")
# # MCC_ab = [["3_0"], ["1_0", "5_0", "6_0"], ["10_0", "2_0", "4_0", "7_0", "8_0", "9_0"]]
# # MCC_ac = [["1_0", "6_0"], ["2_0", "5_0"], ["10_0", "3_0", "4_0", "7_0", "8_0", "9_0"]]
# # MCC_bc = [["1_0", "6_0"], ["2_0", "3_0", "5_0"], ["10_0", "4_0", "7_0", "8_0", "9_0"]]
# # MCC_ab =  [["10_0", "3_0", "5_0", "6_0"], ["1_0", "2_0", "4_0", "7_0", "8_0", "9_0"]]
# # MCC_ac = [["6_0"], ["1_0", "2_0", "4_0", "7_0"], ["10_0", "3_0", "5_0", "8_0", "9_0"]]
# # MCC_bc = [["6_0"], ["1_0", "2_0", "4_0", "7_0"], ["10_0", "3_0", "5_0", "8_0", "9_0"]]
# MCC_ab = [["6_0"], ["10_0", "5_0"], ["1_0", "2_0", "3_0", "4_0", "7_0", "8_0", "9_0"]]
# MCC_ac =  [["2_0", "6_0", "8_0"], ["10_0", "1_0", "3_0", "4_0", "5_0", "7_0", "9_0"]]
# MCC_bc = [["2_0", "8_0"], ["10_0", "1_0", "3_0", "4_0", "5_0", "6_0", "7_0", "9_0"]]


# iMCCs = TreeKnit.fix_consist!([MCC_ac, MCC_bc, MCC_ab], [ta, ta], merge=false)

# print(iMCCs)

# ##trees given to TreeKnit as input 
# t_a = "(5_0:1.0,3_0:1.0,((8_0:1.0,13_0:1.0,2_0:1.0)internal_16:1.0,11_0:1.0,1_0:1.0,10_0:1.0)internal_20:1.0,12_0:1.0,7_0:1.0,9_0:1.0,15_0:1.0,14_0:1.0,6_0:1.0,4_0:1.0)internal_24:1.0;"
# t_b = "(((14_0:1.0,15_0:1.0)internal_7:1.0,12_0:1.0,7_0:1.0,9_0:1.0,6_0:1.0,1_0:1.0,4_0:1.0,13_0:1.0)internal_14:1.0,5_0:1.0,3_0:1.0,8_0:1.0,2_0:1.0,(10_0:1.0,11_0:1.0)internal_4:1.0)internal_24:1.0;"
# t_c = "(((5_0:1.0,3_0:1.0)internal_21:1.0,(8_0:1.0,2_0:1.0)internal_16:1.0,11_0:1.0,1_0:1.0,10_0:1.0)internal_22:1.0,(12_0:1.0,7_0:1.0,9_0:1.0,4_0:1.0,13_0:1.0,6_0:1.0)internal_11:1.0,15_0:1.0,14_0:1.0)internal_24:1.0;"

# ##output trees after running recursively
# nwk_a_output = "((5_0:1.0,3_0:1.0)RESOLVED_4:1.0,((14_0:1.0,15_0:1.0)RESOLVED_2:1.0,(9_0:1.0,7_0:1.0,6_0:1.0,4_0:1.0,12_0:1.0)RESOLVED_6:1.0,((10_0:1.0,11_0:1.0)RESOLVED_1:1.0,(1_0:1.0,(13_0:1.0,(8_0:1.0,2_0:1.0)RESOLVED_5:1.0)internal_16:1.0)RESOLVED_7:1.0)internal_20:1.0)RESOLVED_3:1.0)internal_24:1.0;"
# nwk_b_output = "(((8_0:1.0,2_0:1.0)RESOLVED_3:1.0,(10_0:1.0,11_0:1.0)internal_4:1.0)RESOLVED_6:1.0,((5_0:1.0,3_0:1.0)RESOLVED_4:1.0,((14_0:1.0,15_0:1.0)internal_7:1.0,((13_0:1.0,1_0:1.0)RESOLVED_1:1.0,(9_0:1.0,7_0:1.0,6_0:1.0,4_0:1.0,12_0:1.0)RESOLVED_5:1.0)RESOLVED_7:1.0)internal_14:1.0)RESOLVED_2:1.0)internal_24:1.0;"
# nwk_c_output = "((((1_0:1.0,(8_0:1.0,2_0:1.0)internal_16:1.0)RESOLVED_5:1.0,(10_0:1.0,11_0:1.0)RESOLVED_2:1.0)RESOLVED_3:1.0,(5_0:1.0,3_0:1.0)internal_21:1.0)internal_22:1.0,((14_0:1.0,15_0:1.0)RESOLVED_1:1.0,(13_0:1.0,(9_0:1.0,7_0:1.0,6_0:1.0,4_0:1.0,12_0:1.0)RESOLVED_4:1.0)internal_11:1.0)RESOLVED_6:1.0)internal_24:1.0;"
# t_a_output = node2tree(TreeTools.parse_newick(nwk_a_output, node_data_type=TreeTools.MiscData), label = "a")
# t_b_output = node2tree(TreeTools.parse_newick(nwk_b_output, node_data_type=TreeTools.MiscData), label = "b")
# t_c_output = node2tree(TreeTools.parse_newick(nwk_c_output, node_data_type=TreeTools.MiscData), label = "c")
# ##infered MCCs, 13 and 1 need to be removed in MCC_ab or 2,8 and 10, 11 split
# MCC_ab = [["2_0", "8_0", "10_0", "11_0"], ["1_0", "3_0", "4_0", "5_0", "6_0", "7_0", "9_0", "15_0", "14_0", "13_0", "12_0"]]
# MCC_ac = [["13_0"], ["3_0", "5_0"], ["10_0", "11_0", "12_0", "14_0", "15_0", "1_0", "2_0", "4_0", "6_0", "7_0", "8_0", "9_0"]]
# MCC_bc = [["1_0"], ["10_0", "11_0", "2_0", "3_0", "5_0", "8_0"], ["12_0", "13_0", "14_0", "15_0", "4_0", "6_0", "7_0", "9_0"]]

# #MCC_ab_new = TreeKnit.split_MCCs_topologically(MCC_ab, t_a_output, t_b_output)
# MCCs = TreeKnit.convert_MCC_list_to_set(3, ["a", "b", "c"], [MCC_ab, MCC_ac, MCC_bc])
# iMCC_split = TreeKnit.split_MCCs_topologically(MCCs, [t_a_output, t_b_output, t_c_output])
# new_MCCs_split_bc = TreeKnit.fix_consist!([MCC_ac, MCC_bc, MCC_ab_new], [t_a_output, t_b_output], i=2, merge=false)[2]
# new_MCCs_split_ac = TreeKnit.fix_consist!([MCC_ac, MCC_bc, MCC_ab_new], [t_a_output, t_b_output], i=1, merge=false)[1]
# print(TreeKnit.consistency_rate(MCC_ab_new, new_MCCs_split_ac, MCC_bc, [t_a_output, t_b_output, t_c_output]))
# print(TreeKnit.consistency_rate(MCC_ab_new, MCC_ac, new_MCCs_split_bc, [t_a_output, t_b_output, t_c_output]))

# new_MCCs_split_ab = TreeKnit.fix_consist!([MCC_ab_new, new_MCCs_split_ac , MCC_bc], [t_b_output, t_c_output], i=1, merge=false)[1]
# new_MCCs_split_bc_after = TreeKnit.fix_consist!([MCC_bc, new_MCCs_split_ab, new_MCCs_split_ac], [t_c_output, t_a_output], i=1, merge=false)[1]
# print(TreeKnit.consistency_rate(new_MCCs_split_ab, new_MCCs_split_ac, new_MCCs_split_bc_after, [t_a_output, t_b_output, t_c_output]))

# t_a = "((6_0:1.0,9_0:1.0)internal_1:1.0,((10_0:1.0,7_0:1.0)internal_10:1.0,((4_0:1.0,5_0:1.0)internal_3:1.0,(1_0:1.0,2_0:1.0,(8_0:1.0,3_0:1.0)RESOLVED_1:0.0)internal_6:1.0)internal_14:1.0)internal_21:1.0)internal_25:1.0;"
# t_b = "((1_0:1.0,2_0:1.0,(8_0:1.0,3_0:1.0)RESOLVED_2:0.0)internal_6:1.0,((6_0:1.0,9_0:1.0)internal_1:1.0,(10_0:1.0,(7_0:1.0,(4_0:1.0,5_0:1.0)internal_3:1.0)internal_8:1.0)RESOLVED_1:0.0)internal_21:1.0)internal_24:1.0;"
# t_c = "(((4_0:1.0,5_0:1.0)internal_3:1.0,(1_0:1.0,2_0:1.0,(8_0:1.0,3_0:1.0)internal_4:1.0)RESOLVED_1:0.0)internal_9:1.0,(10_0:1.0,(7_0:1.0,(6_0:1.0,9_0:1.0)internal_1:1.0)internal_18:1.0)internal_21:1.0)internal_23:1.0;"
# t_a_output = node2tree(TreeTools.parse_newick(t_a, node_data_type=TreeTools.MiscData), label = "a")
# t_b_output = node2tree(TreeTools.parse_newick(t_b, node_data_type=TreeTools.MiscData), label = "b")
# t_c_output = node2tree(TreeTools.parse_newick(t_c, node_data_type=TreeTools.MiscData), label = "c")
# MCC_ab = [["4_0", "5_0"], ["10_0", "6_0", "7_0", "9_0"], ["1_0", "2_0", "3_0", "8_0"]]
# MCC_ac = [["6_0", "9_0"], ["10_0", "1_0", "2_0", "3_0", "4_0", "5_0", "7_0", "8_0"]]
# MCC_bc = [["4_0", "5_0"], ["6_0", "9_0"], ["10_0", "1_0", "2_0", "3_0", "7_0", "8_0"]]
# MCCs = TreeKnit.convert_MCC_list_to_set(3, ["a", "b", "c"], [MCC_ab, MCC_ac, MCC_bc])
# iMCCs_fixed = TreeKnit.fix_consist!(MCCs, [t_a_output, t_b_output, t_c_output]; rounds=1, merge=false)
# print(iMCCs_fixed)

MCC_ab = [["5_0"], ["6_0", "10_0"], ["1_0", "2_0", "3_0", "4_0", "7_0", "8_0", "9_0"]]
MCC_ac = [["1_0", "4_0", "9_0"], ["10_0", "2_0", "3_0", "5_0", "6_0", "7_0", "8_0"]]
MCC_bc = [["5_0"], ["10_0", "6_0"], ["2_0", "7_0"], ["3_0", "8_0"], ["1_0", "4_0", "9_0"]]


nwk_a = "((8_0:1.0,3_0:1.0)internal_1:1.0,((9_0:1.0,4_0:1.0,1_0:1.0)RESOLVED_1:1.0,((10_0:1.0,6_0:1.0)RESOLVED_3:1.0,(5_0:1.0,(7_0:1.0,2_0:1.0)internal_4:1.0)internal_8:1.0)RESOLVED_2:1.0)internal_23:1.0)internal_26:1.0;"
nwk_b = "((5_0:1.0,(10_0:1.0,6_0:1.0)RESOLVED_4:1.0)internal_5:1.0,((8_0:1.0,3_0:1.0)RESOLVED_1:1.0,((7_0:1.0,2_0:1.0)RESOLVED_2:1.0,(9_0:1.0,4_0:1.0,1_0:1.0)RESOLVED_5:1.0)RESOLVED_3:1.0)internal_24:1.0)internal_27:1.0;"
nwk_c = "((9_0:1.0,4_0:1.0,1_0:1.0)internal_23:1.0,((8_0:1.0,3_0:1.0)internal_1:1.0,((10_0:1.0,6_0:1.0)RESOLVED_3:1.0,(5_0:1.0,(7_0:1.0,2_0:1.0)RESOLVED_1:1.0)RESOLVED_2:1.0)internal_11:1.0)internal_20:1.0)internal_26:1.0;" 
t_a_output = node2tree(TreeTools.parse_newick(nwk_a, node_data_type=TreeTools.MiscData), label = "a")
t_b_output = node2tree(TreeTools.parse_newick(nwk_b, node_data_type=TreeTools.MiscData), label = "b")
t_c_output = node2tree(TreeTools.parse_newick(nwk_c, node_data_type=TreeTools.MiscData), label = "c")
iMCCs = TreeKnit.convert_MCC_list_to_set(3, ["a", "b", "c"], [MCC_ab, MCC_ac, MCC_bc])
for i in range(1,10)
    new_sets = TreeKnit.fix_consist_sets!(iMCCs, [t_a_output, t_b_output, t_c_output])
    changed_ab, changed_ac = new_sets.mccs[Set(["b", "a"])], new_sets.mccs[Set(["c", "a"])] 
    option1 = [["5_0"], ["10_0", "6_0"], ["2_0", "7_0"], ["1_0", "3_0", "4_0", "8_0", "9_0"]], [["1_0", "4_0", "9_0"], ["10_0", "2_0", "3_0", "5_0", "6_0", "7_0", "8_0"]]
    option2 = [["5_0"], ["6_0", "10_0"], ["1_0", "2_0", "3_0", "4_0", "7_0", "8_0", "9_0"]], [["2_0", "7_0"], ["1_0", "4_0", "9_0"], ["10_0", "3_0", "5_0", "6_0", "8_0"]]
    @assert changed_ab== option1[1] && changed_ac == option1[2] || changed_ab== option2[1] && changed_ac == option2[2]
end

