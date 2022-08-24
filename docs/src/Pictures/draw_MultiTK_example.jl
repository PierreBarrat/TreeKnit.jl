using TreeTools
using TreeKnit
using Plots

nwk1 = "((A:1,B:1):1,C:2)R:0;"
nwk2 = "(A:2,(B:1,C:1):1)R:0;"
nwk3 = "(A:2,B:2,C:2)R:0;"
t1 = node2tree(TreeTools.parse_newick(nwk1, node_data_type=TreeTools.MiscData), label = "a")
t2 = node2tree(TreeTools.parse_newick(nwk2, node_data_type=TreeTools.MiscData), label= "b")
t3 = node2tree(TreeTools.parse_newick(nwk3, node_data_type=TreeTools.MiscData), label= "c")


p1 = TreeKnit.draw_ARG(
        copy(convert(Tree{TreeTools.MiscData}, t1)), nothing,
        nothing, save=true, filename="tree1.png", font_size=18
)
p2 = TreeKnit.draw_ARG(
        copy(convert(Tree{TreeTools.MiscData}, t2)), nothing,
        nothing, save=true, filename="tree2.png", font_size=18
)
p3 = TreeKnit.draw_ARG(
        copy(convert(Tree{TreeTools.MiscData}, t3)), nothing,
        nothing, save=true, filename="tree3.png", font_size=18
)
full_p = plot(p1, p2, p3, layout = (1, 3), legend =true)
savefig(full_p, "MultiTK_example.png")