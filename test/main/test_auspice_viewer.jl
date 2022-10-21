using Test
using TreeTools
using TreeKnit

println("### test auspice view of ARG ###")
nwk1 = "((1:1,2:1):1,((3:1,4:1):1,5:1,6:1):1)"
nwk2 = "((6:1,2:1):1,((3:1,4:1):1,5:1,1:1):1)"

t1 = node2tree(TreeTools.parse_newick(nwk1), label = "a")
t2 = node2tree(TreeTools.parse_newick(nwk2), label = "b")

MCCs = [["1"], ["2"], ["3", "4", "5", "6"]]

outdir = "test_auspice"
if !isdir(outdir)
	mkdir(outdir)
end
t1, t2, rS = TreeKnit.resolve_strict(t1, t2, MCCs )
TreeTools.ladderize!(t1)
TreeKnit.sort_polytomies!(t1, t2, MCCs )
TreeKnit.write_auspice_json(outdir * "/", t1, t2, MCCs )
TreeKnit.write_newick(outdir * "/" * t1.label *".nwk", t1)
TreeKnit.write_newick(outdir * "/" * t2.label*".nwk", t2)

## to test the auspice viewer now run the auspice_visulaization script from inside the test_auspice folder
## bash auspice_visulaization.sh a b 
