using Test
using TreeTools
using TreeKnit

MCCs_ref = [
	["A/NewYork/105/2002"],
	["A/NewYork/177/1999"],
	["A/NewYork/137/1999", "A/NewYork/138/1999"],
	["A/NewYork/52/2004", "A/NewYork/59/2003"],
	["A/NewYork/198/2003", "A/NewYork/199/2003", "A/NewYork/32/2003"],
	["A/NewYork/10/2004", "A/NewYork/11/2003", "A/NewYork/12/2003", "A/NewYork/13/2003", "A/NewYork/14/2003", "A/NewYork/15/2003", "A/NewYork/16/2003", "A/NewYork/17/2003", "A/NewYork/18/2003", "A/NewYork/19/2003", "A/NewYork/2/2003", "A/NewYork/21/2003", "A/NewYork/22/2003", "A/NewYork/23/2003", "A/NewYork/24/2003", "A/NewYork/25/2003", "A/NewYork/26/2003", "A/NewYork/27/2003", "A/NewYork/28/2003", "A/NewYork/29/2003", "A/NewYork/30/2003", "A/NewYork/31/2004", "A/NewYork/33/2004", "A/NewYork/34/2003", "A/NewYork/35/2003", "A/NewYork/36/2003", "A/NewYork/38/2003", "A/NewYork/39/2003", "A/NewYork/4/2003", "A/NewYork/40/2003", "A/NewYork/41/2003", "A/NewYork/42/2003", "A/NewYork/43/2003", "A/NewYork/44/2003", "A/NewYork/45/2003", "A/NewYork/46/2003", "A/NewYork/47/2003", "A/NewYork/48/2003", "A/NewYork/49/2003", "A/NewYork/5/2004", "A/NewYork/50/2003", "A/NewYork/51/2003", "A/NewYork/53/2003", "A/NewYork/54/2003", "A/NewYork/55/2003", "A/NewYork/56/2003", "A/NewYork/6/2004", "A/NewYork/60A/2003", "A/NewYork/61A/2003", "A/NewYork/62A/2003", "A/NewYork/63/2003", "A/NewYork/64/2003", "A/NewYork/65/2003", "A/NewYork/67/2003", "A/NewYork/69/2004", "A/NewYork/7/2003", "A/NewYork/70/2004", "A/NewYork/8/2003"],
]

t1 = read_tree("$(dirname(pathof(TreeKnit)))/../test/NYdata/tree_ha.nwk")
t2 = read_tree("$(dirname(pathof(TreeKnit)))/../test/NYdata/tree_na.nwk")

MCCs = computeMCCs(t1, t2)

outdir = "NY_treeknit_results"
# Setting up directories
mkpath(outdir)
rS = resolve!(t1, t2, MCCs)
TreeTools.ladderize!(t1)
TreeKnit.sort_polytomies!(t1, t2, MCCs)


arg, rlm, lm1, lm2 = SRG.arg_from_trees(t1, t2, MCCs)

TreeKnit.write_mccs(outdir * "/" * "MCCs.dat", MCCs)
out_nwk1, out_nwk2 = TreeKnit.make_output_tree_names("$(dirname(pathof(TreeKnit)))/../test/NYdata/tree_ha.nwk", "$(dirname(pathof(TreeKnit)))/../test/NYdata/tree_na.nwk")
TreeKnit.write_newick(outdir * "/" * out_nwk1, t1)
TreeKnit.write_newick(outdir * "/" * out_nwk2, t2)
TreeKnit.write(outdir * "/" * "arg.nwk", arg)
TreeKnit.write_rlm(outdir * "/" * "nodes.dat", rlm)

@testset "NY" begin
	@test MCCs[1:end-1] == MCCs_ref
end
if MCCs[1:end-1] != MCCs_ref
	@warn "Found different MCCs for the NewYork data. Could indicate a problem..."
end
