using Test
using TreeTools
using RecombTools, RecombTools.SplitGraph

println("##### main #####")

nwk1 = "((1:1,2:1):1,((3:1,4:1):1,5:1,6:1):1)"
nwk2 = "((6:1,2:1):1,((3:1,4:1):1,5:1,1:1):1)"

t1 = node2tree(TreeTools.parse_newick(nwk1))
t2 = node2tree(TreeTools.parse_newick(nwk2))

## Testing runopt
ct1 = deepcopy(t1)
ct2 = deepcopy(t2)
out = RecombTools.runopt(t1, t2; likelihood_sort=false, verbose=true)
##
trees = Dict(1=>deepcopy(t1), 2=>deepcopy(t2));
@testset "Pre-resolve" begin
	MCCs = RecombTools.computeMCCs!(
		trees, OptArgs(likelihood_sort=false);
		preresolve = true,
	)
	@test MCCs[1,2] == [["1"], ["2"], ["6"], ["3", "4", "5"]]
end

trees2 = Dict(1=>deepcopy(t1), 2=>deepcopy(t2));
solutions = [
	[["1"], ["2"], ["3", "4", "5", "6"]],
	[["2"], ["6"], ["1", "3", "4", "5"]],
	[["1"], ["6"], ["2", "3", "4", "5"]],
	[["1"], ["2", "6"], ["3", "4", "5"]],
	[["6"], ["1", "2"], ["3", "4", "5"]],
	[["1"], ["2"], ["6"], ["3", "4", "5"]],
]
@testset "Dyn-resolve" begin
	MCCs = RecombTools.computeMCCs!(
		trees2, OptArgs(likelihood_sort=false);
		preresolve=false,
	)
	@test in(MCCs[1,2], solutions)
end
