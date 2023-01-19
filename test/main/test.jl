using Test
using TreeTools
using TreeKnit, TreeKnit.SplitGraph, TreeKnit.MTK
using JSON3

nwk1 = "((1:1.0,2:1.0):1.0,((3:1.0,4:1.0):1.0,5:1.0,6:1.0):1.5);"
nwk2 = "((6:1,2:1):1,((3:1,4:1):1,5:1,1:1):1);"

t1 = parse_newick_string(nwk1)
label!(t1, "a")

t2 = parse_newick_string(nwk2)
label!(t2, "b")

solutions = [
	[["1"], ["2"], ["3", "4", "5", "6"]],
	[["2"], ["6"], ["1", "3", "4", "5"]],
	[["1"], ["6"], ["2", "3", "4", "5"]],
	[["1"], ["2", "6"], ["3", "4", "5"]],
	[["6"], ["1", "2"], ["3", "4", "5"]],
]

@testset "Dyn-resolve" begin
	MCCs = TreeKnit.computeMCCs(
		t1, t2, OptArgs(likelihood_sort=false, γ=3);
	)
	@test in(MCCs, solutions)
end
MCCs_MTK = MTK.compute_mcc_pairs!([t1, t2], TreeKnit.OptArgs(rounds=1))

@testset "write MCCs as JSON" begin
	MCCs_MTK = MTK.compute_mcc_pairs!([t1, t2], TreeKnit.OptArgs(rounds=1))
	json_output = TreeKnit.write_mccs(MCCs_MTK)
	mccs_json = JSON3.write(MCCs_MTK.mccs[Set([t1.label, t2.label])])
	@test json_output == "{\"MCC_dict\":{\"1\":{\"trees\":[\"a\",\"b\"],\"mccs\":"*mccs_json*"}}}"
end

MCCs = TreeKnit.MCC_set(2, ["a", "b"], [solutions[1]])
file = "$(dirname(pathof(TreeKnit)))/../test/main/auspice_a.json"
json_string = read(file, String)
json_output = TreeKnit.get_auspice_json(t1, [t1, t2], MCCs)
@testset "write auspice JSON" begin
	@test replace(json_string, " " => "", "\n" => "", "\t" => "") == json_output
end

@testset "tree naming" begin
	nwk = ["~/Documents/TreeKnit/TreeKnit.jl/examples/tree_h3n2_ha.nwk", "~/Documents/TreeKnit/TreeKnit.jl/examples/tree_h3n2_na.nwk"]
	fn, ext = TreeKnit.get_tree_names(nwk)
	@test fn == ["tree_h3n2_ha", "tree_h3n2_na"]
	@test TreeKnit.make_output_tree_names(fn, ext) == ["tree_h3n2_ha_resolved.nwk", "tree_h3n2_na_resolved.nwk"]

	nwk = ["~/Documents/TreeKnit/TreeKnit.jl/examples1/tree.nwk", "~/Documents/TreeKnit/TreeKnit.jl/examples2/tree.nwk"]
	fn, ext = TreeKnit.get_tree_names(nwk)
	@test fn == ["tree_examples1", "tree_examples2"]
	@test TreeKnit.make_output_tree_names(fn, ext) == ["tree_examples1_resolved.nwk", "tree_examples2_resolved.nwk"]
end
