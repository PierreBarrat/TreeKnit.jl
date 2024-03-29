using Test
using TreeKnit
using TreeTools

println("##### mcc_tools.jl #####")

t1 = parse_newick_string("((A,B)i1,(((C,D)i2,(E,(F1,F2)i3)i4)i5,G)i6)i7;")
t1_md = convert(Tree{TreeTools.MiscData}, t1)
MCCs = TreeKnit.sort([["A", "B", "E", "F1", "F2"], ["G"], ["C", "D"]], lt=TreeKnit.clt)
real_map = Dict{String, Union{Int, Nothing}}(
	"A" => 3,
	"B" => 3,
	"C" => 2,
	"D" => 2,
	"E" => 3,
	"F1" => 3,
	"F2" => 3,
	"G" => 1,
	"i1" => 3,
	"i2" => 2,
	"i3" => 3,
	"i4" => 3,
	"i5" => 3,
	"i6" => 3,
	"i7" => 3,
)

@testset "map_mccs base" begin
	@test map_mccs(t1, MCCs) == real_map

	leaf_map_1 = map_mccs(t1, MCCs; internals=false)
	for (k,v) in leaf_map_1
		@test real_map[k] == v
	end

	leaf_map_2 = map_mccs(MCCs)
	for (k,v) in leaf_map_2
		@test real_map[k] == v
	end

	t1_md_copy = copy(t1_md)
	map_mccs!(t1_md_copy, MCCs; internals=true)
	for n in nodes(t1_md_copy)
		@test n.data["mcc"] == real_map[n.label]
		@test haskey(n.data, "child_mccs")
	end

	t1_md_copy = copy(t1_md)
	map_mccs!(t1_md_copy, MCCs; internals=false)
	for n in leaves(t1_md_copy)
		@test n.data["mcc"] == real_map[n.label]
	end
	for n in internals(t1_md_copy)
		@test !haskey(n.data, "mcc")
		@test !haskey(n.data, "child_mccs")
	end

	println("-- Error correctly thrown below --")
	@test_throws ErrorException map_mccs!(t1, MCCs)
end



tree = parse_newick_string("((A,B)i1,((C,D)i2,(E,(F1,F2)i6,G)i5)i3)i7;")
tree_misc = convert(Tree{TreeTools.MiscData}, tree)

MCC = [["G"], ["A", "B"], ["C", "D"], ["E", "F1", "F2"]]
real_map = Dict{String, Union{Int, Nothing}}(
	"B" => 2, 
	"A" => 2, 
	"C" => 3, 
	"D" => 3, 
	"G" => 1, 
	"E" => 4, 
	"F1" => 4, 
	"F2" => 4, 
	"i1" => 2, 
	"i2" => 3, 
	"i3" => nothing, 
	"i5" => nothing, ##could be in MCC 1 or 4
	"i6" => 4, 
	"i7" => nothing
)

@testset "map_mccs with unknown internals" begin
    empty_data_tree_map = TreeKnit.map_mccs(tree, MCC)
	@test empty_data_tree_map== real_map 
	misc_data_tree_map = TreeKnit.map_mccs(tree_misc, MCC)
    @test misc_data_tree_map  == real_map
	TreeKnit.map_mccs!(tree_misc, MCC)
	for n in POT(tree_misc)
		@test real_map[n.label] == n.data["mcc"]
	end
end
