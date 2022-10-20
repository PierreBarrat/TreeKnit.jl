using Test
using TreeKnit
using TreeTools

println("##### mcc_tools.jl #####")

t1 = node2tree(parse_newick("((A,B)i1,(((C,D)i2,(E,(F1,F2)i3)i4)i5,G)i6)i7"))
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

	println("-- Error thrown below --")
	@test_throws ErrorException map_mccs!(t1, MCCs)
end



# Not implemented yet
# @testset "mark_shared_branches functions" begin
#     TreeKnit.mark_shared_branches!(constraint, tree3)
#     true_labels = Set(["A", "B", "C", "D", "E", "F1", "F2", "i6"])
#     for node in nodes(tree3)
#         if node.label in true_labels
#             @test node.data["shared_branch"] == true
#         else
#             @test node.data["shared_branch"] == false
#         end
#     end
# end
